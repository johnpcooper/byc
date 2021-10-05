import os
import re

import tkinter as tk
import tkinter.filedialog as tkdia

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from read_roi import read_roi_zip

from byc import constants, utilities, database, files

# Recently updated bycDataSet so that you instantiate
# it by giving it a key through which it can find
# the master index in byc.database.byc_database.master_index_dfs_dict
class bycDataSet(object):
    """
    Upon instantiation, ask the user to choose a master index for the experiment
    they want to create trace dataframes for each cell in the master index.

    Should be run after segmenting cells, then creating measurement rois
    using imagejpc/utilities/save_cell_roi_set, and then creating and saving 
    those measurements using imagejpc/utilities/measure_rois
    """       
    def __init__(self, manual_select=False, **kwargs):

        # example gotten from one of byc.database.byc_database.master_index_dfs_dict.keys()
        example_compartment_name = '20200214_byc'
        self.compartment_name = kwargs.get('compartment_name', None)
        self.master_index_df = kwargs.get('mdf', None)
        if manual_select:
            # set the files that will be used for every cell analyzed
            self.master_index_df = self.select_master_index_df("Choose the .csv for the master index of this experiment")
        else:
            # Use compartment name to fetch a master index df
            # using database.byc_database.
            # If no kwarg for compartment_name, look for a
            # master index df object in kwargs as 'mdf'

            if self.compartment_name != None:
                self.master_index_df = database.byc_database.master_index_dfs_dict[self.compartment_name]
                self.master_index_df = self.clean_master_index_df(self.master_index_df)

        # Read in data for each cell in the master index, add proper
        # time column and other annotations from the master index df
        self.cell_trace_dfs = self.make_cell_trace_dfs(self.master_index_df)

    def set_fp(self, prompt):
        """
        Return the path to a file of interest. Call with the prompt you would like to 
        display in the dialog box.
        """
        # create the dialog box and set the fn
        root = tk.Tk()
        fp = tkdia.askopenfilename(parent=root, title=prompt)
        root.destroy() # very important to destroy the root object, otherwise python 
        # just keeps running in a dialog box

        return fp # return the path to the file you just selected

    def clean_master_index_df(self, master_index_df):

        if 'sub_coord' in master_index_df.columns:
            master_index_df.rename(columns={'sub_coord': 'cell_index'},
                                   inplace=True)
        if 'path' in master_index_df.columns:
            master_index_df.rename(columns={'path': 'compartment_dir'},
                                   inplace=True)
            # add compartment_reldir here
            abspaths = master_index_df.compartment_dir
            relpaths = [utilities.get_relpath(abspath) for abspath in abspaths]
            master_index_df.loc[:, 'compartment_reldir'] = relpaths

        return master_index_df

    def select_master_index_df(self, prompt):
        """
        Return a dataframe read from the master cells index .csv
        """
        # define the path to the index .csv
        master_cells_fp = self.set_fp(prompt)
        # define the filename for the master expt index
        master_index_df = pd.read_csv(master_cells_fp)
        master_index_df = self.clean_master_index_df(master_index_df)

        return master_index_df

    def make_cell_trace_df(self, master_index_df, index, collection_interval=10):
        """
        Find all .csv fluorescent trace measurement files
        for the cell at index <index> in the master_index_df,
        add labels found in the master index for that cell,
        and return the completed dataframe
        """
        mdf = master_index_df
        cellrow = mdf.set_index('cell_index').loc[index, :]
        cellrow['cell_index'] = index
        base_filename = f'{cellrow.date}_byc_xy{str(cellrow.xy).zfill(2)}_cell{str(cellrow.cell_index).zfill(3)}'
        pattern = f'({base_filename})_(.+)_(stack.csv)'
        measurement_dfs = []
        for dirpath, dirnames, filenames in os.walk(cellrow.compartment_dir):   
            matches = []
            for filename in filenames:
                match = re.search(pattern, filename)
                # Define path to file and read in as dataframe
                filepath = os.path.join(dirpath, filename)
                if match:
                    print(f'Found data for cell {index} at {filepath}')
                    df = pd.read_csv(filepath)
                    # Set time data. Absolute hours is hours
                    # since the experiment started. Hours is
                    # hours since beginning of data collection
                    #
                    # Note that I don't add these columns to the 
                    # df until I'm down finding new channel .csvs
                    # because I don't want the time columns to get
                    # labelled with '_<channel_name>'
                    hours = ((df.Slice-1)*collection_interval)/60
                    print(f"Looking for crop ROIs at {cellrow.crop_roi_set_path}")
                    rois = read_roi_zip(cellrow.crop_roi_set_path)
                    first_position = rois[list(rois.keys())[0]]['position']
                    first_position_ind = first_position - 1
                    abs_hours = ((df.Slice-1+first_position_ind)*collection_interval)/60
                    # Add hours to death column if end event for crop
                    # ROIs was 'death'. Otherwise set column to np.nan
                    last_position = rois[list(rois.keys())[-1]]['position']
                    last_position_ind = last_position - 1
                    if cellrow.end_event_type == 'death':
                        hours_to_death = hours[hours.index[-1]] - hours
                    else:
                        hours_to_death = np.full(len(hours), np.nan)
                    matches.append(match)
                    # Name the measurement data with the fluorescent channel.
                    # Add identifying information from master index
                    channel = match.group(2)
                    print(f'Adding data from {channel} channel')
                    newcols = [f'{col}_{channel}' for col in df.columns]
                    df.columns = newcols
                    if pd.DataFrame(df).empty:
                        print(f'Empty dataframe created from data found at \n{filepath}')
                    print(f'Length of channel df {len(df)}')
                    measurement_dfs.append(df)
                    print(f'Length of channel dfs now {len(measurement_dfs)}')
        n_meas_dfs = len(measurement_dfs)
        print(f'Found {n_meas_dfs} channel dfs')
        if n_meas_dfs == 1:
            # Only one channel df found
            cell_df = measurement_dfs[0]
        elif n_meas_dfs == 2:
            # Two channels, need to combine
            cell_df = measurement_dfs[0].join(measurement_dfs[1])
        elif n_meas_dfs == 3:
            # Three channels, need to combine. Haven't found
            # a good way to merge more than two dataframes :(
            cell_df = measurement_dfs[0].join(measurement_dfs[1])
            cell_df = cell_df.join(measurement_dfs[2])
        elif len(measurement_dfs) < 1:
            # No data found
            print(f'No data found for cell with pattern: {pattern}')
            print(f'Base filename: {base_filename}')
            print(f'Checked compartment dir: {cellrow.compartment_dir}')
            return None
        
        # Add annotation information to the cell trace df
        for col in mdf.columns:
            print(f'Adding {col} from master index')
            cell_df.loc[:, col] = cellrow[col]
        # Add time information defined above
        cell_df.loc[:, 'hours'] = hours
        cell_df.loc[:, 'frame_index'] = (hours*60)/collection_interval
        cell_df.loc[:, 'abs_hours'] = abs_hours
        cell_df.loc[:, 'hours_to_death'] = hours_to_death
        cell_df.loc[:, 'total_hours'] = (len(hours)*collection_interval)/60

        if pd.DataFrame(cell_df).empty==True:
            print(f'Warning, empty dataframe from cell with base filename {base_filename}')
        print(f'Created cell_df with length {len(cell_df)}')
        return cell_df

    def make_cell_trace_dfs(self, master_index_df):
        """
        Use make_cell_df() to look up and aggregate
        fluorescent channel trace data for each cell
        recorded in the provied master index
        
        Return a list of these dataframes
        """
        mdf = master_index_df
        cell_dfs = []
        for i in mdf.index:
            cellrow = mdf.loc[i, :]
            if len(cellrow.shape) > 1:
                print(f'Multiple entries for cell at row {i}')
                pass
            else:
                cell_index = cellrow.cell_index
                cell_df = self.make_cell_trace_df(mdf, cell_index)
                cell_dfs.append(cell_df)
        
        return cell_dfs

def set_file_paths(prompt):
    # create the dialog box and set the fn
    root = tk.Tk()
    fps = tkdia.askopenfilenames(parent=root, title=prompt)
    root.destroy() # very important to destroy the root object, otherwise python 
    # just keeps running in a dialog box

    return fps # return the path to the file you just selected

def get_dfs_list():
    fps = set_file_paths("Choose the .csv files for this expt condition")
    
    cell_trace_dfs_list = []
    for i in range(0, len(fps)):
        df = pd.read_csv(fps[i])
        cell_trace_dfs_list.append(df)
        
    return cell_trace_dfs_list

def get_bud_hours(celldf, reference='death'):
    """
    Return a numpy array containing timepoints
    of bud appearances

    Accepted references are:
    death (correspends to 'hours_to_death' column of celldf)
    expt_start (corresponds to 'abs_hours' column)
    cell_obs_start (corresponds to 'hours' column)
    """
    try:
        collection_interval = celldf.collection_interval.unique()[0]
    except:
        collection_interval = 10
        print(f'Warning! Defaulting to collection interval of {collection_interval} minutes')

    # Create empty columns
    bud_rois_inds = np.full(len(celldf), np.nan)
    bud_abs_hours = np.full(len(celldf), np.nan)
    bud_hours =  np.full(len(celldf), np.nan)
    bud_hours_to_death = np.full(len(celldf), np.nan)

    if True in celldf.bud_roi_set_path.values:
        print(f'No path found in celldf with xy {celldf.xy.unique()[0]} and cell_index{celldf.cell_index.unique()[0]}')
        return None
    else:
        bud_rois_path = celldf.bud_roi_set_path.iloc[0]
        print(bud_rois_path)
        bud_rois_inds = files.read_roi_position_indices(bud_rois_path)
        bud_abs_hours = (bud_rois_inds*collection_interval)/60
        bud_hours =  bud_abs_hours - celldf.abs_hours.min()
        bud_hours_to_death = celldf.hours.max() - bud_hours

    vals = [bud_rois_inds,
            bud_abs_hours,
            bud_hours,
            bud_hours_to_death]

    keys = ['bud_rois_inds',
            'bud_abs_hours',
            'bud_hours',
            'bud_hours_to_death']
    
    if reference == 'death':
        return bud_hours_to_death
    elif reference == 'experiment_start':
        return bud_abs_hours
    elif reference == 'cell_obs_start':
        return bud_hours
    else:
        print(f"Bad reference passed: {reference}")
        print("Returning bud roi frame indices")
        return bud_roi_inds

def annotate_absolute_time(cell_trace_df, mdf, time_delta_mins=10):
    """
    Look up the path to the crop rois .zip  in <mdf> using the cell_index
    found in <cell_trace_df>. Use that crop roi to find the first position
    of the crop roi and use that first position to set absolute experiment
    time column for the <cell_trace_df>

    Return the absolute time annotated <cell_trace_df>
    """


def annotate_buds(mdf, buds_mdf, abs_chase_frame, return_bud_roi_df=False):

    mdf.loc[:, 'bud_roi_set_path'] = np.nan
    mdf.loc[:, 'rls'] = np.nan
    mdf.loc[:, 'dist_from_sen'] = np.nan
    mdf.loc[:, 'age_at_chase'] = np.nan
    mdf.loc[:, 'first_crop_frame'] = np.nan
    bud_dfs = []
    for idx in buds_mdf.index:
        bud_roi_set_path = buds_mdf.bud_roi_set_path[idx]
        cell_index = buds_mdf.cell_index[idx]
        # Just in case the master index df has a different
        # cell index order
        crop_roi_set_relpath = mdf[mdf.cell_index==cell_index].crop_roi_set_relpath.iloc[0]
        crop_roi_set_path = os.path.join(constants.byc_data_dir, crop_roi_set_relpath)
        # Need to define death offset because if death was observed, the last
        # annotated 'bud' frame is actually the last frame before death, not
        # a budding event
        if buds_mdf.loc[idx, 'end_event_type'] == 'death':
            death_offset = 1
            print(f'death observed for cell {idx}')
        elif buds_mdf.loc[idx, 'end_event_type'] == 'sen':
            death_offset = 0
            print(f'Death not observed for cell {idx}')
        elif buds_mdf.loc[idx, 'end_event_type'] == 'escape':
            death_offset = 0
        if os.path.exists(bud_roi_set_path) and os.path.exists(crop_roi_set_path):
            mdf.loc[idx, 'bud_roi_set_path'] = bud_roi_set_path
            bud_roi_df = files.read_rectangular_rois_as_df(bud_roi_set_path)
            crop_roi_df = files.read_rectangular_rois_as_df(crop_roi_set_path)
            bud_roi_df.loc[:, 'frame'] = bud_roi_df.position - 1
            age_at_chase = len(bud_roi_df.loc[bud_roi_df.frame<=abs_chase_frame, :])
            dist_from_sen = len(bud_roi_df.loc[bud_roi_df.frame>abs_chase_frame, :]) - death_offset
            rls = len(bud_roi_df) - death_offset
            # Annotate master index df with cell survival, first bud roi,
            # distance from senescence at start of chase, etc.
            mdf.loc[idx, 'bud_roi_set_path'] = bud_roi_set_path
            mdf.loc[idx, 'rls'] = rls
            mdf.loc[idx, 'dist_from_sen'] = dist_from_sen
            mdf.loc[idx, 'age_at_chase'] = age_at_chase
            first_bud_roi_pos = bud_roi_df.loc[0, 'position']
            first_crop_roi_pos = crop_roi_df.loc[0, 'position']
            mdf.loc[idx, 'first_bud_frame'] = first_bud_roi_pos - 1
            mdf.loc[idx, 'first_crop_frame'] = first_crop_roi_pos - 1
            bud_roi_df.loc[:, 'cell_index'] = cell_index
            bud_dfs.append(bud_roi_df)
        else:
            print(f'No bud roi set found at:\n{bud_roi_set_path}')
            print(f'No crop roi set found at:\n{crop_roi_set_path}')
    if return_bud_roi_df:
        return (mdf, pd.concat(bud_dfs, ignore_index=True))
    else:
        return mdf

def t0_normalize_trace_df(cell_trace_df, yvar='Mean_yfp'):
    """
    Normalize <yvar> column in <cell_trace_df> to y value
    at t0 as define by <chase_frame> value found in 
    <cell_trace_df>. Uses background substracted yvar
    where background is minimum yvar value
    
    Return <cell_trace_df> with the new y_norm column
    """
    tracedf = cell_trace_df
    chase_frame = tracedf.chase_frame.iloc[0]
    tracedf.loc[:, f'{yvar}_bg_sub'] = np.nan
    tracedf.loc[:, f'{yvar}_bg_sub'] = tracedf[yvar] - tracedf[yvar].min()
    yt0 = tracedf.loc[chase_frame, f'{yvar}_bg_sub']
    y_norm = tracedf[f'{yvar}_bg_sub']/yt0
    yvar_name = yvar[-3:]
    norm_col_name = f'{yvar_name}_norm'
    tracedf.loc[:, norm_col_name] = np.nan
    tracedf.loc[:, norm_col_name] = y_norm
    
    return cell_trace_df

def create_and_annotate_mdf(exptname, compartmentname,
                            chase_frame, chase_roi_start_frame,
                            **kwargs):
    """
    Look in the compartmentdir found using <exptname> and
    <compartmentname>, create a master index using the
    crop_roi_df csvs found in the compartmentdir, then annotate
    that master index using bud rois found in the compartmentdir
    """
    savemdf = kwargs.get('savemdf', True)
    channels = kwargs.get('channels_collected', 'bf yfp')
    age_state = kwargs.get('age_state', 'old')
    abs_chase_frame = chase_frame + chase_roi_start_frame
    mdf_type = 'crop_rois'
    file_pattern = constants.patterns.crop_roi_df_file
    compartmentdir = files.get_byc_compartmentdir(exptname, compartmentname)
    print(f'Found compartment directory:\n{compartmentdir}')
    mdf, savepath = files.mdf_from_file_pattern(compartmentdir, file_pattern, mdf_type=mdf_type)
    mdf.loc[:, 'chase_frame'] = chase_frame
    mdf.loc[:, 'abs_chase_frame'] = abs_chase_frame
    mdf.loc[:, 'compartment_name'] = compartmentname
    mdf.loc[:, 'channels_collected'] = channels
    mdf.loc[:, 'age_state'] = age_state
    # Add paths to cell tracking ROIs
    mdf = files.path_annotate_master_index_df(mdf)
    # Create buds master index using bud roi dfs found
    # in compartment directory, use them to annotate the
    # master index created above
    mdf_type = 'bud_rois'
    file_pattern = constants.patterns.bud_roi_df_file
    compartmentdir = files.get_byc_compartmentdir(exptname, compartmentname)
    print(f'Found compartment directory:\n{compartmentdir}')
    bud_mdf = files.mdf_from_file_pattern(compartmentdir, file_pattern, mdf_type=mdf_type, return_savepath=False)

    mdf = annotate_buds(mdf, bud_mdf, abs_chase_frame)
    if savemdf:
        mdf.to_csv(savepath)
        print(f'Saved annotated master index df at:\n{savepath}')
    return mdf

def filter_low_dynamic_range_cells(trace_dfs, threshold, **kwargs):
    yvar = kwargs.get('yvar', 'Mean_yfp')
    filtered_trace_dfs = []
    for tracedf in trace_dfs:
        chase_frame = tracedf.chase_frame.unique()[0]
        ymin = tracedf[yvar].min()
        yt0 = tracedf.loc[chase_frame, yvar]

        delta_y = yt0 - ymin
        print(f'Range from {yt0} to {ymin}')
        print(delta_y)
        if delta_y < threshold:
            print('Cell thrown out')
        else:
            filtered_trace_dfs.append(tracedf)

    print(f'Kept {len(filtered_trace_dfs)} of {len(trace_dfs)}')
    return filtered_trace_dfs