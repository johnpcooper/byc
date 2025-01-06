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
        self.channel = kwargs.get('channel', 'yfp')
        # example gotten from one of byc.database.byc_database.master_index_dfs_dict.keys()
        example_compartment_name = '20200214_byc'
        if manual_select:
            # set the files that will be used for every cell analyzed
            self.master_index_df = self.select_master_index_df("Choose the .csv for the master index of this experiment")
        else:
            self.master_index_df = kwargs.get('mdf', None)

        # Read in data for each cell in the master index, add proper
        # time column and other annotations from the master index df
        self.cell_trace_dfs = self.make_cell_trace_dfs(self.master_index_df, channel=self.channel)
        self.compartment_name = kwargs.get('compartment_name', self.get_compartment_name())
        self.compartment_dir = self.get_compartment_dir()

    def get_compartment_name(self, **kwargs):
        mdf = kwargs.get('mdf', self.master_index_df)
        if type(mdf) == pd.core.series.Series:
            compartment_name = mdf.compartment_name
        else:
            compartment_name = mdf.compartment_name.iloc[0]
        return compartment_name

    def get_compartment_dir(self, **kwargs):
        mdf = kwargs.get('mdf', self.master_index_df)

        if type(mdf) == pd.core.series.Series:        
            compartment_dir = os.path.join(constants.byc_data_dir,
                                        mdf.compartment_reldir)
        else:
            compartment_dir = os.path.join(constants.byc_data_dir,
                                        mdf.compartment_reldir.iloc[0])
        return compartment_dir

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
        compartment_dir = os.path.join(constants.byc_data_dir, cellrow.compartment_reldir)
        for dirpath, dirnames, filenames in os.walk(compartment_dir):   
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
                    hours = ((df.frame)*collection_interval)/60
                    crop_roi_set_path = os.path.join(constants.byc_data_dir, cellrow.crop_roi_set_relpath)
                    if type(crop_roi_set_path) != str:
                        crop_roi_set_path = cellrow.crop_roi_set_path
                    print(f"Looking for crop ROIs at {crop_roi_set_path}")
                    rois = read_roi_zip(crop_roi_set_path)
                    first_position = rois[list(rois.keys())[0]]['position']
                    first_position_ind = first_position - 1
                    abs_hours = ((df.frame+first_position_ind)*collection_interval)/60
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

    # def make_cell_trace_dfs(self, master_index_df):
    #     """
    #     Use make_cell_df() to look up and aggregate
    #     fluorescent channel trace data for each cell
    #     recorded in the provied master index
        
    #     Return a list of these dataframes
    #     """
    #     mdf = master_index_df
    #     cell_dfs = []
    #     for i in mdf.index:
    #         cellrow = mdf.loc[i, :]
    #         if len(cellrow.shape) > 1:
    #             print(f'Multiple entries for cell at row {i}')
    #             pass
    #         else:
    #             cell_index = cellrow.cell_index
    #             cell_df = self.make_cell_trace_df(mdf, cell_index)
    #             cell_dfs.append(cell_df)
        
    #     return cell_dfs

    def make_cell_trace_dfs(self, mdf, **kwargs):
        """
        
        """
        if len(mdf.cell_index.unique()) == len(mdf.cell_index):
            pass
        else:
            print(f'There are multiple entries for the same cell in the master index')
                
            return None
        channel = kwargs.get('channel', mdf.channels_collected.iloc[0].split()[0])
        col_name = f'{channel}_df_path'
        mdf.loc[:, col_name] = np.nan
        compartmentname = mdf.compartment_name.iloc[0]
        date = compartmentname[0:8]
        if 'byc' in compartmentname:
            exptname = f'{date}_byc'
        elif 'fylm' in compartmentname:
            exptname = f'{date}_fylm'
        else:
            print(f'Compartment name needs to contain date followed by expt typ (fylm or byc)')
        kwargs = {'return_exptdir': True}
        exptdir, compdir = files.get_byc_compartmentdir(exptname, compartmentname, **kwargs)

        filenames = os.listdir(compdir)
        csvs = [f for f in filenames if f[-4:] == '.csv']
        measurement_csvs = [csv for csv in csvs if f'{channel}_stack' in csv]
        measurement_csv_paths = [os.path.join(compdir, name) for name in measurement_csvs]

        trace_dfs = []
        for cell_index in mdf.cell_index:
            # print(f'Finding measurements for cell {cell_index}')
            pattern = f'cell{str(cell_index).zfill(3)}'
            candidates = [re.search(pattern, path) for path in measurement_csv_paths]
            matches = [candidate for candidate in candidates if candidate != None]
            
            if len(matches) == 1:
                path = matches[0].string
                mdf.loc[cell_index, col_name] = path
                trace_df = pd.read_csv(path)
                trace_df.rename(columns={'Mean': f'Mean_{channel}'}, inplace=True)
                for col in mdf.columns:
                    if col not in list(trace_df.columns):
                        trace_df.loc[:, col] = np.nan
                        new_val = mdf.loc[cell_index, col]
                        trace_df.loc[:, col] = new_val
                        
                trace_dfs.append(trace_df)
                
            elif len(matches) > 1:
                print(f'Multiple measurement dfs for cell {cell_index} {channel} channel')
                for match in matches:
                    print(match.string)
            elif len(matches) == 0:
                print(f'No measurement dfs for cell {cell_index} {channel} channel in compartment dir {compdir}')

        return trace_dfs

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

def add_first_crop_frame(mdf):
    """
    For each entry in the mdf, read in that cell's crop roi set
    as a dataframe and mark 
    
    Return the annotated mdf
    """
    mdf.loc[:, 'first_crop_frame'] = np.nan
    # Just in case the master index df has a different
    # cell index order
    for idx in mdf.index:
        cell_index = mdf.loc[idx, 'cell_index']
        crop_roi_set_relpath = mdf[mdf.cell_index==cell_index].crop_roi_set_relpath.iloc[0]
        # print(f'Crop roi set relpath:\n{crop_roi_set_relpath}')
        crop_roi_set_path = os.path.join(constants.byc_data_dir, crop_roi_set_relpath)
        print(f'Looking for crop roi set path for cell {cell_index} at \n{crop_roi_set_path}')
        if os.path.exists(crop_roi_set_path):
            mdf.loc[idx, 'crop_roi_set_path'] = crop_roi_set_path
            crop_roi_df = files.read_rectangular_rois_as_df(crop_roi_set_path)
            first_crop_roi_pos = crop_roi_df.loc[0, 'position']
            first_crop_roi_frame = first_crop_roi_pos - 1
            
        else:
            first_crop_roi_frame = np.nan

        mdf.loc[idx, 'first_crop_frame'] = first_crop_roi_frame
        
    return mdf


def annotate_buds(mdf, buds_mdf,
                  return_bud_roi_df=False,
                  daughters=False,
                  observed_since_start_cutoff_frame=18):

    default_chase_frame = 150

    mdf.loc[:, 'bud_roi_set_path'] = np.nan
    mdf.loc[:, 'rls'] = np.nan
    mdf.loc[:, 'dist_from_sen'] = np.nan
    mdf.loc[:, 'age_at_chase'] = np.nan

    if len(mdf) != len(buds_mdf):
        print('Warning, mdf and buds_mdf not the same length')
        return None
    else:
        pass
    bud_dfs = []
    for idx in buds_mdf.index:
        # Make sure it's actually a relative path. It's often misannotated
        # in older data
        bud_roi_set_relpath = utilities.get_relpath(buds_mdf.bud_roi_set_relpath[idx])
        bud_roi_set_path = os.path.join(constants.byc_data_dir, bud_roi_set_relpath)
        cell_index = buds_mdf.cell_index[idx]
        print(f'Annotating bud appearance frames for cell {cell_index}')
        if 'abs_chase_frame' in mdf.columns:
            # If not here, this is probably just a bud ROIs master index
            # so just make up a number since it won't be relevant for these
            # data 
            abs_chase_frame = mdf[mdf.cell_index==cell_index].abs_chase_frame.iloc[0]
            print(f'Found absolute chase frame of {abs_chase_frame}')
        else:
            abs_chase_frame = default_chase_frame
            print(f'No abs_chase_frame column found, using default value of frame {abs_chase_frame}')
        # Need to define death offset because if death was observed, the last
        # annotated 'bud' frame is actually the last frame before death, not
        # a budding event
        if buds_mdf.loc[idx, 'end_event_type'] == 'death':
            death_offset = 1
            print(f'death observed for cell {cell_index}')
        elif buds_mdf.loc[idx, 'end_event_type'] == 'sen':
            death_offset = 0
            print(f'Death not observed for cell {cell_index}')
        elif buds_mdf.loc[idx, 'end_event_type'] == 'escape':
            # If the cell escapes, the last frame is still 
            # the last frame at which it was een alive, so 
            # the last frame is not an actual budding event
            death_offset = 1
        if os.path.exists(bud_roi_set_path):
            mdf.loc[idx, 'bud_roi_set_path'] = bud_roi_set_path
            print(f'Read in bud_roi_set at\n{bud_roi_set_path}')
            bud_roi_df = files.read_rectangular_rois_as_df(bud_roi_set_path)
            bud_roi_df.loc[:, 'frame'] = bud_roi_df.position - 1
            if daughters == False:
                age_at_chase = len(bud_roi_df.loc[bud_roi_df.frame<=abs_chase_frame, :])
                dist_from_sen = len(bud_roi_df.loc[bud_roi_df.frame>abs_chase_frame, :]) - death_offset
            else:
                # If we're annotating daughters "dist_from_sen", then dist_from_sen means
                # the age/dist_from_sen of the mother cells when the daughter cell was
                # born. In that case, mother dist_from_sen when daughter is born should
                # be calculate as number of buds appearing after first crop frame of the daughter
                # because that first crop frame will be when the daughter was born
                first_crop_frame = mdf.loc[idx, 'first_crop_frame']
                print(f'Annotating mother dist from sen when daughter was born')
                age_at_chase = len(bud_roi_df.loc[bud_roi_df.frame<=first_crop_frame, :])
                dist_from_sen = len(bud_roi_df.loc[bud_roi_df.frame>first_crop_frame, :]) - death_offset
            print(f'Gen. from senescence = {dist_from_sen}')
            rls = len(bud_roi_df) - death_offset
            # Annotate master index df with cell survival, first bud roi,
            # distance from senescence at start of chase, etc.
            mdf.loc[idx, 'bud_id'] = mdf.loc[idx, 'compartment_name'] +'-' + '-'.join(bud_roi_df.frame.apply(lambda x: str(x)))
            if bud_roi_df.frame.min() <= observed_since_start_cutoff_frame:
                mdf.loc[idx, 'observed_since_start'] = True
            else:
                mdf.loc[idx, 'observed_since_start'] = False            
            mdf.loc[idx, 'bud_roi_set_path'] = bud_roi_set_path
            mdf.loc[idx, 'rls'] = rls
            mdf.loc[idx, 'dist_from_sen'] = dist_from_sen
            print(f'Set value for dist_from_sen {dist_from_sen}')
            mdf.loc[idx, 'age_at_chase'] = age_at_chase
            first_bud_roi_pos = bud_roi_df.loc[0, 'position']
            mdf.loc[idx, 'first_bud_frame'] = first_bud_roi_pos - 1
            bud_roi_df.loc[:, 'cell_index'] = cell_index
            bud_dfs.append(bud_roi_df)
            
            print(f"Final dist from sen={mdf.loc[idx, 'dist_from_sen']}")
        else:
            print(f'No bud roi set found at:\n{bud_roi_set_path}')
    if return_bud_roi_df:
        print(f'Returning mdf and bud_mdf')
        return (mdf, pd.concat(bud_dfs, ignore_index=True))
        
    else:
        return mdf

def t0_normalize_trace_df(cell_trace_df, yvar='Mean_yfp', **kwargs):
    """
    Normalize <yvar> column in <cell_trace_df> to y value
    at t0 as define by <chase_frame> value found in 
    <cell_trace_df>. Uses background substracted yvar
    where background is minimum yvar value
    
    Return <cell_trace_df> with the new y_norm column
    """
    norm_col_name = kwargs.get('norm_col_name', None)
    delta_t = kwargs.get('delta_t', 10)
    tracedf = cell_trace_df
    tracedf.sort_values(by='frame', inplace=True)
    chase_frame = int(tracedf.chase_frame.iloc[0])
    tracedf.loc[:, f'{yvar}_bg_sub'] = np.nan
    tracedf.loc[:, f'{yvar}_bg_sub'] = tracedf[yvar] - tracedf[yvar].min()
    if chase_frame < len(tracedf):
        yt0 = tracedf.loc[chase_frame, f'{yvar}_bg_sub']
    else:
        print(f'Chase frame {chase_frame} outside of trace data of length {len(tracedf)}')
        yt0 = tracedf.loc[tracedf.index.min(), f'{yvar}_bg_sub']
    y_norm = tracedf[f'{yvar}_bg_sub']/yt0
    if norm_col_name is None:
        yvar_name = yvar[-3:]
        norm_col_name = f'{yvar_name}_norm'
    tracedf.loc[:, norm_col_name] = np.nan
    tracedf.loc[:, norm_col_name] = y_norm
    if 'frame_rel' not in list(tracedf.columns):
        try:
            tracedf.loc[:, 'frame_rel'] = tracedf.frame
        except:
            tracedf.loc[:, 'frame_rel'] = tracedf.frame
    tracedf.loc[:, 'hours'] = (tracedf.frame_rel*delta_t)/60
    tracedf.loc[:, 'hours_rel'] = tracedf.hours - (chase_frame*delta_t)/60
    
    return cell_trace_df

def create_and_annotate_mdf(exptname, compartmentname,
                            channels=None, daughters=False,
                            **kwargs):
    """
    Look in the compartmentdir found using <exptname> and
    <compartmentname>, create a master index using the
    crop_roi_df csvs found in the compartmentdir, then annotate
    that master index using bud rois found in the compartmentdir
    """
    # Allow user to pass an mdf that may have been manually
    # edited etc.
    mdf = kwargs.get('mdf', None)
    add_buds_to_mdf = kwargs.get('add_buds_to_mdf', True)
    savepath =kwargs.get('savepath', None)
    savemdf = kwargs.get('savemdf', True)
    if channels == None:
        print(f'Defaulting to channels bf, rfp, yfp')
        channels = ['bf', 'yfp', 'rfp']
    else:
        pass
    channels_str = ' '.join(channels)
    print(f'Analyzing channels {channels}')
    age_state = kwargs.get('age_state', 'old')
    chase_frame_dict = kwargs.get('chase_frame_dict', None)
    compartmentdir = files.get_byc_compartmentdir(exptname, compartmentname)
    print(f'Found compartment directory:\n{compartmentdir}')
    if mdf is None:
        try:
            mdf_type = 'crop_rois'
            file_pattern = constants.patterns.crop_roi_df_file
            mdf, savepath = files.mdf_from_file_pattern(compartmentdir, file_pattern, mdf_type=mdf_type)            
        except Exception as E:
            print(f'Could not generate a master index df using {mdf_type} csvs\n')
            mdf_type = 'bud_rois'
            file_pattern = constants.patterns.bud_roi_df_file
            mdf, savepath = files.mdf_from_file_pattern(compartmentdir, file_pattern, mdf_type=mdf_type)
            mdf.loc[:, 'compartment_name'] = compartmentname
            mdf = annotate_buds(mdf, mdf)
            return mdf
    else:
        print(f"Using user defined mdf and savepath")
    # Make sure that the relpaths set in the master index are truly
    
    # relative
    for col in mdf.columns:
        if 'rel' in col:
            mdf.loc[:, col] = [utilities.get_relpath(val) for val in mdf[col].values]
    mdf = add_first_crop_frame(mdf)
    # Annotate when the chases start
    mdf.loc[:, 'chase_frame'] = np.nan
    mdf.loc[:, 'abs_chase_frame'] = np.nan
    if chase_frame_dict is None:
        print('Please pass a chase_frame_dict with {<first_crop_frame>: <chase_frame>}')
        return None
    else:
        for first_frame, chase_frame in chase_frame_dict.items():
            abs_chase_frame = chase_frame + first_frame
            mdf.loc[mdf.first_crop_frame==first_frame, 'abs_chase_frame'] = abs_chase_frame
            mdf.loc[mdf.first_crop_frame==first_frame, 'chase_frame'] = chase_frame

    mdf.loc[:, 'compartment_name'] = compartmentname
    mdf.loc[:, 'age_state'] = age_state

    mdf.loc[:, 'channels_collected'] = channels_str

    # Add paths to cell tracking ROIs
    mdf = files.path_annotate_master_index_df(mdf, channels=channels)
    if len(mdf) > 0:
        print(f'Successfully path annotated master index')
    # Create buds master index using bud roi dfs found
    # in compartment directory, use them to annotate the
    # master index created above
    print(f'Attempting to annotate bud information')
    mdf_type = 'bud_rois'
    file_pattern = constants.patterns.bud_roi_df_file
    compartmentdir = files.get_byc_compartmentdir(exptname, compartmentname)
    print(f'Found compartment directory:\n{compartmentdir}')
    bud_mdf = files.mdf_from_file_pattern(compartmentdir, file_pattern, mdf_type=mdf_type, return_savepath=False)
    for col in mdf.columns:
        if 'rel' in col:
            mdf.loc[:, col] = [utilities.get_relpath(val) for val in mdf[col].values]
    if add_buds_to_mdf:
        print(f'{len(mdf)} cells in mdf')
        print(f'{len(bud_mdf)} cells in buds_mdf')
        mdf = annotate_buds(mdf, bud_mdf, daughters=daughters)
    else:
        print(f'Not annotating budding events on master index')
    # print(mdf.dist_from_sen.unique())
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

def make_fits_table(fits_df, **kwargs):
    """
    """
    # Drop cells that are clearly bad data based on above plots
    drops = kwargs.get('drops_cell_indices', [])
    agg_idx_temp = ['age_at_chase',
               'rls',
               'dist_from_sen',
               'first_bud_frame',
               'cell_index']
    # agg_idx_temp = list(fits_df.columns)
    agg_idx = kwargs.get('agg_idx', agg_idx_temp)
    fits_df = fits_df[~(fits_df.cell_index.isin(drops))].reset_index()
    
    fits_table = pd.pivot_table(index=agg_idx, data=fits_df).reset_index()
    
    return fits_table

def merge_dfs(df1, df2, **kwargs):
    """
    Take two dataframes, <df1> usually a master_index_df,
    <df2> usually a fits table, and add all data from df1
    to df2 that's not already in df2. 

    Unique cells identified by 'cell_index' and 'dist_from_sen'.
    In current byc data sets, a new 'cell_index' only means
    trace index. Cell 20 and 21 could be same cell chase at
    different timepoints. They can be identified as the same
    cell using their set of bud appearance timepoints

    In older datasets (e.g. 20210430), I was creating a young and
    old compartmentdir to distinguish different chase timepoints
    for the same cell which created a mess

    Return df2
    """
    idx = kwargs.get('idx', ['cell_index', 'dist_from_sen'])
    df1.set_index(idx, inplace=True, drop=False)
    df2.set_index(idx, inplace=True, drop=False)
    notfounds = []
    for cell_dist_dex in df2.index:
        if cell_dist_dex in list(df1.index):
            for col in df1.columns:
                if col not in df2.columns:
                    df2.loc[:, col] = np.nan
                value = df1.loc[cell_dist_dex, col]
                df2.loc[cell_dist_dex, col] = value
        else:
            print(f'Cell with index, dist_from_sen {cell_dist_dex} not in mdf')
            notfounds.append(1)
    print(f'Did not find {len(notfounds)} of {len(df2)} in <df1>')
    df2.index = range(0, len(df2))
    df1.reset_index(inplace=True, drop=True)
    return df2