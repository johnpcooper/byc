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

    Should be run after segemnting cells, then creating measurement rois
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
        cellrow = mdf.loc[index, :]
        base_filename = f'{cellrow.date}_byc_xy{str(cellrow.xy).zfill(2)}_cell{str(cellrow.cell_index).zfill(3)}'
        pattern = f'({base_filename})_(.+)_(stack.csv)'
        for dirpath, dirnames, filenames in os.walk(cellrow.compartment_dir):   
            measurement_dfs = []
            matches = []
            for filename in filenames:
                match = re.search(pattern, filename)
                if match:
                    # Define path to file and read in as dataframe
                    filepath = os.path.join(dirpath, filename)
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
                    newcols = [f'{col}_{channel}' for col in df.columns]
                    df.columns = newcols
                    measurement_dfs.append(df)

        if len(measurement_dfs) == 1:
            # Only one channel df found
            cell_df = measurement_dfs[0]        
        elif len(measurement_dfs) > 1:
            # Multiple channels, need to combine
            for df in measurement_dfs[1:]:
                cell_df = measurement_dfs[0].join(df)
        else:
            # No data found
            print(f'No data found for cell with pattern: {pattern}')
            print(f'Checked compartment dir: {cellrow.compartment_dir}')
            return None
        
        # Add annotation information to the cell trace df
        for col in mdf.columns:
            cell_df[col] = cellrow[col]
        # Add time information defined above
        cell_df.loc[:, 'hours'] = hours
        cell_df.loc[:, 'abs_hours'] = abs_hours
        cell_df.loc[:, 'hours_to_death'] = hours_to_death
        cell_df.loc[:, 'total_hours'] = (len(hours)*collection_interval)/60
            
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
                cell_df = self.make_cell_trace_df(mdf, i)
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




