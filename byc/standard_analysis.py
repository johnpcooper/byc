import os

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.filedialog as tkdia

from byc import constants, utilities, database

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
        self.compartment_name = kwargs.get('compartment_name', example_compartment_name)
        
        if manual_select:
            # set the files that will be used for every cell analyzed
            self.master_index_df = self.select_master_index_df("Choose the .csv for the master index of this experiment")
        elif self.compartment_name != None:
            self.master_index_df = database.byc_database.master_index_dfs_dict[self.compartment_name]
            self.master_index_df = self.clean_master_index_df(self.master_index_df)
        # Get fluorescent channels collected
        if 'channels_collected' in self.master_index_df.columns:
            self.channel_names = self.master_index_df.channels_collected.iloc[0].split()
            self.fluor_channel_names = [channel for channel in self.channel_names if channel != 'bf']
        else:
            # default to dsred and yfp if channels not
            # recorded in master index
            self.fluor_channel_names = constants.default_channels
            print(f'WARNING: defaulting to fluor channels: {self.fluor_channel_names}')
            print("Make sure to set 'channels_collected' column in master index")
        # Get image collection interval in minutes. Default to 10
        if 'collection_interval' in self.master_index_df.columns:
            self.collection_interval = self.master_index_df.collection_interval.iloc[0]
        else:
            print(f'WARNING: defaulting to 10 minute collection interval')
            print('Can be changed by adding collection_interval column in master_index')
            self.collection_interval = 10
        # Split up the master index df into a list of dfs, each
        # one being the cell's slice from the master index
        self.cells_dfs = self.set_cells_dfs(self.master_index_df)
        # set lists containing senescent slice for each cell
        self.senescent_slices = self.set_senenescent_slices(self.cells_dfs)
        # set lists containing dataframes with measurements for each cell in cells_dfs
        self.cell_trace_dfs_dict = self.get_raw_cell_trace_dfs_dict(self.cells_dfs, self.senescent_slices)
        self.cell_trace_dfs_dict = self.name_fluor_columns(self.cell_trace_dfs_dict)
        self.cell_trace_dfs = self.get_processed_cell_trace_dfs()
        self.cell_trace_paths = self.get_cell_trace_paths()

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


    def set_cells_dfs(self, master_index_df):
        """
        Return a list of sub-dataframes of the master_index_df. The sub-dataframes
        contain data for only one cell according to its cell_index
        """
        # create the cells_dfs list
        cells_dfs = []

        # add to the cells_dfs a dataframe at index i for every
        # unique value in the master_index_df['cell_index'] or 
        # 'sub_coord' if that's in the master index
        if 'sub_coord' in master_index_df.columns:
            cell_index_colname = 'sub_coord'
        elif 'cell_index' in master_index_df.columns:
            cell_index_colname = 'cell_index'
        else:
            print("Couldn't find cell_index or sub_coord column in master_index_df")
        for i in master_index_df[cell_index_colname].unique():
            # set the logic mask to a sub-dataframe of master_index_df containing
            # only values whose 'sub_coord' value is value
            print(f"Slicing master index df for cell {i}")
            logic_mask = (master_index_df[cell_index_colname] == i)
            cells_dfs.append(master_index_df[logic_mask])
            print("cells_dfs is now %d elements long" % len(cells_dfs))

        return cells_dfs
    
    def set_senenescent_slices(self, cells_dfs):        
        """
        Return a list of number indicating where along the expt each cell went
        senescent (ended its last division). If senescence wasn't observed
        for the cell, then False will be added to the list.
        """        
        senescent_slices = []
        
        for cell_index in range(0, len(cells_dfs)):
            
            cell_df = cells_dfs[cell_index]
            
            sen_value = cell_df['sen_start'][cell_df.index.min()]
            first_frame = cell_df['Pos'][cell_df.index.min()]
            last_frame = cell_df['Pos'][cell_df.index.max()]
            print("Checking senescence data for cell %s" % cell_index)
            print("sen_value is %s" % sen_value)
            print("first_frame is %s" % first_frame)
            print("last_frame is %s" % last_frame)
            
            if sen_value == "FALSE":
                
                adj_sen_value = "FALSE"
            
            elif sen_value != "FALSE":
            
                sen_distance_from_start = int(sen_value) - int(first_frame)
                print("sen_distance_from_start is %s" % sen_distance_from_start)

                adj_sen_value = sen_distance_from_start              
            
            senescent_slices.append(adj_sen_value)
            
        return senescent_slices

    def get_raw_cell_trace_dfs_dict(self, cells_dfs, senescent_slices):
        """
        Return a dictionary of cell channel trace dfs.
        The keys of the dict are self.flour_channel_names (set in __init__)
        and values are lists of individual cell trace dataframes for
        that channel_name
        """
        channel_names = self.fluor_channel_names
        channel_dfs_lists = [[] for name in channel_names]
        cell_traces_dict = dict(zip(channel_names, channel_dfs_lists))

        for cell_index in range(0, len(cells_dfs)):

            cell_df = cells_dfs[cell_index]
            # Accounting
            if 'sub_coord' in cell_df.columns:
                cell_index_colname = 'sub_coord'
            elif 'cell_index' in cell_df.columns:
                cell_index_colname = 'cell_index'
            else:
                print("Couldn't find cell_index or sub_coord column in master_index_df")

            if 'path' in cell_df.columns:
                path = os.path.abspath(cell_df.path.iloc[0])
            elif 'compartment_reldir' in cell_df.columns:
                path = os.path.join(constants.byc_data_dir,
                                    cell_df.compartment_reldir.iloc[0])
            else:
                print(f"""
                    No comparment directory data in cell {cell_index}'s master index

                    Either add a 'path' column containing the full path to the 
                    master index or make sure you are creating all ROIs for each
                    cell using the imagejpc plugin, 'save_cell_roi_set.py' which
                    constructs master index with proper formatting
                    """)            
            expt_date = cell_df['date'][cell_df.index.min()]
            expt_type = cell_df['expt_type'][cell_df.index.min()]
            xy = cell_df['xy'][cell_df.index.min()]
            cell_index = cell_df[cell_index_colname][cell_df.index.min()]
            
            stack_title = str("%s_%s_xy%s_cell%s" % (expt_date, expt_type, str(xy).zfill(2), str(cell_index).zfill(3)))

            # read the .csv containing measurements
            cell_trace_channel_paths = [os.path.join(path, f'{stack_title}_{name}_stack.csv') for name in channel_names]
            path_existances = [os.path.exists(path) for path in cell_trace_channel_paths]
            if False in path_existances:
                print("Old cell index formatting. Looking for single digit cell index and coordinate")
                stack_title = str("%s_%s_xy%s_cell%s" % (expt_date, expt_type, str(xy), str(cell_index)))
                # read the .csv containing measurements
                cell_trace_channel_paths = [os.path.join(path, f'{stack_title}_{name}_stack.csv') for name in channel_names]
            
            # Read in the dataframe for each channel collected and measured
            cell_trace_channel_dfs = []
            for path in cell_trace_channel_paths:
                cell_trace_channel_dfs.append(pd.read_csv(path))
            # create the time to senescence column
            if not (senescent_slices[cell_index] == "FALSE"):               
                for channel_df in cell_trace_channel_dfs:
                    channel_df['slices_to_senescence'] = channel_df['Slice'] - int(senescent_slices[cell_index]) 
                    channel_df['cell_index'] = cell_index
            else:
                print("Senescence not observed for cell %s" % cell_index)
            
            cell_trace_dict = dict(zip(channel_names, cell_trace_channel_dfs))
            # Add cell traces to dictionary that gets returned at
            # end of function
            for channel_name, channel_df in cell_trace_dict.items():
                cell_traces_dict[channel_name].append(channel_df)

        return cell_traces_dict

    def correct_time(self, dataframe):        
        """
        Return a DataFrame of the DataFrame passed to correct_time() with
        columns added for minutes and hours according to the collection_interval argument
        """
        collection_interval = self.collection_interval
        if "slices_to_senescence" in dataframe.columns:        
            dataframe.loc[:, 'minutes_to_senescence'] = dataframe['slices_to_senescence'] * collection_interval
            dataframe.loc[:, 'hours_to_senescence'] = dataframe['minutes_to_senescence'] / 60
        else:
            pass
        # Slice really refers to frame in time. It's the name
        # given by by imagej when making fluorescence measurements
        dataframe.loc[:, 'minutes'] = (dataframe['Slice']-1) * collection_interval
        dataframe.loc[:, 'hours'] = dataframe['minutes'] / 60
        
        return dataframe

    def cell_diameter(self, dataframe):
        """
        Return a DataFrame with added column for cell diameter in um
        """        
        dataframe.loc[:, 'cell_diameter(um)'] = 0.44*(np.sqrt(dataframe['Area'] / np.pi))
        return dataframe

    def name_fluor_columns(self, cell_trace_dfs_dict):

        for channel in cell_trace_dfs_dict.keys():
            print(channel)
            for df in cell_trace_dfs_dict[channel]:
                assert channel in df.Label.iloc[0], "Measurement .csv matched to wrong channel"
                df.rename(columns={'Mean': f'{channel}_mean',
                                   'Median': f'{channel}_median',
                                   'RawIntDen': f'{channel}_int'},
                          inplace=True)
                # Add corrected time and a cell diameter
                # column to the trace dataframe. Might re
                # add mean filters etc. trace tools at some
                # point. But for now keeping it simple
                df = self.correct_time(df)
                df = self.cell_diameter(df)

        return cell_trace_dfs_dict

    def get_processed_cell_trace_df(self, cell_index, master_index_df):
    
        dfs_dict = self.cell_trace_dfs_dict
        cell_channel_dfs = [dfs_dict[channel][int(cell_index)] for channel in self.fluor_channel_names]
        
        if len(cell_channel_dfs) == 0:
            print("No channel_dfs found")
            return None
        # If there's only one fluorescent channel collected
        elif len(cell_channel_dfs) == 1:
            final_cell_df = cell_channel_dfs[0]
        # If there are more than one fluoroscent channels
        # collected, as is normal.
        elif len(cell_channel_dfs) > 1:
            # Expanded dataframe merge
            base_df = cell_channel_dfs[0]
            for channel_df in cell_channel_dfs:
                cols_to_add = [col for col in channel_df.columns if col not in base_df.columns]
                for col in cols_to_add:
                    base_df.loc[:, col] = channel_df.loc[:, col]
            final_cell_df = base_df

        if final_cell_df.empty == True:
            m1 = f'cell_index: {cell_index}'
            m2 = f'compartment_dir: {master_index_df.compartment_dir[cell_index]}'
            print(f'Warning, returning empty cell_trace_df\n{m1}\n{m2}')

        return final_cell_df

    def get_processed_cell_trace_dfs(self):
        """
        Return a list of cell trace dataframes.
        All channels combined with added correct 
        time etc.
        """
        cell_trace_dfs = []
        if 'sub_coord' in self.master_index_df.columns:
            self.master_index_df.rename(columns={'sub_coord': 'cell_index'}, inplace=True)

        for cell_index in self.master_index_df.cell_index.unique():
            cell_trace_df = self.get_processed_cell_trace_df(cell_index, self.master_index_df)
            cell_trace_dfs.append(cell_trace_df)
    
        return cell_trace_dfs

    def get_cell_trace_paths(self):
        """
        Create a .csv filepath for each cell trace
        df in self.cell_trace_dfs, save the cell
        trace df at that path, and return the list
        of paths
        """
        master_index = self.master_index_df
        cells_dfs = self.cells_dfs
        cell_trace_dfs_dict = self.cell_trace_dfs_dict
    
        processed_traces = self.cell_trace_dfs
        trace_filepaths = []

        for cell_index in range(0, len(cells_dfs)):
            cell_master_df = cells_dfs[cell_index]
            cell_trace_df = processed_traces[cell_index]        
        # Find data to create a filename for this cell's trace
        # and save it as a csv
        # If the master index is old school, the compartment_dir
        # column will be called 'path'. 
            if 'path' in self.master_index_df.columns:
                datadir = os.path.abspath(master_index.loc[cell_index, 'path'])
            elif 'compartment_reldir' in master_index.columns:
                datadir = os.path.join(constants.byc_data_dir,
                                       master_index.loc[cell_index, 'compartment_reldir'])
            else:
                print(f"""
                    No comparment directory data in cell {cell_index}'s' master index

                    Either add a 'path' column containing the full path to the 
                    master index or make sure you are creating all ROIs for each
                    cell using the imagejpc plugin, 'save_cell_roi_set.py' which
                    constructs master index with proper formatting
                    """)

            date = self.master_index_df.loc[cell_index, 'date']
            expt_type= self.master_index_df.loc[cell_index, 'expt_type']
            try:
                compartment_name = datadir[datadir.rindex('\\')+1:]
            except:
                compartment_name = datadir[datadir.rindex('/')+1:]
            filename = f'{datadir}\\{date}_{expt_type}_{compartment_name}_cell{str(cell_index).zfill(3)}.csv'
            filepath = os.path.join(datadir, filename)
            cell_trace_df.to_csv(filepath, index=False)
            trace_filepaths.append(filepath)
    
        return trace_filepaths

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

