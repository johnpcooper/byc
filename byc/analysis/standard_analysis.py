import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.filedialog as tkdia

plt.style.use('default')
matplotlib.rcParams['font.sans-serif'] = 'Arial'

class Cell_Data(object):
    
    """ Doc Strings """
    
    def __init__(self):
        
        # set the files that will be used for every cell analyzed
        self.master_cells_df = self.set_master_cells_df("Choose the .csv for the master index of this experiment")
        self.cells_dfs = self.set_cells_dfs(self.master_cells_df)
        
        # set lists containing senescent slice for each cell
        self.senescent_slices = self.set_senenescent_slices(self.cells_dfs)
        
        # set lists containing dataframes with measurements for each cell in cells_dfs
        self.raw_dsred_trace_dfs, self.raw_yfp_trace_dfs = self.set_raw_cell_trace_dfs(self.cells_dfs, self.senescent_slices)
        #self.processed_dsred_trace_dfs, self.processed_trace_dfs = self.set_processed_cell_trace_dfs(self.cells_dfs, self.senescent_slices)
    
    def set_fp(self, prompt):

        """ Return the path to a file of interest. Call with the prompt you would like to 
            display in the dialog box."""

        # create the dialog box and set the fn
        root = tk.Tk()
        fp = tkdia.askopenfilename(parent=root, title=prompt)
        root.destroy() # very important to destroy the root object, otherwise python 
        # just keeps running in a dialog box

        return fp # return the path to the file you just selected

    def set_master_cells_df(self, prompt):

        """ Return a dataframe read from the master cells index .csv"""

        # define the path to the index .csv
        master_cells_fp = self.set_fp(prompt)
        # define the filename for the master expt index
        master_cells_df = pd.read_csv(master_cells_fp)

        return master_cells_df

    def set_cells_dfs(self, master_cells_df):

        """ Return a list of sub-dataframes of the master_cells_df. The sub-dataframes
           contain data for only one cell according to its sub_coord """

        # create the cells_dfs list
        cells_dfs = []

        # add to the cells_dfs a dataframe at index i for every
        # unique value in the master_cells_df['sub_coord']
        for i in master_cells_df['sub_coord'].unique():
            # set the logic mask to a sub-dataframe of master_cells_df containing
            # only values whose 'sub_coord' value is value
            print("value is ", i)
            logic_mask = (master_cells_df['sub_coord'] == i)
            cells_dfs.append(master_cells_df[logic_mask])
            print("cells_dfs is now %d elements long" % len(cells_dfs))

        return cells_dfs
    
    def set_senenescent_slices(self, cells_dfs):
        
        """ Return a list of number indicating where along the expt each cell went
            senescent (ended its last division). If senescence wasn't observed
            for the cell, then False will be added to the list. """
        
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

    def set_raw_cell_trace_dfs(self, cells_dfs, senescent_slices):

        """ Return two lists of dataframes (dsred, yfp) containing traces for each cell 
            in cells_dfs. The path for each trace is constructed based on how roi_management names
            and saves tif stacks. """

        cell_traces_dsred = []
        cell_traces_yfp = []

        for cell_index in range(0, len(cells_dfs)):

            cell_df = cells_dfs[cell_index]

            path = str(cell_df['path'][cell_df.index.min()] + "\\")
            expt_date = cell_df['date'][cell_df.index.min()]
            expt_type = cell_df['expt_type'][cell_df.index.min()]
            xy = cell_df['xy'][cell_df.index.min()]
            sub_coord = cell_df['sub_coord'][cell_df.index.min()]

            stack_title = str("%s_%s_xy%s_cell%s" % (expt_date, expt_type, xy, sub_coord))

            # read the .csv containing measurements
            cell_trace_dsred_path = str(path + stack_title + "dsred.csv")
            cell_trace_yfp_path = str(path + stack_title + "yfp.csv")
            
            # create each dataframe for dsred and yfp
            cell_trace_dsred = pd.read_csv(cell_trace_dsred_path)
            cell_trace_yfp = pd.read_csv(cell_trace_yfp_path)
            # create the time to senescence column
            if not (senescent_slices[cell_index] == "FALSE"):
                
                cell_trace_dsred['slices_to_senescence'] = cell_trace_dsred['Slice'] - int(senescent_slices[cell_index])
                cell_trace_yfp['slices_to_senescence'] = cell_trace_yfp['Slice'] - int(senescent_slices[cell_index])
                
            else:
                print("Senescence not observed for cell %s" % cell_index)
                
            cell_traces_dsred.append(cell_trace_dsred)
            cell_traces_yfp.append(cell_trace_yfp)

        return cell_traces_dsred, cell_traces_yfp
    
    def set_processed_cell_trace_dfs(self, cells_dfs):
        pass

# define a list that will be filled with indices of cells that are in the same slice
def sliding_median_window(dataframe, window_width):
    
    """ Return the input DataFrame with a column ("sliding_median") added containing a 
        sliding mean sampling window with width specified. """
    
    dataframe.loc[:, 'sliding_median'] = np.zeros(len(dataframe['Slice']))
    
    for i in dataframe.index:
        if i == 0:
            dataframe.loc[:, 'sliding_median'][i] = dataframe['Mean'][i: window_width + 1].median()
        
        elif 0 < i < window_width:
            dataframe.loc[:, 'sliding_median'][i] = dataframe['Mean'][0: i + window_width + 1].median()
                    
        elif i >= window_width:
            dataframe.loc[:, 'sliding_median'][i] = dataframe['Mean'][i - window_width: i + window_width + 1].median()
            
        else:
            print("Somehow reached an index value out of range. This is the current index:", i)
            
    return dataframe

# define a list that will be filled with indices of cells that are in the same slice
def sliding_mean_window(dataframe, window_width):
    
    """ Return the input DataFrame with a column ("sliding_mean") added containing a 
        sliding mean sampling window with width specified. """
    
    dataframe.loc[:, 'sliding_mean'] = np.zeros(len(dataframe['Slice']))
    
    for i in dataframe.index:
        if i == 0:
            dataframe.loc[:, 'sliding_mean'][i] = dataframe['Mean'][i: window_width + 1].mean()
        
        elif 0 < i < window_width:
            dataframe.loc[:, 'sliding_mean'][i] = dataframe['Mean'][0: i + window_width + 1].mean()
                    
        elif i >= window_width:
            dataframe.loc[:, 'sliding_mean'][i] = dataframe['Mean'][i - window_width: i + window_width + 1].mean()
            
        else:
            print("Somehow reached an index value out of range. This is the current index:", i)
            
    return dataframe

# define a list that will be filled with indices of cells that are in the same slice
def diameter_sliding_window(dataframe, window_width):
    
    """ Return the input DataFrame with a column ("cell_diamter_sliding") added containing a 
        sliding mean sampling window with width specified. """
    
    dataframe.loc[:, 'cell_diameter_sliding'] = np.zeros(len(dataframe['Slice']))
    
    for i in dataframe.index:
        if i == 0:
            dataframe.loc[:, 'cell_diameter_sliding'][i] = dataframe['cell_diameter(um)'][i: window_width + 1].mean()
        
        elif 0 < i < window_width:
            dataframe.loc[:, 'cell_diameter_sliding'][i] = dataframe['cell_diameter(um)'][0: i + window_width + 1].mean()
                    
        elif i >= window_width:
            dataframe.loc[:, 'cell_diameter_sliding'][i] = dataframe['cell_diameter(um)'][i - window_width: i + window_width + 1].mean()
            
        else:
            print("Somehow reached an index value out of range. This is the current index:", i)
            
    return dataframe

def correct_time(dataframe, collection_interval):        
    
    """ Return a DataFrame of the DataFrame passed to correct_time() with
        columns added for minutes and hours according to the collection_interval argument """
    
    if "slices_to_senescence" in dataframe.columns:        
        dataframe.loc[:, 'minutes_to_senescence'] = dataframe['slices_to_senescence'] * collection_interval
        dataframe.loc[:, 'hours_to_senescence'] = dataframe['minutes_to_senescence'] / 60
        
    else:
        print("Senescence not observed for this cell")
    
    dataframe.loc[:, 'minutes'] = dataframe['Slice'] * collection_interval
    dataframe.loc[:, 'hours'] = dataframe['minutes'] / 60    
    
    return dataframe

def cell_diameter(dataframe):
    """ Return a DataFrame with added column for cell diameter in um """
    
    dataframe.loc[:, 'cell_diameter(um)'] = 0.44*(np.sqrt(dataframe['Area'] / np.pi))
    return dataframe

def complete_dataframe(dsred_df, yfp_df, cell_number, window_width, collection_interval):
    dsred_df = dsred_df.copy()
    yfp_df = yfp_df.copy()
    
    # Add time in minutes and hours
    final = correct_time(dsred_df, collection_interval)
    # Add sliding mean window column
    dsred_df = sliding_mean_window(dsred_df, window_width)
    yfp_df = sliding_mean_window(yfp_df, window_width)
    # add sliding median window column
    dsred_df = sliding_median_window(dsred_df, window_width)
    yfp_df = sliding_median_window(yfp_df, window_width)
    # Add cell diameter in um column
    dsred_df = cell_diameter(dsred_df)
    yfp_df = cell_diameter(yfp_df)
    # Add a cell diameter sliding column
    dsred_df = diameter_sliding_window(dsred_df, window_width)
    yfp_df = diameter_sliding_window(yfp_df, window_width)
    
    # Put dsred and yfp dataframes together
    final.loc[:, 'dsred_mean'] = dsred_df['Mean']
    final.loc[:, 'yfp_mean'] = yfp_df['Mean']
    
    final.loc[:, 'dsred_sliding_mean'] = dsred_df['sliding_mean']
    final.loc[:, 'yfp_sliding_mean'] = yfp_df['sliding_mean']
    
    final.loc[:, 'dsred_sliding_median'] = dsred_df['sliding_median']
    final.loc[:, 'yfp_sliding_median'] = yfp_df['sliding_median']
    
    final.loc[:, 'yfp/dsred'] = yfp_df['Mean'] / dsred_df['Mean']
    final.loc[:, 'yfp/dsred_sliding'] = yfp_df['sliding_mean'] / dsred_df['sliding_mean']
    
    return final

def set_processed_trace_dfs(cells_dfs, window_width, collection_interval, raw_dsred_trace_dfs, raw_yfp_trace_dfs):
    
    """ Return a list of processed dataframes for each cell trace. The processed
        dataframes have added columns like adjusted time, ratios, etc. """
    
    processed_trace_dfs = []
    
    for cell_index in range(0, len(cells_dfs)):
        print("Processing trace dataframe for cell %s of %s" % (cell_index, len(cells_dfs) - 1))
        processed_trace_df = complete_dataframe(raw_dsred_trace_dfs[cell_index],
                                                raw_yfp_trace_dfs[cell_index],
                                                cell_index,
                                                window_width,
                                                collection_interval)

        processed_trace_dfs.append(processed_trace_df)
        
    return processed_trace_dfs

def set_processed_traces(window_width, collection_interval):
    
    cd = Cell_Data()
    master_index = cd.master_cells_df
    cells_dfs = cd.cells_dfs
    raw_dsred_trace_dfs = cd.raw_dsred_trace_dfs
    raw_yfp_trace_dfs = cd.raw_yfp_trace_dfs
    
    processed_traces = set_processed_trace_dfs(cells_dfs, window_width, collection_interval, raw_dsred_trace_dfs, raw_yfp_trace_dfs)
    
    return processed_traces

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

