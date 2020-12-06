import os
import numpy as np
import pandas as pd
from scipy.stats import shapiro
import matplotlib.pyplot as plt
import matplotlib
import tkinter as tk
import tkinter.filedialog as tkdia
from functools import reduce


def select_files(prompt):
    # return a list of files that user selects in the dialogue box
    root = tk.Tk()
    files = tkdia.askopenfilenames(parent=root, title=prompt)
    files_list = root.tk.splitlist(files)
    root.destroy()
    return sorted(files_list)

def select_file():
    # return a list of files that user selects in the dialogue box
    root = tk.Tk()
    file_path = tkdia.askopenfilename(parent=root, title='Choose file')
    root.destroy()
    return file_path

def make_dfs(prompt):
    # return a list of dataframes created by reading the list of filenames
    # that you pass the function
    fns = select_files(prompt)
    dfs = []
    for i in range(0, len(fns)):
        dfs.append(pd.read_csv(fns[i]))
    return dfs, fns

def set_steady_state_dfs_list(master_df, max_n_fovs):
    
    """ Return a list of Dataframes, one for each distinct condition in the dataset
        e.g. plasmid, clone, genetic background """
    
    if max_n_fovs == None:
        max_n_fovs = 20

    dfs_list = []
    for dataset_index in range(0, len(master_df)):
    
        info = master_df.loc[dataset_index, :]
        channel_names = info.fluor_channel_names.split()
        fovs = list(np.arange(1, max_n_fovs+1))

        fov_dfs_list = []
        for fov in fovs:
            channel_dfs = []
            for channel_name in channel_names:
                condition_descriptors = [info.expt_date,
                                         info.plasmid,
                                         info.genotype,
                                         str(info.tet_concn).zfill(3)+'uM-Tet',
                                         'clone' + str(info.clone),
                                         str(fov).zfill(3),
                                         channel_name]
                filename = '_'.join(str(desc) for desc in condition_descriptors) + '.csv'
                filepath = os.path.join(info.path, filename)
                
                if os.path.exists(filepath):
                    print(f"Found data at {filepath}")
                    channel_df = pd.read_csv(filepath)
                    channel_df = channel_df.rename(columns={'Mean': str(f'{channel_name}_mean'), ' ': 'cell_index'})
                    channel_df = channel_df.rename(columns={'RawIntDen': str(f'{channel_name}_int'), ' ': 'cell_index'})
                    channel_dfs.append(channel_df)
                    fov_merged_df = reduce(lambda x, y: pd.merge(x, y, on='cell_index'), channel_dfs)
                    fov_dfs_list.append(fov_merged_df)
                else:
                    pass           

        final_df = pd.concat(fov_dfs_list, ignore_index=True, sort=False)
            
        # add identifying information to final dataset:
        for i in range(0, len(master_df.columns)):
            column_name = list(master_df.columns)[i]
            value = info[i]
            final_df.loc[:, f'{column_name}'] = value
        
        dfs_list.append(final_df)
        
    return dfs_list


def make_expt_df(master_index_path, bg_channel='yfp', filter_cells=False):
    """
    Find a {exptdate}_master_index.csv file at master_index_df_path,
    read in all steady state imaging measurement .csvs found 
    recorded in that master_index as dataframes. Then take some ratios,
    clean, etc. that data and return a concatenated dataframe of all
    those measurement csvs
    """
    master_df = pd.read_csv(master_index_path)
    all_dfs_list = set_steady_state_dfs_list(master_df, max_n_fovs=20)
    all_data_df = pd.concat(all_dfs_list, sort=False, ignore_index=True)
    
    # Set a no-plasmid slice to define background values later
    no_plasmid = all_data_df[all_data_df.plasmid == 'no-plasmid']
    
    # Get the fluor channel names manually recorded in master index
    fluor_channels = master_df.fluor_channel_names.iloc[0].split(' ')
    
    # Set background autofluorescence values for each channel and normalize
    for channel in fluor_channels:
        
        # Normalize signal to median of integrated fluorescence within cell 
        # for no plasmid, BY4741 cells
        channel_bg = no_plasmid[f'{channel}_mean'].median()
        all_data_df.loc[:, f'{channel}_norm'] = all_data_df[f'{channel}_mean']/channel_bg
        
        # Normalize signal to median of integrated fluorescence within cell 
        # for no plasmid, BY4741 cells
        channel_bg = no_plasmid[f'{channel}_int'].median()
        all_data_df.loc[:, f'{channel}_int_norm'] = all_data_df[f'{channel}_int']/channel_bg
        
    # Set coloumns containing ratios of each channel normalized
    # and raw to each other channel
    for channel in fluor_channels:
        
        other_channels = [name for name in fluor_channels if name != channel]
        
        for channel2 in other_channels:
            
            # Ratios of normalized mean within cell channel signals
            ratios = all_data_df[f'{channel}_norm'] / all_data_df[f'{channel2}_norm']
            all_data_df.loc[:, f'{channel}_{channel2}'] = ratios
            # Ratios of raw mean within cell channel signals
            rawratios = all_data_df[f'{channel}_mean'] / all_data_df[f'{channel2}_mean']
            all_data_df.loc[:, f'{channel}_{channel2}_raw'] = rawratios
            # Ratios of normalized integrated within cell channel signals
            ratios = all_data_df[f'{channel}_int_norm'] / all_data_df[f'{channel2}_int_norm']
            all_data_df.loc[:, f'{channel}_{channel2}_int'] = ratios
            # Ratios of raw integrated within cell channel signals
            rawratios = all_data_df[f'{channel}_int'] / all_data_df[f'{channel2}_int']
            all_data_df.loc[:, f'{channel}_{channel2}_int_raw'] = rawratios

    # Filter out background expr cells
    if filter_cells:

        std = all_data_df.loc[all_data_df.plasmid=='no-plasmid', f'{bg_channel}_mean'].std()
        med = all_data_df.loc[all_data_df.plasmid=='no-plasmid', f'{bg_channel}_mean'].median()
        thresh = med + 2*std
        print(f'No-plasmid {bg_channel} median + 2*stdev = {thresh}')
        all_data_df = all_data_df[all_data_df[f'{bg_channel}_mean'] > thresh]
        
    else:
        pass
        
    return all_data_df

def set_flow_cyto_dfs_list(master_df):
    
    """ Return a list of Dataframes, one for each distinct condition in the dataset
        e.g. plasmid, clone, genetic background """
    
    n_datasets_found = 0
    dfs_list = []
    for dataset_index in range(0, len(master_df)):
    
        info = master_df.loc[dataset_index, :]
        channel_names = info.fluor_channel_names.split()
        dataset_id = f'{info.plasmid}_{info.genotype}_{info.clone}'
        
        try:
            filepath = f'{info.path}\\{info.file_name}'
            df = pd.read_csv(filepath)
            print(f"Found .csv for {dataset_id} at {filepath}\n")
            n_datasets_found += 1
            
        except:                
            print(f"No file found at {filename} for dataset {dataset_id}\n")

        # add identifying information to final dataset:
        for i in range(0, len(master_df.columns)):
            column_name = list(master_df.columns)[i]
            value = info[i]
            df.loc[:, f'{column_name}'] = value
        
        dfs_list.append(df)
        
    if len(dfs_list) != len(master_df):
        print("WARNING, did not find a .csv for all rows in master_df")
    else:
        print("Found .csv for all rows in master_df")
        
    return dfs_list

def set_proportional_weights_by_plasmid(df):
    
    """ Return the dataframe passed to this function with a new column called 'weight'. 
        The weight for a row (cell) is 1 - the number of cells with that cell's unique
        plasmid / total number of cells in the data set. 
        
        This allows evenly selecting from each plasmid (or potentially other) group 
        when using df.sample(n=some_number, weights=df.weight). 
        
        WARNING: currently this function applies weight by df.index.levels[0]"""
    
    df.loc[:, 'weight'] = 0
    df = df.set_index(['plasmid', 'cell_index'])

    for plasmid_level in df.index.levels[0]:
        print(plasmid_level)

        n_cells = len(df.loc[plasmid_level, :])
        print(f"Number of cells in {plasmid_level} group = {n_cells}")
        proportion = n_cells / len(df)
        print(f"Fraction of all cell in {plasmid_level} = {proportion}")
        weight = 1 / proportion
        print(f"weight={weight}")
        df.loc[plasmid_level, 'weight'] = weight
        
    return df.reset_index()

def set_proportional_weights(df, by=['plasmid', 'genotype', 'clone']):
    
    """ Set the index of the dataframe using the list 
        provided in by. Count the number of cells in 
        each unique group according to that index and set
        a weight value = 1 / proportion of total cells
        in that unique group"""
    
    df.loc[:, 'weight'] = 0
    df = df.set_index(by)

    for group in df.index.unique():
        print(group)
        n_cells = len(df.loc[group, :])
        print(f"Number of cells in {group} group = {n_cells}")
        if n_cells != 0:
            proportion = n_cells / len(df)
            print(f"Fraction of all cells in {group} = {proportion}")
            weight = 1 / proportion
            print(f"weight={weight}")
            df.loc[group, 'weight'] = weight
        else:
            print(f'No cells found in {group}')
            df.loc[group, 'weight'] = 0
        
    return df.reset_index()


