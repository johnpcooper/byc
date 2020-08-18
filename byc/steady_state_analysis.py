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
    
    dfs_list = []
    for dataset_index in range(0, len(master_df)):
    
        info = master_df.loc[dataset_index, :]
        channel_names = info.fluor_channel_names.split()
        fovs = list(np.arange(1, max_n_fovs+1))

        fov_dfs_list = []
        for fov in fovs:
            channel_dfs = []
            try: # Find all the fov csv files that match the conditions specified in info.
                for channel_name in channel_names:

                    filename = f'{info.path}\\{info.expt_date}_{info.plasmid}_{info.genotype}_C{info.clone}_00{fov}_{channel_name}.csv'
                    channel_df = pd.read_csv(filename)
                    channel_df = channel_df.rename(columns={'Mean': str(f'{channel_name}_mean'), ' ': 'cell_index'})
                    channel_dfs.append(channel_df)

                fov_merged_df = reduce(lambda x, y: pd.merge(x, y, on='cell_index'), channel_dfs)
                fov_dfs_list.append(fov_merged_df)
            except:
                pass
                #print(f"No FOV {fov} for dataset {dataset_index}")
            
        try: # this try except statement assumes that this will only fail if pd.concat has
             # 0 objects to concatenate and this is because no files were found
            final_df = pd.concat(fov_dfs_list, ignore_index=True, sort=False)
        except:
            print(f"No file at : {filename}")
        # add identifying information to final dataset:
        for i in range(0, len(master_df.columns)):
            column_name = list(master_df.columns)[i]
            value = info[i]
            final_df.loc[:, f'{column_name}'] = value
        
        dfs_list.append(final_df)
        
    return dfs_list

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

    for plasmid_level in df.index.levels[0]:
        print(plasmid_level)

        n_cells = len(df.loc[plasmid_level, :])
        print(f"Number of cells in {plasmid_level} group = {n_cells}")
        proportion = n_cells / len(df)
        print(f"Fraction of all cell in {plasmid_level} = {proportion}")
        weight = 1 / proportion
        print(f"weight={weight}")
        df.loc[plasmid_level, 'weight'] = weight
        
    return df

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
        
    return df
