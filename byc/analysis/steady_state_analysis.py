import pandas as pd
import numpy as np
import tkinter as tk
import tkinter.filedialog as tkdia

def select_files(message):
    # return a list of files that user selects in the dialogue box
    root = tk.Tk()
    files = tkdia.askopenfilenames(parent=root, title=message)
    files_list = root.tk.splitlist(files)
    root.destroy()
    return sorted(files_list)

def select_file():
    # return a list of files that user selects in the dialogue box
    root = tk.Tk()
    file_path = tkdia.askopenfilename(parent=root, title='Choose files')
    root.destroy()
    return file_path

def make_dfs(fns):
    # return a list of dataframes created by reading the list of filenames
    # that you pass the function
    dfs = []
    for i in range(0, len(fns)):
        dfs.append(pd.read_csv(fns[i]))
    return dfs

def calc_bg_pix(df):
    # return a background value for an individual df
    
    # calculate total intensities of bg and cells
    bg_I = df.iloc[-1]['IntDen']
    cells_I = df.iloc[0:-1]['IntDen'].sum()
    
    bg_A = df.iloc[-1]['Area']
    cells_A = df.iloc[0:-1]['Area'].sum()
    
    bg_final = (bg_I - cells_I) / (bg_A - cells_A)
    return bg_final

def calc_bg_list_pix(dfs):
    # return a list of background values for each df
    # in the list of dfs passed to the function
    
    bg_list = []
    for i in range(0, len(dfs)):
        bg_list.append(calc_bg_pix(dfs[i]))
        
    return bg_list

def calc_bg_um(df):
    # return a background value for an individual df
    
    # calculate total intensities of bg and cells
    bg_I = df.iloc[-1]['IntDen']
    cells_I = df.iloc[0:-1]['IntDen'].sum()
    
    bg_A = df.iloc[-1]['Area']
    cells_A = df.iloc[0:-1]['Area'].sum()
    
    bg_final = (bg_I - cells_I) / (bg_A - cells_A)
    return bg_final

def calc_bg_list_um(dfs):
    # return a list of background values for each df
    # in the list of dfs passed to the function
    
    bg_list = []
    for i in range(0, len(dfs)):
        bg_list.append(calc_bg_um(dfs[i]))
        
    return bg_list

def add_norm_col(df, bg):
    # add a column to the df passed to this function that contains
    # bg normalized mean I values
    
    df['norm_mean'] = df['Mean'] - bg
    
def normalize(dfs, bgs):
    for i in range(0, len(dfs)):
        add_norm_col(dfs[i], bgs[i])

def drop_bg(df):
    # return the input df with the last row removed
    df_pruned = df.drop(len(df) - 1)
    return df_pruned

def drop_all_bgs(dfs):
    # return the input dfs with their last row removed
    
    pruned_dfs = []
    for i in range(0, len(dfs)):
        dropped = drop_bg(dfs[i])
        pruned_dfs.append(dropped)
        
    return pruned_dfs

def run_pix(channel_1, channel_2):
    
    """Return two dataframes (channel_1_zipper_dfs and channel_2_zipped_dfs).
       These dfs are read from csvs chosen by the user. """
    
    if type(channel_1) and type(channel_2) == str:
        
        channel_1_fns = select_files("Choose {} files for this condition".format(channel_1))
        channel_2_fns = select_files("Choose {} files for this condition".format(channel_2))

        channel_1_dfs = make_dfs(channel_1_fns)
        channel_2_dfs = make_dfs(channel_2_fns)
        # calclucate a background value (defined in function above)
        # for each FOV's dataframe
        channel_1_bgs = calc_bg_list_pix(channel_1_dfs)
        channel_2_bgs = calc_bg_list_pix(channel_2_dfs)
        # apply the backgorund subraction value
        normalize(channel_1_dfs, channel_1_bgs)
        normalize(channel_2_dfs, channel_2_bgs)
    
    channel_1_pruned_dfs = drop_all_bgs(channel_1_dfs)
    channel_2_pruned_dfs = drop_all_bgs(channel_2_dfs)
    
    channel_1_zipped_dfs = pd.concat(channel_1_pruned_dfs, ignore_index=True)
    channel_2_zipped_dfs = pd.concat(channel_2_pruned_dfs, ignore_index=True)
    
    return (channel_1_zipped_dfs, channel_2_zipped_dfs, channel_1, channel_2)

def finalize_pix(channel_1_name, channel_2_name):
    
    """ Return a completed df which combines data from channel_1 and channel_2 
        for a single population of cells from the fovs passed in run_pix """
    
    (channel_1_df, channel_2_df, channel_1_name, channel_2_name) = run_pix(channel_1_name, channel_2_name)
    
    df_final = channel_1_df
    
    df_final['bg_sub_mean_{}'.format(channel_1_name)] = channel_1_df['norm_mean']
    df_final['raw_mean_{}'.format(channel_1_name)] = channel_1_df['Mean']
    df_final['raw_int_{}'.format(channel_1_name)] = channel_1_df['IntDen']
    
    df_final['bg_sub_mean_{}'.format(channel_2_name)] = channel_2_df['norm_mean']
    df_final['raw_mean_{}'.format(channel_2_name)] = channel_2_df['Mean']
    df_final['raw_int_{}'.format(channel_2_name)] = channel_2_df['IntDen']
    
    return df_final

