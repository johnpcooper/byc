import numpy as np
import pandas as pd
import tkinter as tk
import tkinter.filedialog as tkdia
import os
import re
import ast
from skimage import io
from tifffile import imsave
from skimage.util import img_as_uint
from read_roi import read_roi_zip

def select_files(prompt):
    # return a list of files that user selects in the dialogue box
    root = tk.Tk()
    files = tkdia.askopenfilenames(parent=root, title=prompt)
    files_list = root.tk.splitlist(files)
    root.destroy()
    return sorted(files_list)

def select_file(prompt):
    # return a list of files that user selects in the dialogue box
    root = tk.Tk()
    file_path = tkdia.askopenfilename(parent=root, title=prompt)
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

def select_directory(prompt):
    # Choose the directory holding all the fields of view that you'll align
    root = tk.Tk()
    input_dir = tkdia.askdirectory(parent=root,
                                 title=prompt)
    root.destroy()
    
    return input_dir

def read_roi_position_indices(path):

    roi_set = read_roi_zip(path)
    keys = list(roi_set.keys())
    
    # define an array of bud start positions for this cell
    bud_positions = []
    for i in range(0, len(keys)):
        bud_positions.append(roi_set[keys[i]]['position']-1)
        
    return np.array(bud_positions)

def find_condition(descriptor_list, conditions):
    """
    Search conditions (a list of directory names) for 
    one that contains all the strings in descriptor list
    (gotten from values in some master_index) and return
    the condition (directory name) that matches every
    descriptor in descriptor list
    
    If 0 or > 1 matches found, return None
    """
    matched_conditions = []
    # Should be regex but this works
    for condition in conditions:
        match_results = [str(descriptor) in condition for descriptor in descriptor_list]
        if False not in match_results:
            matched_conditions.append(condition)

    if len(matched_conditions) == 1:
        return matched_conditions[0]
    else:
        print(f'Failed to find a single directory that matched descriptor list: \n{descriptor_list}')
        print(f'Found {len(matched_conditions)} matched directories:\n{matched_conditions}')
        return None
    
def find_fov_paths(conditionpath):
    """
    Return the list of fov directory names found
    in the directory at conditionpath
    """
    assert os.path.isdir(conditionpath), f'{conditionpath} is not a directory'    
    
    items = os.listdir(conditionpath)
    matches = [name for name in items if 'Pos' in name]
    
    fov_paths = [os.path.join(conditionpath, name) for name in matches]
    
    if len(fov_paths) > 0:
        errmsg = f'Found a non-directory fov path in fov_paths:\n{fov_paths}'
        assert False not in [os.path.isdir(path) for path in fov_paths], errmsg
        return fov_paths
    else:
        print(f'Found no fov paths at conditionpath:\n{conditionpath}')
        return None
    
def rename_channels(fov_path, channels_dict, base_filename, exptdir):
    """
    Move and rename each channel .tif file according
    to arguments passed.
    
    Return nothing
    """
    channels = [path for path in os.listdir(fov_path) if '.tif' in path]
    fov = fov_path[fov_path.rindex('Pos') + len('Pos'):]
    found_channels = []
    for fov_channel_filename in channels:
        
        for key in channel_dict.keys():

            if channel_dict[key] in fov_channel_filename:

                src = os.path.join(fov_path, fov_channel_filename)
                dst = os.path.join(exptdir, f'{base_filename}_{fov.zfill(3)}_{key}.tif')
                shutil.copyfile(src, dst)

                found_channels.append(key)
    
    # Check if we found a file for every channel dictated in the master_index
    if len(found_channels) == len(channel_dict):
        pass
    else:
        print(f"Couldn't find a file for every channel in:\n{channels_dict.keys}. Only found files for the following: {found_channels}")

def rename_steady_state():
    """
    Ask the user to choose a master_index for their steady
    state experiment. Then rename files tifs based
    on that master index
    """

    master_index_df = pd.read_csv(files.select_file("Choose master index .csv"))
    exptdir = master_index_df.path.iloc[0]
    conditions = os.listdir(exptdir)

    for condition_index in master_index_df.index:

        row = master_index_df.loc[condition_index, :]
        # Should make it easier to change which columns
        # are used to match conditions to directories,
        # but for now this list works
        descriptor_list = [row.expt_date,
                           row.plasmid,
                           row.genotype,
                           f'clone{row.clone}']

        # Check if there's a directory that matches the 
        # descriptors provided for this condition in 
        # the master index
        condition = find_condition(descriptor_list, conditions)
        if condition:
            base_filename = f'{row.expt_date}_{row.plasmid}_{row.genotype}_C{row.clone}'
            base_path = f'{exptdir}\\{base_filename}'
            channel_dict = dict(zip(row.channel_names.split(), row.raw_channel_names.split()))
            conditionpath = os.path.join(exptdir, condition)
            print(f'Processing files found in condition dir:\n{conditionpath}')
        else:
            print(f'No directories found matching condition {condition_index} in master index')

        # Look for FOV directories in the conditionpath
        fov_paths = find_fov_paths(conditionpath)
        if fov_paths != None:
            for fov_path in fov_paths:
                rename_channels(fov_path, channel_dict, base_filename, exptdir)
        else:
            print(f"No FOV directories found at conditionpath {conditionpath}")
            
def rename_steady_state_legacy(max_n_fovs):
    """
    This verision of rename_steady_state() works with data
    collected on or before 20200312
    """

    master_index_df = pd.read_csv(select_file("Choose master index for this experiment"))
    fovs_list = range(1, max_n_fovs+1)

    for measurement_index in master_index_df.index:

        info = master_index_df.loc[measurement_index, :]
        source_dir = info.path
        base_filename = f'{info.expt_date}_{info.plasmid}_{info.genotype}_C{info.clone}'
        print(f'Evaluating data for measurement {base_filename}')
        base_path = f'{source_dir}\\{base_filename}'

        channel_dict = dict(zip(info.channel_names.split(), info.raw_channel_names.split()))

        for fov in fovs_list:
            try:
                fov_path = f'{base_path}_{fov}\\Pos0'
                fov_channel_filenames = os.listdir(fov_path)
                print(f'Found a directory for {base_filename}_00{fov}')

                found_channels = []
                for fov_channel_filename in fov_channel_filenames:

                    for key in channel_dict.keys():
                        
                        if channel_dict[key] in fov_channel_filename:

                            print(f'Matching {fov_channel_filename} to {key}')
                            src = f'{fov_path}\\{fov_channel_filename}'
                            dst = f'{source_dir}\\{base_filename}_00{fov}_{key}.tif'
                            os.rename(src, dst)

                            found_channels.append(key)
                            
                        else:
                            pass

                # Loop to check if we found a file for every channel dictated in the master_index
                if len(found_channels) == len(channel_dict):
                    pass
                else:
                    print(f"Couldn't find a file for every channel. Only found files for the following: {found_channels}")



            except:
                print(f'No fov {fov} for {base_filename}')
                print(f'Attempted target path: {fov_path}')
                pass
            

def set_xy_dir_names(expt_dir, expt_name, **kwargs):
    """ 
    Return a list of xy position directory names (not paths). Rename the list of xy positions
    that have been created in micromanager (should be Pos0, Pos1, Pos2...) to look like
    expt_name_xy01...
    """
    rename = kwargs.get('rename', True)
    xy_dir_names = []
    fov_dir_pattern = 'Pos'
    # Iterate through the expt_dir and look for directories
    # that match the fov_dir_pattern. When found, rename the
    # directory and append its name to the list of xy_dir_names
    for dirpath, dirnames, filenames in os.walk(expt_dir):

        match = re.search(fov_dir_pattern, dirpath)
        if match:            
            xy_str = match.string[match.string.rindex(fov_dir_pattern)+len(fov_dir_pattern):]
            xy = int(xy_str)
            xy_dir_name = f'{expt_name}_xy{str(xy).zfill(2)}'
            new_xy_path = os.path.join(expt_dir, xy_dir_name)
            xy_dir_names.append(xy_dir_name)
            if rename:
                os.rename(dirpath, new_xy_path)

    return xy_dir_names

def reshape_timepoints(xy_dir_names, expt_dir, n_channels):    
    """ 
    For each xy FOV, combine individual channel .tifs for each timepoint. Output
    Shape is (height, width, n_channels). This allows files to be read and aligned using 
    byc.process.
    """    
    for xy_dir_name in xy_dir_names:

        all_fns_list = os.listdir(os.path.join(expt_dir, xy_dir_name))
        fns_list = []
        # Refine files in xy dir to only .tifs, there will likely
        # be some .txt etc. metadata in the directory that we don't need
        fns_list = [fn for fn in all_fns_list if fn.endswith('.tif')]

        timepoint = 0
        for i in range(0, len(fns_list), n_channels):

            timepoint_channel_fns = fns_list[i:i + n_channels]
            timepoint_channel_paths = [f'{expt_dir}//{xy_dir_name}//{tp_channel_fn}' for tp_channel_fn in timepoint_channel_fns]
            print(f'Timepoint {timepoint:03} channels {timepoint_channel_fns}')

            timepoint_channels = [io.imread(fn) for fn in timepoint_channel_paths]
            timepoint_stack = io.concatenate_images(timepoint_channels)

            imsave(f'{expt_dir}//{xy_dir_name}//{xy_dir_name}_t{timepoint:03}.tif',
                   img_as_uint(timepoint_stack),
                   shape=(timepoint_stack.shape[0], timepoint_stack.shape[1], n_channels))

            for file_path in timepoint_channel_paths:
                os.remove(file_path)

            timepoint += 1

def get_byc_display_and_comments_dict(expt_dir):
    """
    Return a dictionary generated by literal 
    interpretion of the 'display_and_comments.txt'
    file found in expt_dir
    """
    filepath = os.path.join(expt_dir, 'display_and_comments.txt')
    with open(filepath, 'r') as file:
        contents = file.read()
        dictionary = ast.literal_eval(contents)

    return dictionary
            
def rename_byc(**kwargs):
    
    expt_dir = kwargs.get('expt_dir', None)
    if expt_dir == None:
        expt_dir = select_directory("Choose the directoy holding the micromanager output of your byc experiment")
    else:
        pass
    expt_name = kwargs.get('expt_name', expt_dir.split(sep='/')[-1])
    # Get channels from metadata
    metadata = get_byc_display_and_comments_dict(expt_dir)
    channels = kwargs.get('channels', [channel['Name'].lower() for channel in metadata['Channels']])
    assert type(channels) == list, "channels must be a list of strings"
    n_channels = len(channels)
    
    xy_dir_names = set_xy_dir_names(expt_dir, expt_name)
    reshape_timepoints(xy_dir_names, expt_dir, n_channels)

def set_file_paths(prompt):
    # create the dialog box and set the fn
    root = tk.Tk()
    fps = tkdia.askopenfilenames(parent=root, title=prompt)
    root.destroy() # very important to destroy the root object, otherwise python 
    # just keeps running in a dialog box

    return fps # return the path to the file you just selected

def get_dfs_list():
    """
    Ask the user to choose a list of .csv files that will be read in as
    dataframes. Return (the list of dataframes, list of csv paths)
    """
    fps = set_file_paths("Choose the .csv files for this expt condition")

    if all(['.csv' == path[-4:] for path in fps]):    
        cell_trace_dfs_list = []
        for i in range(0, len(fps)):
            df = pd.read_csv(fps[i])
            cell_trace_dfs_list.append(df)
            
        return cell_trace_dfs_list, fps

    else:
        print('Selection included non-csvs, no dfs and paths lists constructed')
        return None




