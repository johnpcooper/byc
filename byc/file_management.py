import pandas as pd
import tkinter as tk
import tkinter.filedialog as tkdia
import os
from skimage import io
from tifffile import imsave
from skimage.util import img_as_uint

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

def rename_steady_state(max_n_fovs):

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
                print(f'Attempted target directory: {fov_path}')
                pass
            

def set_xy_dir_names(expt_dir, expt_name):
    """ Return a list of xy position directory names (not paths). Rename the list of xy positions
        that have been created in micromanager (should be Pos0, Pos1, Pos2...) to looks like
        expt_name_xy001..."""

    xy = 0
    xy_dir_names = []
    for pos_dir in os.listdir(expt_dir):
        xy_dir_name = f'{expt_name}_xy{xy:02}'
        try:        
            os.rename(f'{expt_dir}//{pos_dir}', f'{expt_dir}//{expt_name}_xy{xy:02}')
        except:
            print(f"File already exists {xy_dir_name}")

        xy_dir_names.append(xy_dir_name)
        xy += 1
        
    return xy_dir_names

def reshape_timepoints(xy_dir_names, expt_dir, n_channels):
    
    """ Rename and combine individual channel .tifs for each timepoint in the experiment. Output
        Shape is (height, width, n_channels). This allows files to be ready and aligned using 
        byc alignment code. """
    
    for xy_dir_name in xy_dir_names:

        fns_list = os.listdir(f'{expt_dir}//{xy_dir_name}')

        timepoint = 0
        for i in range(0, len(fns_list), n_channels):

            timepoint_channel_fns = fns_list[i:i + n_channels]
            timepoint_channel_paths = [f'{expt_dir}//{xy_dir_name}//{tp_channel_fn}' for tp_channel_fn in timepoint_channel_fns]
            print(f'Timepoint {timepoint:03} channels {timepoint_channel_fns}')

            timepoint_collection = io.imread_collection(timepoint_channel_paths)
            timepoint_stack = io.concatenate_images(timepoint_collection)

            imsave(f'{expt_dir}//{xy_dir_name}//{xy_dir_name}_t{timepoint:03}.tif',
                   img_as_uint(timepoint_stack),
                   shape=(timepoint_stack.shape[0], timepoint_stack.shape[1], n_channels))

            for file_path in timepoint_channel_paths:
                os.remove(file_path)

            timepoint += 1
            
def rename_byc():
    
    expt_dir = select_directory("Choose the directoy holding the micromanager output of your byc experiment")
    pos_dir_list = os.listdir(expt_dir)
    
    channels = ['Brightfield', 'YFP', 'RFP']
    fluor_names = ['bf', 'yfp', 'dsred']
    channel_dict = dict(zip(fluor_names, channels))
    expt_name = '20200214_byc'
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

