import os
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import tkinter as tk
import tkinter.filedialog as tkdia
import ast
from skimage import io
from tifffile import imsave
from skimage.util import img_as_uint
from read_roi import read_roi_zip

from byc import constants

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
    """
    Read the .zip or .roi Fiji ROI file at <path>
    and return the 0 based index for each position
    in that file (eg position 1 becomes index 0)

    Return the position indices as a numpy array
    """
    roi_set = read_roi_zip(path)
    keys = list(roi_set.keys())
    
    # define an array of bud start positions for this cell
    bud_positions = []
    for i in range(0, len(keys)):
        bud_positions.append(roi_set[keys[i]]['position']-1)
        
    return np.array(bud_positions)

def add_abspath_to_df(df, colname='crop_roi_set_relpath'):
    """
    Find <colname> in df and append constants.byc_data_dir to 
    the beginning of the <colname> value to make it a full path
    relative to the installation directory of byc

    Return nothing, new column will be added to dataframe in place
    """
    fullpaths = constants.byc_data_dir + df[colname].astype(str)
    fullpaths = [os.path.abspath(path) for path in fullpaths]
    screenedpaths = []
    nonexisting_paths = []
    for p in fullpaths:
        if os.path.exists(p):
            screenedpaths.append(p)
        else:
            if p not in nonexisting_paths:
                nonexisting_paths.append(p)
            screenedpaths.append(p)
    if len(nonexisting_paths)>0:
        print(f'Did not find {len(nonexisting_paths)} following paths')
        for p in nonexisting_paths:
            print(p)
    
    df.loc[:, colname.replace('relpath', 'abspath')] = screenedpaths

def add_cell_channel_crop_stack_paths(df, channels=['bf', 'yfp']):
    
    dates = df.date.astype(int).astype(str)
    types = df.expt_type
    xys = df.xy.astype(int).astype(str).str.zfill(2)
    cell_idxs = df.cell_index.astype(int).astype(str).str.zfill(3)
    compartment_reldir = df.compartment_reldir.iloc[0]

    for channel in channels:
        filenames = dates + '_' + types + '_xy' + xys + '_cell' + cell_idxs + f'_{channel}_stack.tif'
        filerelpaths = [os.path.join(compartment_reldir, fn) for fn in filenames]
        fileabspaths = [os.path.join(constants.byc_data_dir, filerelpath) for filerelpath in filerelpaths]
        fileabspaths = [os.path.abspath(path) for path in list(fileabspaths)]

        df.loc[:, f'{channel}_crop_stack_path'] = fileabspaths

def add_cell_channel_xy_source_stack_paths(df, channels=['bf', 'yfp']):
    
    dates = df.date.astype(int).astype(str)
    types = df.expt_type
    xys = df.xy.astype(int).astype(str).str.zfill(2)
    compartment_reldir = df.compartment_reldir.iloc[0]

    for channel in channels:
        filenames = dates + '_' + types + '_xy' + xys + f'_{channel}_stack.tif'
        filerelpaths = [os.path.join(compartment_reldir, fn) for fn in filenames]
        fileabspaths = [os.path.join(constants.byc_data_dir, filerelpath) for filerelpath in filerelpaths]
        fileabspaths = [os.path.abspath(path) for path in list(fileabspaths)]

        df.loc[:, f'{channel}_stack_path'] = fileabspaths

def add_cell_measurement_roi_paths(df, channels=['bf', 'yfp'], source_col_suffix='_crop_stack_abspath'):
    
    if channels==None:
        channels = constants.all_channel_names

    to_replace = '_crop_rois.zip'
    # For each row in the datafrme, we should be able to find
    # one measurement_roi at its path unless this experiment
    # was automatically segmented. Need to merge the lists of 
    # filepaths
    channel_meas_paths = {}
    newcolnames = []
    for channel in channels:
        source_col = f'{channel}{source_col_suffix}'
        to_replace = f'{channel}_stack.tif'
        replace_with = f'{channel}_stack_measurement_rois.zip'
        newcolname = f'{channel}_stack_measurement_rois_path'
        newcolnames.append(newcolname)
        newvals = df[source_col].str.replace(to_replace, replace_with)
        channel_meas_paths[channel] = newvals
        df.loc[:, newcolname] = newvals

    # figure out which stack path measurement rois version has data,
    # if any. If neither refers to an existing rois.zip file, assume
    # the cell was segmented automatically
    df.loc[:, 'measurement_rois_path'] = np.nan
    for meas_path_name in newcolnames:
        # Set a generic "measurement_rois_path" columns using only generated measurement
        # roi paths that exist. If the path doesn't exist for any channel,
        # then the value will remain np.nan
        channel_bools = df.apply(lambda row : os.path.exists(row[meas_path_name]), axis=1)
        df.loc[channel_bools, 'measurement_rois_path'] = df.loc[channel_bools, meas_path_name]

def read_rectangular_rois_as_df(rois_path):
    """
    Use the read_roi package to read in the
    Fiji roi file at rois_path, store information
    from that roi file in a dataframe, and 
    return the dataframe
    """
    rois_dict = read_roi_zip(rois_path)
    # Create a dataframe using the crop rois file
    keys = [key for key in rois_dict.keys()]
    roi_dfs = []
    for key in keys:
        roi_dict = rois_dict[key]
        roi_df = pd.DataFrame(roi_dict, index = [key])
        roi_dfs.append(roi_df)

    rois_df = pd.concat(roi_dfs).reset_index(drop=True)
    # Add position_max so that we can which range of 
    # stack frame positions should use this ROI to find 
    # a cell etc.
    for index in rois_df.index:
        if index != np.max(rois_df.index):
            position_max = rois_df.loc[index+1, 'position'] - 1
        else:
            position_max = rois_df.loc[index, 'position']
            
        rois_df.loc[index, 'position_max'] = int(position_max)
        
    return rois_df

def read_roi_as_df(path):
    """
    Iterate through individual ROIs in the ROI file. Typically
    a stack of ROIs in time. Each ROI is a dictionary with
    a key that is the serial number of that ROI and a bunch of
    values which are also dicts
    """
    roi = read_roi_zip(path)
    # print(path)
    frame_dfs = []
    roi_index = 0
    for key, val in roi.items():
        # print(roi_index)
        # Cycle through the dictionaries for this ROI
        # and add the information to a dataframe (roi_df)
        frame_df = pd.DataFrame(columns=val.keys())
        val_types = [type(item) for item in val.values()]
        # print(val['type'])
        if list in val_types:
            for k, v in val.items():
            # Different dicts in the ROI have different dimensions
            # that need to be accounted for when population roi_df
                if type(v) == list and k == 'paths':
                    pass
                elif  type(v) == list and k != 'paths':
                    frame_df[k] = v
            for k, v in val.items():
                if type(v) == int or type(v) == str or type(v) == float:
                    frame_df.loc[:, k] = v

            frame_df.loc[:, 'roi_index'] = roi_index
            frame_df.loc[:, 'path'] = path
            frame_df.loc[:, 'frame'] = frame_df.position-1
            frame_dfs.append(frame_df)
            roi_index += 1

    roi_df = pd.concat(frame_dfs, ignore_index=True)
    return roi_df

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
        
        for key in channels_dict.keys():

            if channels_dict[key] in fov_channel_filename:

                src = os.path.join(fov_path, fov_channel_filename)
                dst = os.path.join(exptdir, f'{base_filename}_{fov.zfill(3)}_{key}.tif')
                shutil.copyfile(src, dst)

                found_channels.append(key)
    
    # Check if we found a file for every channel dictated in the master_index
    if len(found_channels) == len(channels_dict):
        pass
    else:
        print(f"Couldn't find a file for every channel in:\n{channels_dict.keys}. Only found files for the following: {found_channels}")


def get_channel_names(sampledir):
    """
    Return three strings joined on a space from three lists:
    fluor_channel_names, channel_names, raw_channel_names
    """
    metadict = get_byc_display_and_comments_dict(sampledir)
    namesdict = {'Brightfield': 'bf',
                 'YFP': 'yfp',
                 'RFP': 'rfp',
                 'GFP': 'gfp',
                 'mKO': 'mko',
                 'BFP': 'bfp'}
    fluor_channel_names = []
    channel_names = []
    raw_channel_names = []
    channeldictslist = metadict['Channels']
    for channeldict in channeldictslist:
        raw_name = channeldict['Name']
        name = namesdict[raw_name]
        raw_channel_names.append(raw_name)
        channel_names.append(name)
        if name != 'bf':
            fluor_channel_names.append(name)
            
    return ' '.join(fluor_channel_names), ' '.join(channel_names), ' '.join(raw_channel_names)

def make_ss_mdf(exptname, **kwargs):
    """
    Create and save master index made by scanning the directory
    matching `exptname` in constants.steady_state_data_dir
    and looking for features in those micromanager output
    directories defined in patterns
    """
    return_mdf = kwargs.get('return_mdf', True)
    save_mdf = kwargs.get('save_mdf', True)
    exptdir = os.path.join(constants.steady_state_data_dir, exptname)
    datadir = f'{exptdir}/data'
    if os.path.exists(datadir):
        print(f'Found dataset directory at \n{datadir}')
    else:
        print(f'No such file exists\n{datadir}')
        
    samplenames = os.listdir(datadir)
    # Filter non-directory items
    sampledirs = [os.path.join(datadir, name) for name in samplenames]
    sampledirs = [d for d in sampledirs if os.path.isdir(d)]
    # Extract features from sample names

    patterns = [constants.patterns.date,
                constants.patterns.strain_name,
                constants.patterns.plasmid_name,
                constants.patterns.tet_concn,
                constants.patterns.estradiol_concn,
                constants.patterns.genotype,
                constants.patterns.clone_number,
                constants.patterns.culture_condition,
                constants.patterns.timepoint_mins]
    patternkeys = ['date',
                   'strain',
                   'plasmid',
                   'tet_concn',
                   'estradiol_concn',
                   'genotype',
                   'clone',
                   'culture_condition',
                   'minutes']
    channelkeys = ['fluor_channel_names',
                   'channel_names',
                   'raw_channel_names']
    colnames = patternkeys + channelkeys
    patternsdict = dict(zip(patternkeys, patterns))
    # Create an empty master index to fill with sample information
    mdf = pd.DataFrame(index=range(len(sampledirs)), columns=colnames)
    # Iterate through each unique sample directory, detect features
    # from filename and add those features to the master index
    for i, sampledir in enumerate(sampledirs):
        # Add imaging channels collected to master index
        fluorchannels, channels, rawchannels = get_channel_names(sampledir)
        for idx, val in enumerate([fluorchannels, channels, rawchannels]):
            mdf.loc[i, channelkeys[idx]] = val
        samplename = os.path.basename(sampledir)
        matches = []
        groups = []
        for key, p in patternsdict.items():
            match = re.search(p, samplename)
            matches.append(match)
            if match is not None:
                # If we find a clone number, we only want the integer
                # of that clone, not the whole string
                if key == 'clone':
                    group = match.groups()[1]
                else:
                    group = match.group()
            else:
                # If we're looking for a plasmid_name and none are
                # found, we check if this is the 'no-plasmid' part
                # of the experiment and annotate as such
                if key == 'plasmid' and 'no-plasmid' in samplename:
                    group = 'no-plasmid'
                else:
                    group = np.nan
            mdf.loc[i, key] = group
        # Add measurement directory name to master index
        mdf.loc[i, 'measdirname'] = sampledir
    # Add path to the data directory for this
    # entire experiment
    mdf.loc[:, 'path'] = datadir
    if save_mdf:
        m = re.search(constants.patterns.date, exptname)
        date = m.group()
        savepath = os.path.join(exptdir, f'{date}_master_index.csv')
        mdf.to_csv(savepath, index=False)
    if return_mdf:
        return mdf

def rename_steady_state(master_index_df=None, **kwargs):
    """
    Ask the user to choose a master_index for their steady
    state experiment. Then rename files tifs based
    on that master index
    """
    if pd.DataFrame(master_index_df).empty == True:
        master_index_df = pd.read_csv(select_file("Choose master index .csv"))
    
    exptdir = master_index_df.path.iloc[0]
    conditions = os.listdir(exptdir)
    
    for condition_index in master_index_df.index:
        row = master_index_df.loc[condition_index, :]
        channel_dict = dict(zip(row.channel_names.split(), row.raw_channel_names.split()))
        conditionpath = os.path.join(row.path, row.measdirname)
        fov_paths = find_fov_paths(conditionpath)
        # Get rid of "_1" added to end of each sample name
        # by micro-manager during acquisition
        base_filename = row.measdirname
        if base_filename[-2:] == '_1':
            base_filename = base_filename[0:-2]
        if fov_paths != None:
            for fov_path in fov_paths:
                rename_channels(fov_path, channel_dict, base_filename, exptdir)
        else:
            print(f"No FOV directories found at conditionpath {conditionpath}")

def rename_steady_state_legacy(tet_treated=False, estradiol_treated=False, return_index=False):
    """
    Ask the user to choose a master_index for their steady
    state experiment. Then rename files tifs based
    on that master index

    condition_descriptor_cols should be a list of column names
    found in '<expt_date>_master_index.csv'

    Default conditions to look for in expt_dir are below
    where descriptor_list is set
    """

    # updating this so that a new master index df is instantiated
    # based on contents of the steady state exptdir
    master_index_df = pd.read_csv(select_file("Choose master index .csv"))
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
        if tet_treated:
            tet_concn = str(row.tet_concn).zfill(3) + 'uM-Tet'
            descriptor_list.append(tet_concn)
        if estradiol_treated:
            estradiol_concn = str(row.estradiol_concn).zfill(3) + 'nM-Estradiol'
            descriptor_list.append(estradiol_concn)
        # Check if there's a directory that matches the 
        # descriptors provided for this condition in 
        # the master index
        condition = find_condition(descriptor_list, conditions)
        if condition:
            print(f'Condition found: {condition}')
            base_filename = condition
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
    root.destroy() # important to destroy the root object, otherwise python 
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

def get_byc_compartmentdir(exptname, compartmentname, **kwargs):
    """
    Return path to byc experiment compartment directory given
    <exptname> and <compartmentname>
    
    Example:
    
    exptname = '20210903_byc'
    compartmentname = '20210903_JPC096_NLS-YFP-ODC(47)x3_BY4741_old_chase'
    dir = get_byc_compartmentdir(exptname, compartmentname)
    """
    return_exptdir = kwargs.get('return_exptdir', False)
    exptdir = os.path.join(constants.byc_data_dir, exptname)
    if os.path.exists(exptdir):
        pass
    else:
        print(f'{exptdir}\n does not exist, list of existing expt directories below\n')
        for name in os.listdir(constants.byc_data_dir):
            if 'byc' in name:
                print(name)
        return None
    compartmentdir = os.path.join(exptdir, compartmentname)
    if os.path.exists(compartmentdir):
        pass
    else:
        print(f'{compartmentdir} does not exist')
        for name in os.listdir(exptdir):
            print(name)
        return None    
    
    if return_exptdir:
        return exptdir, compartmentdir
    else:
        return compartmentdir

def mdf_from_file_pattern(compartmentdir, file_pattern, **kwargs):
    """
    Find all files that match <file_pattern>, typically
    drawn from the constants.patterns class, and make
    them into a concatenated DataFrame if they're all .csvs
    """
    save_mdf = kwargs.get('save_mdf', True)
    mdf_type = kwargs.get('mdf_type', 'crop_rois')
    return_savepath = kwargs.get('return_savepath', True)
    savepath = f'{compartmentdir}_{mdf_type}.csv'
    filenames = os.listdir(compartmentdir)
    df_fns = []
    for filename in filenames:
        m = re.search(file_pattern, filename)
        if m:
            df_fns.append(filename)
    filepaths = [os.path.join(compartmentdir, fn) for fn in df_fns]
    bools = ['.csv' == filepath[-4:] for filepath in filepaths]
    if False in bools:
        print(f'Found non .csv file type(s) in paths:')
        print(filepaths)
        return None
    else:
        pass
    print(f'Found {len(filepaths)} potential {mdf_type} roi df .csv files')
    # If there aren't any single cell index csvs of the type specified,
    # return an empty dataframe
    file_exist_bools = [os.path.exists(path) for path in filepaths]
    if True not in file_exist_bools:
        print(f'Found no csvs at {filepaths}')
        print(F'Returning empty dataframe')
        return pd.DataFrame(None)
    else:
        pass
    dfs = [pd.read_csv(path) for path in filepaths]
    # try:
    mdf = pd.concat(dfs, ignore_index=True).sort_values(by='cell_index')
    if 'crop_roi_set_path' in mdf.columns:
        print(f'Adding relative path to mdf')
        # Need to add relative path to crop roi set
        mdf.loc[:, 'crop_roi_set_relpath'] = mdf.crop_roi_set_path.apply(lambda x: str(x).replace(constants.byc_data_dir, ''))
    if 'bud_roi_set_path' in mdf.columns:
        print(f'Adding relative path to mdf')
        # Need to add relative path to crop roi set
        mdf.loc[:, 'bud_roi_set_relpath'] = mdf.bud_roi_set_path.apply(lambda x: str(x).replace(constants.byc_data_dir, ''))

    if 'compartment_dir' in mdf.columns:
        mdf.loc[:, 'compartment_reldir'] = mdf.compartment_dir.apply(lambda x: str(x).replace(constants.byc_data_dir, ''))
    if save_mdf:
        mdf.to_csv(savepath, index=False)
        print(f'Saved master index df at:\n{savepath}')
    if return_savepath:
        return (mdf, savepath)
    else:
        return mdf
    # except Exception as e:
    #     print(f'Could not concatanate .csvs into a single dataframe\nError: {e}')
    #     return None

def measurement_rois_path_from_crop_rois_path(cell_crop_rois_path, xy, channels):
    """
    Find and replace in the cell_crop_rois_path to define a path
    to the cell's measurement roi .zip file
    
    Return the path if exactly one is found, if more than one is 
    found return the last one found
    """
    measpath = None
    measpathtemp = cell_crop_rois_path.replace('byc_cell', f'byc_xy{str(xy).zfill(2)}_cell')
    measpathtemp = measpathtemp.replace('crop_rois', 'measurement_rois')

    for channel in channels:
        oldstr = '_measurement_rois'
        newstr = f'_{channel}_stack_measurement_rois'
        measpath = measpathtemp.replace(oldstr, newstr)
        if os.path.exists(measpath):
            print(f'Found measurement rois at {measpath}')
            measpath_final = measpath
        else:
            measpath_final = None
    
    return measpath_final

def path_annotate_master_index_df(mdf, **kwargs):
    """
    For each cell in the master index, add a absolute paths
    to the cell's crop rois set, bud rois set, source channel
    stack tifs, and measurement roi sets (cell outlines annotated
    on individual cell stacks)
    """
    channels = kwargs.get('channels', str(mdf.channels_collected.iloc[0]).split())
    print(f'Adding paths for channels {channels}')
    col_names = ['crop_roi_set_path',
                 'bud_roi_set_path',
                 'outline_roi_set_path']
    for channel in channels:
        col_names.append(f'{channel}_stack_path')
        
    # Create nan columns that will be filled with individual cell
    # information below
    for col_name in col_names:
        mdf.loc[:, col_name] = np.nan
    
    for cell_index in mdf.cell_index.unique():
        # Get paths to the cell's crop roi set, bud roi set,
        # and image used to annotate the crop roi
        cell_row = mdf.loc[mdf.cell_index==cell_index, :]

        crop_rois_relpath = cell_row.crop_roi_set_relpath.iloc[0]
        cell_crop_rois_path = os.path.join(constants.byc_data_dir, crop_rois_relpath)

        cell_bud_rois_path = cell_crop_rois_path.replace('crop_rois.zip', 'bud_rois.zip')

        xy = cell_row.xy.iloc[0]
        date = cell_row.date.iloc[0]
        expt_type = cell_row.expt_type.iloc[0]
        compdir = os.path.join(constants.byc_data_dir, cell_row.compartment_reldir.iloc[0])
        stacks_dict = {}
        
        print(f'Collecting paths to stacks for following channels:{channels}')
        for channel in channels: 
            channel_stack_filename = f'{date}_{expt_type}_xy{str(xy).zfill(2)}_{channel}_stack.tif'
            channel_stack_filepath = os.path.join(compdir, channel_stack_filename)
            stacks_dict[channel] = channel_stack_filepath
        # Get the path to the cell's measurement roi set
        args = [cell_crop_rois_path,
                xy,
                channels]
        outline_rois_path = measurement_rois_path_from_crop_rois_path(*args)
        # Add information found above to the master index
        vals = [cell_crop_rois_path,
                cell_bud_rois_path,
                outline_rois_path]
        # Add column name and path to stack tif
        # for each channel collected. This is specific
        # to each cell because they can come from different
        # position (xy) stakcs
        for key in stacks_dict.keys():
            vals.append(stacks_dict[key])
        for i, col_name in enumerate(col_names):
            mdf.loc[mdf.cell_index==cell_index, col_name] = vals[i]
        
    return mdf

def get_likely_compartment_dirs(**kwargs):
    """
    Return a list of full compartment directory paths found
    in byc.constants.byc_data_dir that match the 
    'byc.constants.byc_data_dir' pattern. 
    """
    # the pathlib module Path allows creation of a Path object with 
    # a parents attribute. The 2nd parent up should be byc_data_dir
    # if the fits_table is in a compartment dir
    rootdir = os.path.abspath(constants.byc_data_dir)
    dirs = [x[0] for x in os.walk(rootdir)]
    dirs = [d for d in dirs if os.path.abspath(str(Path(d).parents[1])) == rootdir]
    
    return dirs

def get_fits_table_paths(**kwargs):
    """
    Return a list of paths to .csv files containing
    exponential fitting, distance from senescence, and
    other analysis with a row for each cell in an individual
    flow compartment of a BYC experiment
    """
    dirs = get_likely_compartment_dirs()
    ft_paths = [os.path.join(d, 'fits_table.csv') for d in dirs]
    ft_paths = [p for p in ft_paths if os.path.exists(p)]
    
    return ft_paths
    
def get_allfitsdf_paths(**kwargs):
    """
    Return a list of paths to .csv files containing
    exponential fitting, distance from senescence, full
    fluorescence traces, and other analysis with a row for
    each timepoint per cell in an individual
    flow compartment of a BYC experiment
    """
    dirs = get_likely_compartment_dirs()
    df_paths = [os.path.join(d, 'allfitsdf.csv') for d in dirs]
    df_paths = [p for p in df_paths if os.path.exists(p)]
    
    return df_paths