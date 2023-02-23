import os
import re

import tkinter as tk
import tkinter.filedialog as dia

import pandas as pd
import numpy as np

from byc import constants, files

def set_fp(prompt):

    """ Return the path to a file of interest. Call with the prompt you would like to 
        display in the dialog box."""

    # create the dialog box and set the fn
    root = tk.Tk()
    fp = dia.askopenfilename(parent=root, title=prompt)
    root.destroy() # very important to destroy the root object, otherwise python 
    # just keeps running in a dialog box

    return fp # return the path to the file you just selected

def filename_from_path(path):
    """
    If path points to a file,
    return the name of the file
    without upstream path
    """
    if os.path.isfile(path):
        filename = path[path.rindex('\\')+1:]
        return filename
    elif os.path.isdir(path):
        print(f"Path:\n{path}\nis a directory, note a file")
        return None

def dirname_from_dirpath(path):
    """
    If path points to a directory,
    return the name of the directory
    without upstream path
    """
    if os.path.isdir(path):
        dirname = path[path.rindex('\\')+1:]
        return dirname
    elif os.path.isfile(path):
        print(f"Path:\n{path}\nis a file, note a directory")
        return None

def get_master_index_paths(exptdir):
    """
    Look for master_index.csv files in the
    exptdir
    """    
    master_index_paths = []
    if os.path.isdir(exptdir):
        files = os.listdir(exptdir)
        for filename in files:
            if filename[-4:] == '.csv' and 'master_index' in filename:
                filepath = os.path.join(exptdir, filename)
                master_index_paths.append(filepath)
            else:
                pass
    else:
        print(f"exptdir:\n{exptdir}\nis not a directory")
        return None

    if len(master_index_paths) > 0:
        return master_index_paths
    else:
        print(f"Could not find any master_index.csv type files in exptdir:\n{exptdir}")
        return None

def make_blank_master_index(exptdir, date, expttype='byc', write=True):
    """
    Create a .csv from a standard, empty
    master_index_df
    """
    if expttype == 'byc':
        df = pd.DataFrame(columns=constants.master_index_cols)
        writepath = os.path.join(exptdir, f'{date}_master_index.csv')

    elif expttype == 'steadystate':
        df = pd.DataFrame(columns=constants.master_index_cols)
        writepath = os.path.join(exptdir, f'{date}_master_index.csv')

    if write:        
        df.to_csv(writepath, index=False)
        print(f"Saved blank master_index at:\n{writepath}")

    return df, writepath

def get_relpath(abspath, **kwargs):
    """
    Return the path of <abspath> relative to the byc data
    directory (defined in constants.byc_data_dir)

    Right now this function just looks for the word "data"
    in the path and the relative path is everything after that

    Was having issues with matching actual paths to identify
    which part of string is the byc_data_dir between OSs
    """
    # upstream part of path should either be in legacy
    # directory or in <byc source path>\data\<relpath>
    keyword = kwargs.get('include_everything_after_string', 'data')
    if keyword in abspath:
        relpath = abspath[abspath.rindex(keyword) + len(keyword)+1:]
    else:
        print(f'{keyword} not found in path:\n{abspath}')
        relpath = None
    return relpath

def get_master_index_tags(match):
    """
    Return a list of words in the filename
    that are not "master_index" or ".csv"
    """
    tags = []

    if match != None:
        modifiers = [match.group(1), match.group(3)]

        for mod in modifiers:
            split_mod = mod.split('_')

            for sm in split_mod:
                if sm != '':
                    tags.append(sm)

    if len(tags) == 0:
        print(f'No tags found in string: {match.string}')

    return tags

def get_all_master_index_paths(rootdir=constants.byc_data_dir, get_tags=False):
    """
    Return the list of paths to all files that match the
    master index pattern in constants.byc_data_dir (recursive
    search)
    """
    pattern = constants.patterns.master_index_file
    tags_lists = []
    filepaths = []

    for dirpath, dirnames, filenames in os.walk(rootdir):   

        for filename in filenames:
            match = re.search(pattern, filename)

            if match:
                filepath = os.path.join(dirpath, filename)

                if os.path.exists(filepath):
                    filepaths.append(filepath)
                    tags= get_master_index_tags(match)
                    tags_lists.append(tags)

    if len(filepaths) == 0:
        print(f"No master index paths found in rootdir:\n{rootdir}")

    if get_tags:
        return filepaths, tags_lists
    else:
        return filepaths

def get_filepaths_with_pattern(pattern, rootdir=constants.byc_data_dir):

    filepaths = []

    for dirpath, dirnames, filenames in os.walk(rootdir):   

        for filename in filenames:
            match = re.search(pattern, filename)

            if match:
                filepath = os.path.join(dirpath, filename)

                if os.path.exists(filepath):
                    filepaths.append(filepath)

    return filepaths

def get_all_master_index_dfs(**kwargs):
    """
    Return a list of dataframes made by reading
    csvs at each path returned by get_all_master_index_paths
    """
    master_index_paths = kwargs.get('master_index_paths',
                                    get_all_master_index_paths())

    dfs = [pd.read_csv(path) for path in master_index_paths if pd.read_csv(path).empty == False]
    return dfs

def exptname_from_compartment(compartmentdir_or_name):
    """
    Find experiment name (e.g. '20220610_byc') in <compartmentdir_or_name>
    as long as <compartmentdir_or_name> is in format:
    <DATE>_byc_<expt_details_etc>
    """
    splitname = compartmentdir_or_name.split('_')
    exptname = '_'.join(splitname[0:2])

    return exptname

def get_cell_indices_in_compartment(compartmentdir, **kwargs):
    filenames = os.listdir(compartmentdir)
    indices = '|'.join([str(num).zfill(3) for num in range(1000)])
    cellindexpattern = f"cell({indices})"

    cellindexpattern = kwargs.get('cellindexpattern', cellindexpattern)

    cell_indices = [re.search(cellindexpattern, fn).groups()[0] for fn in filenames if re.search(cellindexpattern, fn)]
    cell_indices = np.array([np.int32(num) for num in cell_indices])
    cell_indices = np.unique(cell_indices)

    return cell_indices

def annotate_channel_paths_in_mdf(mdf, **kwargs):

    channel_stack_kw = kwargs.get('channel_stack_kw', '_stack.tif')
    indices = '|'.join([str(num).zfill(3) for num in range(1000)])
    cellindexpattern = f"cell({indices})"
    indices = '|'.join([str(num).zfill(2) for num in range(100)])
    fovpattern = f"xy({indices})"
    possible_channels = [
    'bf',
    'bfp',
    'gfp',
    'yfp',
    'mko',
    'rfp'
    ]
    groups = '|'.join(possible_channels)
    channelpattern = f'({groups})'    
    
    compdir = mdf.loc[0, 'compartment_dir']
    allfilenames = os.listdir(compdir)
    for fn in allfilenames:
        channel_stack_conditions = [
            fn[-len(channel_stack_kw):] == channel_stack_kw,
            re.search(cellindexpattern, fn),
            re.search(fovpattern, fn)
        ]
        if channel_stack_conditions[0] and channel_stack_conditions[1] and channel_stack_conditions[2]:
            # This is a cell crop channel stack
            cell_index = np.int32(re.search(cellindexpattern, fn).groups()[0])
            fov = np.int32(re.search(fovpattern, fn).groups()[0])
            channel_name = re.search(channelpattern, fn).group()
            channel_stack_path = os.path.join(compdir, fn)
            channel_stack_relpath = get_relpath(channel_stack_path)
            mdf.loc[cell_index, f'{channel_name}_crop_stack_path'] = channel_stack_path
            mdf.loc[cell_index, f'{channel_name}_crop_stack_relpath'] = channel_stack_relpath

def annotate_bud_and_crop_paths_in_mdf(mdf, **kwargs):

    compdir = mdf.loc[0, 'compartment_dir']
    exptname = mdf.loc[0, 'exptname']

    cell_indices_strs = [val.zfill(3) for val in mdf.index.astype(str)]

    # Set bud roi path location information
    bud_roi_fns = [f'{exptname}_cell{cell_index_str}_bud_rois.zip' for cell_index_str in cell_indices_strs]
    bud_roi_paths = [os.path.join(compdir, fn) for fn in bud_roi_fns]
    bud_roi_relpaths = [get_relpath(path) for path in bud_roi_paths]
    mdf.loc[:, 'bud_roi_set_path'] = bud_roi_paths
    mdf.loc[:, 'bud_roi_set_relpath'] = bud_roi_relpaths
    # Set crop roi path location information
    crop_roi_fns = [f'{exptname}_cell{cell_index_str}_crop_rois.zip' for cell_index_str in cell_indices_strs]
    crop_roi_paths = [os.path.join(compdir, fn) for fn in crop_roi_fns]
    crop_roi_relpaths = [get_relpath(path) for path in crop_roi_paths]
    mdf.loc[:, 'crop_roi_set_path'] = crop_roi_paths
    mdf.loc[:, 'crop_roi_set_relpath'] = crop_roi_relpaths

def annotate_bud_and_crop_df_info_in_mdf(mdf, **kwargs):
    """
    bud and crop df .csvs contain information entered during data annotation
    and roi .zip file saving using imagejpc
    """

    compdir = mdf.loc[0, 'compartment_dir']
    exptname = mdf.loc[0, 'exptname']

    cell_indices_strs = [val.zfill(3) for val in mdf.index.astype(str)]

    # Set bud roi path location information
    bud_roi_fns = [f'{exptname}_cell{cell_index_str}_bud_rois.zip' for cell_index_str in cell_indices_strs]
    bud_roi_paths = [os.path.join(compdir, fn) for fn in bud_roi_fns]
    bud_roi_relpaths = [get_relpath(path) for path in bud_roi_paths]
    mdf.loc[:, 'bud_roi_set_path'] = bud_roi_paths
    mdf.loc[:, 'bud_roi_set_relpath'] = bud_roi_relpaths
    # Set crop roi path location information
    crop_roi_fns = [f'{exptname}_cell{cell_index_str}_crop_rois.zip' for cell_index_str in cell_indices_strs]
    crop_roi_paths = [os.path.join(compdir, fn) for fn in crop_roi_fns]
    crop_roi_relpaths = [get_relpath(path) for path in crop_roi_paths]
    mdf.loc[:, 'crop_roi_set_path'] = crop_roi_paths
    mdf.loc[:, 'crop_roi_set_relpath'] = crop_roi_relpaths

def annotate_bud_and_crop_df_info_in_mdf(mdf, **kwargs):
    """
    bud and crop df .csvs contain information entered during data annotation
    and roi .zip file saving using imagejpc

    Return nothing
    """
    compartment_reldir = get_relpath(mdf.loc[0, 'compartment_dir'])
    compdir = os.path.join(constants.byc_data_dir, compartment_reldir)
    filenames = os.listdir(compdir)
    bud_rois_df_kw = '_bud_rois_df.csv'
    crop_rois_df_kw = '_crop_rois_df.csv'
    crop_rois_df_paths = [os.path.join(compdir, fn) for fn in filenames if crop_rois_df_kw in fn]
    bud_rois_df_paths = [os.path.join(compdir, fn) for fn in filenames if bud_rois_df_kw in fn]
    crop_rois_df_relpaths = [get_relpath(path) for path in crop_rois_df_paths]
    bud_rois_df_relpaths = [get_relpath(path) for path in bud_rois_df_paths]
    mdf.loc[:, 'crop_rois_df_path'] = crop_rois_df_paths
    mdf.loc[:, 'bud_rois_df_path'] = bud_rois_df_paths
    mdf.loc[:, 'crop_rois_df_relpath'] = crop_rois_df_relpaths
    mdf.loc[:, 'bud_rois_df_relpath'] = bud_rois_df_relpaths

    crop_rois_dfs = [pd.read_csv(path) for path in crop_rois_df_paths]
    bud_rois_dfs = [pd.read_csv(path) for path in bud_rois_df_paths]

    crop_rois_alldf = pd.concat(crop_rois_dfs).reset_index()
    bud_rois_alldf = pd.concat(bud_rois_dfs).reset_index()
    # End event type is only accurately annotated in the
    # bud_rois_df so don't annotate it here
    for col in crop_rois_alldf.columns:
        if col not in mdf.columns and col != 'end_event_type':
            mdf.loc[:, col] = crop_rois_alldf.loc[:, col]
    # Contains aggregate is only accurately annotated in the 
    # crop_rois_df so don't annotate it here
    for col in bud_rois_alldf.columns:
        if col not in mdf.columns and col != 'contains_aggregate':
            mdf.loc[:, col] = bud_rois_alldf.loc[:, col]



def generate_mdf(exptname, compartmentname, **kwargs):
    """
    2nd generation master index dataframe generation
    Scan the directory matching <compartmentname> for 
    files corresponding to individual cell stacks, bud
    rois, and crop rois etc.

    Return the generated master_index_df (mdf)
    """
    channels = kwargs.get('channels', ['bf', 'yfp', 'rfp'])
    compdir = files.get_byc_compartmentdir(exptname, compartmentname)
    compdir = os.path.abspath(compdir)
    cell_indices = get_cell_indices_in_compartment(compdir)
    mdf = pd.DataFrame({
        'cell_index': cell_indices
    })
    mdf.loc[:, 'channels_collected'] = ' '.join(channels)
    mdf.loc[:, 'compartment_dir'] = compdir
    mdf.loc[:, 'compartment_reldir'] = get_relpath(compdir)
    mdf.loc[:, 'exptname'] = exptname
    mdf.loc[:, 'date'] = mdf.exptname.str[0:8]
    mdf.loc[:, 'compartment_name'] = mdf.compartment_dir.apply(lambda x: os.path.basename(x))
    mdf.loc[:, 'exptname'] = exptname
    mdf.set_index('cell_index', inplace=True)
    annotate_bud_and_crop_paths_in_mdf(mdf)
    annotate_bud_and_crop_df_info_in_mdf(mdf)
    print(f'Adding paths for channels {channels}')
    files.add_cell_channel_crop_stack_paths(mdf, channels=channels)
    files.add_cell_channel_xy_source_stack_paths(mdf, channels=channels)

    return mdf

def check_bud_and_crop_roi_dfs(compartmentdir):
    allfns = os.listdir(compartmentdir)
    crop_zips = [fn for fn in allfns if fn[-4:]=='.zip' and 'crop' in fn]

    crop_dfs = [fn for fn in allfns if re.search(constants.patterns.crop_roi_df_file, fn) != None]
    buds_dfs = [fn for fn in allfns if re.search(constants.patterns.bud_roi_df_file, fn) != None]

    def extract_cell_index(string):
        rdx = string.rindex('cell')
        number = string[rdx + 4: rdx + 7]

        return int(number)


    crop_cell_idxs = [extract_cell_index(fn) for fn in crop_dfs]
    bud_cell_idxs = [extract_cell_index(fn) for fn in buds_dfs]

    all_annotated_cell_indices = [extract_cell_index(fn) for fn in crop_zips]

    missed_crop_dfs = [idx for idx in all_annotated_cell_indices if idx not in crop_cell_idxs]
    missed_bud_dfs = [idx for idx in all_annotated_cell_indices if idx not in bud_cell_idxs]


    print(f'Missing crop df .csv files for cells {missed_crop_dfs}')
    print(f'Missing bud df .csv files for cells {missed_bud_dfs}')

def check_for_crop_stack_exist(mdf):
    """
    Return False if any of the crop_stack_paths in the <mdf>
    do not exist
    """
    all_channel_crop_stack_exist_bools = []
    for channel in mdf.channels_collected.iloc[0].split(' '):
        pathcolname = f'{channel}_crop_stack_path'
        exist_bools = [os.path.exists(path) for path in mdf[pathcolname]]
        all_channel_crop_stack_exist_bools = all_channel_crop_stack_exist_bools + exist_bools

    if False in all_channel_crop_stack_exist_bools:
        return False
    else:
        return True