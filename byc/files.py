import os
import re
import shutil
import glob

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
    print(path)
    frame_dfs = []
    roi_index = 0
    for key, val in roi.items():
        print(roi_index)
        # Cycle through the dictionaries for this ROI
        # and add the information to a dataframe (roi_df)
        frame_df = pd.DataFrame(columns=val.keys())
        val_types = [type(item) for item in val.values()]
        print(val['type'])
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

def find_measdirname_features(measdirpath, feature_patterns):
    """
    Take a micromanager multi-d acquisition output directory, e.g.
    "20210421_pJC031_BY4741_constantOD_250uM-Tet_time000_1", find
    values for different feature categories defined in feature_patterns,
    and return a dictionary of those feature names and their regex
    match objects in a dictionary
    """
    measdirname = os.path.basename(os.path.dirname(measdirpath))
    matches_dict = {}
    sampledf = pd.DataFrame({'sample': measdirname}, index=[0])
    for key, val in feature_patterns.items():
        feature_match = re.search(val, measdirname)
        if feature_match:
            matches_dict[key] = feature_match
    # Add filename information
    matches_dict['path'] = measdirpath.replace(measdirname, '')
    matches_dict['measdirname'] = measdirname
    # Add default information about which channels were collected
    matches_dict['fluor_channel_names'] = ' '.join(constants.default_fluor_channels)
    matches_dict['channel_names'] = ' '.join(constants.default_channel_names)
    matches_dict['raw_channel_names'] = ' '.join(constants.default_raw_channel_names)
    # print(measdirpath)
    # print(matches_dict)
    return matches_dict

def measdf_from_features_dict(features_dict):
    """
    Iterate through each regex match object in the features
    dict, if the feature needs to be quantified, extract the 
    number. If not, extract only the match.group(). Then put
    these feature names and values into a single row dataframe
    for the measurement
    
    Return the dataframe
    """

    for key, val in features_dict.items():
        group0_vars = ['tet_concn',
                       'estradiol_concn']
        group1_vars = ['minutes',
                       'clone_number']
        # Some information annotated in find_measdirname_features
        # is already a string etc., so only extract groups if
        # the value in the features_dict is actually an re.Match
        if type(val) == re.Match:
            if key in group1_vars:
                number = val.groups()[1]
                features_dict[key] = np.float(number)
            elif key in group0_vars:
                number = val.groups()[0]
                features_dict[key] = np.float(number)
            else:
                value = val.group()
                features_dict[key] = value

    measdf = pd.DataFrame(features_dict, index=[0])
    print(measdf)
    return measdf

def make_ss_mdf(exptname, **kwargs):
    """
    Create and save master index made by scanning the directory
    matching `exptname` in constants.steady_state_data_dir
    and looking for features in those micromanager output
    directories defined in featuer_patterns
    """
    write_mdf = kwargs.get('write_mdf', True)
    return_mdf = kwargs.get('return_mdf', True)
    ssdir = constants.steady_state_data_dir
    exptdir = os.path.join(ssdir, exptname)
    datadir = os.path.join(exptdir, 'data')
    measdirpaths = glob.glob(f"{datadir}/*/")
    measdirpaths = [p for p in measdirpaths if os.path.isdir(p)]
    # print(f'Looking for features in the following paths\n{measdirpaths}')
    patterns = constants.patterns

    feature_patterns = {'expt_date': patterns.date,
                        'plasmid': patterns.plasmid_name,
                        'genotype': patterns.genotype,
                        'tet_concn': patterns.tet_concn,
                        'estradiol_concn': patterns.estradiol_concn,
                        'culture_condition': patterns.culture_condition,
                        'minutes': patterns.timepoint_mins,
                        'clone_number': patterns.clone_number}

    # Make indivudal rows of the master index (each row a different measurement)
    measdfs = []
    for measdirpath in measdirpaths:
        features_dict = find_measdirname_features(measdirpath, feature_patterns)
        measdf = measdf_from_features_dict(features_dict)
        measdfs.append(measdf)
    # Create the master index
    mdf = pd.concat(measdfs, ignore_index=True)
    
    if write_mdf:
        filename = f'{mdf.expt_date.iloc[0]}_master_index.csv'
        writepath = os.path.join(exptdir, filename)
        mdf.to_csv(writepath, index=False)
        print(f'Saved master index at \n{writepath}')
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
        base_filename = row.measdirname.replace('_1', '')
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
        print(f'<exptdir> does not exist, list of existing expt directories below\n')
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
    
    dfs = [pd.read_csv(path) for path in filepaths]
    try:
        mdf = pd.concat(dfs, ignore_index=True).sort_values(by='cell_index')
        if save_mdf:
            mdf.to_csv(savepath, index=False)
            print(f'Saved master index df at:\n{savepath}')
        if return_savepath:
            return (mdf, savepath)
        else:
            return mdf
    except Exception as e:
        print(f'Could not concatanate .csvs into a single dataframe\nError: {e}')
        return None

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
            pass
    
    return measpath_final

def path_annotate_master_index_df(mdf):
    """
    For each cell in the master index, add a absolute paths
    to the cell's crop rois set, bud rois set, source channel
    stack tifs, and measurement roi sets (cell outlines annotated
    on individual cell stacks)
    """
    channels = str(mdf.channels_collected.iloc[0]).split()
    col_names = ['crop_rois_path',
                 'bud_rois_path',
                 'outline_rois_path']
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
        # Create nan columns for the info to be added
        for i, col_name in enumerate(col_names):
            mdf.loc[mdf.cell_index==cell_index, col_name] = vals[i]
        
    return mdf    