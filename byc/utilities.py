import os
import re

import tkinter as tk
import tkinter.filedialog as dia

import pandas as pd

from byc import constants

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

def get_relpath(abspath):
    """
    Return the relative portion of the abspath
    by replacing expected upstream abs path
    (constants.byc_dat_dir etc.) with ''.
    """
    # upstream part of path should either be in legacy
    # directory or in <byc source path>\data\<relpath>
    abspath = os.path.abspath(abspath)
    possibles = [constants.legacy_byc_data_dir,
                 constants.byc_data_dir]
    relpath = None
    for query in possibles:
        #print(f"Checking:\n{abspath}\nfor:\n{query}")
        if abspath.find(query) == 0:
            relpath = abspath.replace(query, '')
        else:
            pass
        
    if relpath != None:
        newabspath = os.path.join(constants.byc_data_dir, relpath)
        if os.path.exists(newabspath):
            return relpath
        else:
            print(f"constants.byc_data_dir + found relpath:\n{newabspath}\ndoes note exist!")
            return relpath
    else:
        print(f"No match found after looking for:\b{possibles}\nin:\n{abspath}")
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

    for dirpath, dirnames, filenames in os.walk(constants.byc_data_dir):   

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

def get_all_master_index_dfs(**kwargs):
    """
    Return a list of dataframes made by reading
    csvs at each path returned by get_all_master_index_paths
    """
    master_index_paths = kwargs.get('master_index_paths',
                                    get_all_master_index_paths())

    dfs = [pd.read_csv(path) for path in master_index_paths if pd.read_csv(path).empty == False]
    return dfs



