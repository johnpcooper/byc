import os

import tkinter as tk
import tkinter.filedialog as dia

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
        print(f"Couldn't find any master_index.csv type files in exptdir:\n{exptdir}")