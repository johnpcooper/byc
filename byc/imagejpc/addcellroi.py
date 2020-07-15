"""
This script gets run by imagejpc/macros/addCell.ijm
as part of the plugini imagejpc/utilities/save_cell_roi_set.py"

Might also add byc.database.dataBase integration at some point
"""
import argparse
import os
import pandas as pd
import sys

from byc import constants, database, utilities

def find_cell_row(cell_roi_df, master_index_df):
    """
    Search the master index for a row that matches the cell_index,
    xy, and date contained in the cell_roi_df and return that row
    
    If no matching index found, return None
    """
    cell_roi_cols = ['cell_index', 'xy', 'date']
    
    if 'sub_coord' in master_index_df.columns:
        master_index_cols = ['sub_coord', 'xy', 'date']
    elif 'cell_index' in master_index_df.columns:
        master_index_cols = ['cell_index', 'xy', 'date']
    else:
        print(f"No cell_index or sub_coord column in master_index_df")
        return None

    matching_rows = []
    for row in master_index_df.index:

        query = [int(cell_roi_df.loc[0, colname]) for colname in cell_roi_cols]
        master = [int(master_index_df.loc[row, colname]) for colname in master_index_cols]

        if query == master:
            matching_rows.append(row)

    if len(matching_rows) > 1:
        print(f"Found multiple rows matching cell_roi_df. Check master index for duplicates at:/n{master_index.path.iloc[0]}")
        return None
    elif len(matching_rows) == 0:
        print("Found no rows in master_index_df matching cell_roi_df")
        return None
    else:
        matching_row = matching_rows[0]
        return matching_row

def merge_cell_to_index(cell_roi_df, master_index_df, path, row, write=True):
    """
    Return the master_index_df with added cell_roi_df
    data. Overwrites values that are already in 
    existing master index for common columns
    """
    cols_to_add = cell_roi_df.columns

    for col in cols_to_add:
        master_index_df.loc[row, col] = cell_roi_df.loc[0, col]
    
    if write:        
        master_index_df.to_csv(path, index=False)
        print(f"Added cell info to master index and saved at\n{path}")
    else:
        return master_index_df

def get_exptdir(active_imp_path):

    exptdir = os.path.abspath(os.path.join(os.path.dirname(active_imp_path), '..'))
    if os.path.exists(exptdir):
        return exptdir
    else:
        print(f"exptdir for active_imp_path:\n{active_imp_path}\ndoes not exist at\n{exptdir}")
        return None

def find_cell_master_index(cell_roi_df, active_imp_path):
    """
    Return the master_index_df and its path that has
    a single row matching the cell_index, xy, and date
    values in cell_roi_df
    """
    # Get the path to the directory holding the directory
    # holding the active imp
    exptdir = get_exptdir(active_imp_path)
    # Get all files that look like master indexes in that
    # exptdir (e.g. data/20200221_byc/20200221_byc_master_index.csv)
    master_index_paths = utilities.get_master_index_paths(exptdir)
    matched_paths = []
    matched_dfs = []
    matched_rows = []
    if master_index_paths != None:
        for i, master_index_df in enumerate([pd.read_csv(path) for path in master_index_paths]):
            # Scan through the master dfs found to look for 
            # ones that have a row that matches 
            row = find_cell_row(cell_roi_df, master_index_df)        
            if row != None and type(row) == int:
                matched_paths.append(master_index_paths[i])
                matched_dfs.append(master_index_df)
                matched_rows.append(row)
            else:
                pass
    else:
        print(f"No master_index_dfs found for cell_roi_df in path:\n{exptdir}")
        return None, None, None

    assert len(matched_dfs) == len(matched_dfs) == len(matched_rows), "paths and dfs list lengths do not match"
    if len(matched_dfs) == 1:
        print(f"Found cell match at row {matched_rows[0]} in master index:\n{matched_paths[0]}")
        return matched_dfs[0], matched_paths[0], matched_rows[0]
    elif len(matched_dfs) == 0:
        print(f"No master_index_dfs found for cell_roi_df in path:\n{exptdir}")
        return None, None, None
    elif len(matched_dfs) > 1:
        print(f"Multiple master_index_dfs found for cell_roi_df in path:\n{exptdir}")
        return None, None, None
    else:
        print("Unknown error")
        return None, None, None

def record_cell_roi_set(cell_index,
                      active_imp_path,
                      roi_set_save_path,
                      end_event_type,
                      roi_set_type,
                      active_imp_dir,
                      **kwargs):
    
    write = kwargs.get('write', True)
    imp_filename = utilities.filename_from_path(active_imp_path)
    exptdir = get_exptdir(active_imp_path)
    # Extract some information from the active image filename
    xy = int(imp_filename[imp_filename.rindex('xy') + 2: imp_filename.rindex('xy') + 4])
    date = imp_filename[0:8]

    # Create path names relative to byc source/data directory
    
    active_imp_rel_path = utilities.get_relpath(active_imp_path)
    roi_set_rel_path = utilities.get_relpath(roi_set_save_path)
    if active_imp_rel_path == None and roi_set_rel_path == None:
        print(f"Could not find byc_data_dir ({constants.byc_data_dir}) in path:/n{active_imp_rel_path}")
        print(f"\nMake sure you're analyzing data in byc_data_dir")
    else:
       pass


    values = [cell_index,
              active_imp_path,
              roi_set_save_path,
              end_event_type,
              xy,
              date,
              active_imp_rel_path,
              roi_set_rel_path]

    keys = ['cell_index',
            'active_imp_abs_path',
            '{}_roi_set_abs_path'.format(roi_set_type),
            'end_event_type',
            'xy',
            'date',
            '{}_active_imp_rel_path'.format(roi_set_type),
            '{}_roi_set_rel_path'.format(roi_set_type)]

    cell_roi_df = pd.DataFrame(dict(zip(keys, values)), index=[0])
    # Find the master_index.csv df for this cell and add it to the 
    # master index.
    master_index_df, path, row = find_cell_master_index(cell_roi_df, active_imp_path)
    mdf = pd.DataFrame(master_index_df)
    # If there's a master index with an existing row for this cell
    # then add data to that row
    if mdf.empty == False and path != None and row != None:
        merge_cell_to_index(cell_roi_df, master_index_df, path, row)
    # If there's a master index in exptdir but not a row for this cell
    # then add a new row for that cell with its data
    elif mdf.empty == False and path != None and row == None:
        row = mdf.index.max() + 1
        merge_cell_to_index(cell_roi_df, master_index_df, path, row)
        print("Appended cell to existing master index, user needs to update other info")
    # If there's no dataframe and no row found, create a blank one and
    # merge the cell's data with it
    elif mdf.empty == True and path == None and row == None:
        # Create a new master, add the cell to it and save
        master_index_df, path = utilities.make_blank_master_index(exptdir, date)
        row = 0
        merge_cell_to_index(cell_roi_df, master_index_df, path, row)
    if write:
        fn = f'{date}_byc_xy{str(xy).zfill(2)}cell{str(cell_index).zfill(3)}_{roi_set_type}_rois_df.csv'
        writepath = os.path.join(active_imp_dir, fn)
        cell_roi_df.to_csv(writepath, index=False)
        print("Saved cell_roi_df at {}".format(writepath))
        return cell_roi_df
    else:
        return cell_roi_df

if __name__ == "__main__":

    if len(sys.argv) == 7:
        record_cell_roi_set(sys.argv[1],
                              sys.argv[2],
                              sys.argv[3],
                              sys.argv[4],
                              sys.argv[5],
                              sys.argv[6])
    else:
        print("wrong number of arg variables. Should be 6")