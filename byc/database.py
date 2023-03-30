import re
from pathlib import Path

import os
import pandas as pd
import numpy as np

from byc import constants, utilities, files, plasmids, trace_tools
from byc import standard_analysis as sa
from byc.constants import patterns

class DataBase():
    """
    This class is the portal through which all data
    in this project can be accessed. Ultimately, 
    constants.byc_data_dir will contain flow, steady
    state imaging, and other kinds of data sets. 

    Each experiment directory (eg. '20200320_byc') will 
    contain a master index file (eg. '20200320_byc_master_index.csv').
    This file is used by database.DataBase to find and annotate
    data in that experiment directory
    """
    _byc_data_dir = constants.byc_data_dir
    # Going to get rid of ss and fc data dirs, put every expt. directory
    # in constants.byc_data_dir, and then distinguish expt_type by 
    # looking in master index in primary subdirectories.    
    _byc_trace_database_path = constants.database_paths['cell_trace_db_path']

    def __init__(self):

        expt_names = os.listdir(self.byc_data_dir)
        self.trace_database_path = self._byc_trace_database_path
        self.trace_database_df = pd.read_csv(self.trace_database_path)
        self.master_index_dfs_dict = self.get_master_index_dfs_dict()

    @property
    def byc_data_dir(self):
        return self._byc_data_dir

    @byc_data_dir.setter
    def byc_data_dir(self, byc_data_dir):
        self._byc_data_dir = byc_data_dir

    @property
    def trace_database_df(self):
        return self._trace_database_df

    @trace_database_df.setter
    def trace_database_df(self, trace_database_df):
        # Can set constraints here for making sure 
        # self.trace_database_df is set appropriately
        self._trace_database_df = trace_database_df

    def get_master_index_dfs_dict(self):
        """ 
        Return a dictionary where each value is
        a dataframe created from each master index
        file identified in constants.byc_data_dir
        and each key is an identifier based on the
        name of that master index.

        For a typical byc experiment, I create a
        separate master index for each flow
        compartment imaged. This is also how the
        imagejpc.addcellroi script works

        This will serve as the foundation for a more final
        master index database, where each row is a unique
        master index with its relpath recorded and other 
        information extracted from the master index df
        itself recorded in its columns (say plasmid, expr_type
        tet_perfusion_frames, etc.)
        """
        paths, tags_list = utilities.get_all_master_index_paths(get_tags=True)
        dfs = utilities.get_all_master_index_dfs()
        keys = ['_'.join(tags) for tags in tags_list]

        master_index_dfs_dict = dict(zip(keys, dfs))
        return master_index_dfs_dict

    def initialize_cell_trace_database(self, overwrite=True):
        """
        Completely clear the database .csv file and save a blank
        one with columns passed here. Reset self.trace_database_df 
        to the initialized df. 
        """
        # 
        cols = ['expt_name', 'trace_path', 'trace_relpath',
                'chase_index', 'date']
        trace_db_df = pd.DataFrame(columns = cols)
        if overwrite:
            writepath = self._byc_trace_database_path
        else:
            writepath = self._byc_trace_database_path

        trace_db_df.to_csv(write_path, index=False)
        self.trace_database_df = pd.read_csv(writepath)
        
    def add_expt_to_trace_database(self, cell_trace_df_paths, expt_name, chase_index, date):
        # Use the list of file names chosen by the user
        # to update a .csv with each row a different 
        # byc_expt_name (from byc_expt_dirs_dict)

        # This will broadcast the expt_name column (which starts as just one
        # row) to the length of the longest following column added
        # relpaths = utilities.get_rel
        relpaths = [utilities.get_relpath(path) for path in cell_trace_df_paths]
        trace_df = pd.DataFrame({'expt_name': expt_name,
                                 'trace_path': cell_trace_df_paths,
                                 'trace_relpath': relpaths,
                                 'chase_index': chase_index,
                                 'data': date})
        self.trace_database_df = pd.concat([self.trace_database_df, trace_df], ignore_index=True, sort=False)
        self.trace_database_df.to_csv(self.trace_database_path, index=False)

    def get_cell_trace_dfs(self, expt_name):
        # read the list of cell trace .csv paths in 
        # into a list of dataframes and return it.
        try:
            expt_df = self.trace_database_df[self.trace_database_df.expt_name == expt_name]
            traces_list = [pd.read_csv(os.path.join(constants.byc_data_dir, path)) for path in expt_df.trace_relpath]
        except Exception as e:
            print(f"Couldn't find trace .csvs for expt {expt_name} in constants.byc_data_dir\nChecking in absolute path")
            print(e)
            try:
                traces_list = [pd.read_csv(os.path.abspath(path)) for path in expt_df.trace_path]
            except:
                print(f'Could not find .csvs for expt {expt_name} in absolute paths')
                traces_list = None

        return traces_list

    def update_expt_trace_dfs(self, new_dfs_list, expt_name, **kwargs):
        """
        For each cell in new_dfs_list, overwrite its old trace .csv 
        with the one passed this function. Overwrite old trace csvs
        by default. If overwrite=False, put v2 on end of new filenames
        but don't add these new paths to the trace database.
        """

        overwrite = kwargs.get('overwrite', True)
        expt_df = self.trace_database_df[self.trace_database_df.expt_name == expt_name]
        # Since the expt_df slice from the main database df doesn't necessarily start
        # at 0. 
        indices = expt_df.index

        if overwrite:
            writepaths = [expt_df.trace_path[index] for index in indices]
        else:
            writepaths = [f'{expt_df.trace_path[index][:-4]}_v2.csv' for index in indices]

        message = f"Number of dfs ({len(new_dfs_list)})did not match number of writepaths in database ({len(writepaths)})"
        assert len(writepaths) == len(new_dfs_list), message 

        for i in range(len(new_dfs_list)):
            df = new_dfs_list[i]
            writepath = writepaths[i]
            df.to_csv(writepath, index=False)

def write_final_fits_dfs(drops=[], fits_df=None, mdf=None, **kwargs):
    """
    After preliminary analysis and data fitting, write a final
    'fits_table.csv' and 'allfitsdf.csv' file for the BYC compartment

    This is after all the manual work of inspecting fits to throw out
    bad cells/traces, segmentation (manual or auto), etc.

    Return allfitsdf, fits_table
    """
    fits_table_filename = kwargs.get('fits_table_filename', 'fits_table.csv')
    fits_and_traces_table_filename = kwargs.get('fits_and_traces_table_filename', 'allfitsdf.csv')
    return_dfs = kwargs.get('return_dfs', True)
    compindex = kwargs.get('compindex', 0)
    # Get rid of cells with 'cell_index' in the <drops> list
    fits_df = fits_df[~(fits_df.cell_index.isin(drops))]
    # `fits_df` contains fluorescence measurements at each timepoint
    # so make a pivot table of it called `fits_table`
    agg_idx = ['age_at_chase',
            'rls',
            'dist_from_sen',
            'first_bud_frame',
            'cell_index']
    for col in agg_idx:
        if col not in list(mdf.columns):
            print(f'Did not find {col} in master index df. Setting to NaN')
            mdf.loc[:, col] = np.nan
    # Label each cell in fits_df with information found in the
    # master index df (<mdf>)
    fits_df = sa.merge_dfs(mdf, fits_df)
    # Some datasets won't have all of the desired columns above, so if the column
    # is not in the fits_df, set all values of that column to np.nan
    final_agg_idx = []
    for col in agg_idx:
        # If the value is an np.nan, don't use it in the pivot table
        # index used to make the fits_table
        if False not in np.unique(np.isnan(fits_df[col])):
            print(f'Only found np.nan in {col} column of fits_df')
        else:
            final_agg_idx.append(col)
    print(f'Using {final_agg_idx} to make pivot table of fits_df')
    fits_table = pd.pivot_table(index=final_agg_idx, data=fits_df).reset_index()
    # Label `fits_table` because non-number data gets eliminated when
    # creating the pivot table
    fits_table = sa.merge_dfs(mdf, fits_table, idx=final_agg_idx)
    # Save the `fits_table` and `fits_df` files in the first compartmentdir
    # found in the master index df (<mdf>)
    date = mdf.date.unique()[0]
    compartment_name = mdf.compartment_name.unique()[0]
    compartment_dir = files.get_byc_compartmentdir(f'{date}_byc', compartment_name)
    compdir = os.path.join(constants.byc_data_dir, compartment_dir)
    fits_table_savepath = os.path.join(compdir, fits_table_filename)
    fits_df_savepath = os.path.join(compdir, fits_and_traces_table_filename)
    fits_table.to_csv(fits_table_savepath, index=False)
    print(f'Saved exponential fits table at {fits_table_savepath}')
    fits_df.to_csv(fits_df_savepath, index=False)
    print(f'saved expanded exponential fits df at {fits_df_savepath}')
    
    if return_dfs:
        return fits_df, fits_table
    else:
        pass

# I REALLY NEED TO CHANGE THIS FUNCTION TO A CLASS SO IT'S EASIER TO DEBUG
def get_byc_fits_df(**kwargs):
    """
    Return exponential fits dataframe read in using 'fits_table.csv' files 
    found in in all compartment directories of all byc experiments
    in constants.byc_data_dir
    
    If kwarg 'return_traces_df' is True, return 

    (fits_df, allfitsdf)

    allfitsdf contains fluorescence vs. time measurements, in addition to
    just the rate constants from fit.   
    """
    compartment_name_var = 'compartment_name'
    return_traces_df = kwargs.get('return_traces_df', False)

    fits_table_paths = files.get_fits_table_paths()
    fits_tables = [pd.read_csv(p) for p in fits_table_paths]
    for i, table in enumerate(fits_tables):
        # Add a the compartment name from the fits_table_path
        compartment_dir_Path = Path(fits_table_paths[i]).parents[0]
        expt_dir_Path = compartment_dir_Path.parents[0]
        date = re.search(patterns.date, expt_dir_Path.name).group()
        compartment_name = compartment_dir_Path.name
        table.loc[:, compartment_name_var] = compartment_name
        table.loc[:, 'date'] = date
        print(f'Found date {date}, compartment name {compartment_name}')

    allfitsdf_paths = files.get_allfitsdf_paths()
    allfitsdfs = [pd.read_csv(p) for p in allfitsdf_paths]
    for i, df in enumerate(allfitsdfs):
        # Add a the compartment name from the fits_table_path
        compartment_name = Path(allfitsdf_paths[i]).parents[0].name
        df.loc[:, compartment_name_var] = compartment_name


    fits_df = pd.concat(fits_tables, ignore_index=True)
    allfitsdf = pd.concat(allfitsdfs, ignore_index=True)
    # Drop rows without a compartment_name
    fits_df = fits_df.dropna(axis=0, subset=[compartment_name_var])
    if return_traces_df:
        return (fits_df, allfitsdf)
    else:
        return fits_df

def label_from_strain_database(features_df):
    """
    Look up the strain_name found in <features_df> in the table
    at constants.strains_db_path, annotate plasmid, substrate, 
    genetic background strain etc. in the features df
    
    Return annotated features_df
    """
    strains_df = pd.read_excel(constants.strains_db_path)
    strain_name = features_df.strain_name.iloc[0]
    if type(strain_name) == str:
        row = strains_df.set_index('Name').loc[strain_name, :]
        plasmid = row.Plasmid
        substrate = plasmids.substrate_name_from_plasmid_name(plasmid)
        background = row.Background

        features_df.loc[:, 'plasmid'] = plasmid
        features_df.loc[:, 'substrate'] = substrate
        features_df.loc[:, 'background'] = background
    else:
        print(f'Did not find suitable strain name: {strain_name}')
        
    return features_df

def label_from_manual_byc_index(df, **kwargs):
    """
    Pass a byc database dataframe (like those generated with 
    database.get_byc_fits_df()) or a buddf. For each row in the dataframe, 
    use the found compartment name to label the obesrvation using 
    the compartment name existing at constants.compartment_index_path
    """
    comp_index_path = (kwargs.get('comp_index_path', constants.compartment_index_path))
    if os.path.exists(comp_index_path):
        compdf = pd.read_csv(comp_index_path)
        print(f'Found file at {comp_index_path}')
    else:
        print(f'No file found at {comp_index_path}')
        return None
    
    compdf.set_index('compartment_name', inplace=True)
    for column in list(compdf.columns):
        df.loc[:, column] = np.nan
        
    for name in list(compdf.index):
        print(f'Labeling data from compartment {name}')
        for column in list(compdf.columns):
            print(f'Adding columns {column}')
            df.loc[df.compartment_name==name, column] = compdf.loc[name, column]

    return df

    
def features_df_from_compartment_name(compartment_name, **kwargs):
    """
    Using feature patterns (ie strain name, plasmid, substrate,
    genotype, etc) found in constants.patterns_dict, find
    features in <compartment_name> and create a dataframe with
    those features
    """
    return_dict = kwargs.get('return_dict', False)
    patterns_dict = constants.patterns_dict
    match_dict = {}
    for item in patterns_dict.items():
        col_name = item[0]
        pattern = item[1]
        match = re.search(pattern, compartment_name)
        if match != None:
            match_dict[col_name] = match.group()
        else:
            match_dict[col_name] = np.nan
    features_df = pd.DataFrame(match_dict, index=[0])
    # If there's a strain name in the compartment_name,
    # use it to look up the strain in the strains database
    # and label
    features_df = label_from_strain_database(features_df)

    if return_dict:
        return features_df, match_dict
    else:
        return features_df

def label_fits_df(fits_df):
    """
    For each unique compartment_name in the fits_df, label each row using 
    features found in that compartment_name
    
    Return the annotated fits_df
    """
    for column in constants.patterns_dict.keys():
        fits_df.loc[:, column] = np.nan

    features_dfs = []
    for compartment_name in fits_df.compartment_name.unique():
        print(compartment_name)
        features_df = features_df_from_compartment_name(compartment_name)
        features_dfs.append(features_df)
        for column in features_df.columns:
            value = features_df[column].iloc[0]
            fits_df.loc[fits_df.compartment_name==compartment_name, column] = value
            
    return fits_df

def set_cell_serial(df):
    """
    Label each row of the <df> dataframe using columns 'comparment_name',
    'cell_index', and 'dist_from_sen' joined by dashes. This allows
    distinguishing of cell traces that are from that same cell at different
    timepoints and both cells have the same 'cell_index'.

    I don't do this in datasets created after ~20211001. Instead each new
    trace gets a new 'cell_index' and same cells can be identified using a
    'bud_serial' which is all bud appearance frames joined by dashes

    Return the 'cell_serial" labelled <df>
    """
    df.loc[:, 'cell_serial'] = np.nan
    for idx in list(df.index.unique()):
        comp_name = df.loc[idx, 'compartment_name']
        cell_index = df.loc[idx, 'cell_index']
        dist_from_sen = df.loc[idx, 'dist_from_sen']
        cell_serial = '-'.join([comp_name, str(cell_index), str(dist_from_sen)])
        print(cell_serial)
        df.loc[idx, 'cell_serial'] = cell_serial
        
    return df

def set_bud_id(df, include_compartment=True):
    """
    
    """
    df.loc[:, 'bud_id'] = np.nan
    df.loc[:, 'compartment-bud_id'] = np.nan
    for idx in df.index.unique():
        if include_compartment:
            comp_name = df.loc[idx ,'compartment_name']
        if 'bud_rois_path' in df.columns:
            path = df.loc[idx, 'bud_rois_path']
        elif 'bud_roi_set_path' in df.columns:
            path = df.loc[idx, 'bud_roi_set_path']
        else:
            print(f'No column with bud roi set path')
            return None

        if not type(path) == float:
            # print(f'Reading in bud rois from\n{path}')
            if os.path.exists(path):
                bud_roi_df = files.read_rectangular_rois_as_df(path)
                bud_roi_serial = '-'.join([str(val) for val in bud_roi_df.position.values - 1])

                if include_compartment:
                    df.loc[idx, 'compartment-bud_id'] = f'{comp_name}-{bud_roi_serial}'
                df.loc[idx, 'bud_id'] = bud_roi_serial
            else:
                df.loc[idx, 'compartment-bud_id'] = np.nan
                df.loc[idx, 'bud_id'] = np.nan
        else:
                df.loc[idx, 'compartment-bud_id'] = np.nan
                df.loc[idx, 'bud_id'] = np.nan
    return df

class BudDataBase():
    """
    Database of all bud annotations in byc.constants.byc_data_dir
    which defaults to "<byc_source_dir>/data"

    Class has methods and attributes referring to where different files are found etc.

    Instantiate using a 'fits_df' generated using database.get_byc_fits_df()
    """
    _byc_data_dir = constants.byc_data_dir
    _byc_trace_database_path = constants.database_paths['cell_trace_db_path']

    def __init__(self, fits_df, **kwargs):
        
        
        self.time_delta_mins = 10
        self.fits_df = fits_df
        self.bud_roi_dfs = self.get_bud_roi_dfs()
        self.buddf = self.get_buddf()

    def get_bud_roi_dfs(self, use_rel_paths=False):
        if use_rel_paths:
            bud_roi_paths = [os.path.join(constants.byc_data_dir, relpath) for relpath in self.fits_df.bud_roi_set_relpath.unique()]
            xys = [xy for xy in self.fits_df.xy]
        else:
            bud_roi_paths = self.fits_df.bud_rois_path
            xys = [xy for xy in self.fits_df.xy]
        self.fits_df.loc[:, 'bud_rois_path'] = bud_roi_paths
        # Need to account for cells where there are no bud_rois annotated
        # (ie the value will be np.nan). I still want to keep the rates for 
        # downstream in the analysis because all np.nan cells are just
        # time 0 BYC cells that weren't followed
        bud_roi_dfs = []
        print(f'Len bud_roi_paths = {len(bud_roi_paths)}')
        print(f'Len xys = {len(xys)}')
        for i, path in enumerate(bud_roi_paths):
            if not type(path) == float:
                if os.path.exists(path):
                    bud_roi_df = files.read_rectangular_rois_as_df(path)
                    bud_roi_df.loc[:, 'bud_roi_set_path'] = path
                    bud_roi_df.loc[:, 'xy'] = xys[i]
                else:
                    print(f'Found non-float path that does not exist: {path}')
            else:
                bud_roi_df = pd.DataFrame(None)
            bud_roi_dfs.append(bud_roi_df)

        return bud_roi_dfs

    def get_buddf(self):
        """
        """
        bud_roi_cell_inds  = list(self.fits_df.cell_index)
        for i, df in enumerate(self.bud_roi_dfs):
            cell_index = bud_roi_cell_inds[i]
            if not df.empty:
                print(f'Annotating cell cycle features for cell {cell_index}')
                # Make sure that bud appearance slices ('position') are in
                # correct order so that deltas make sense
                df.sort_values(by='position', ascending=True, inplace=True)
                df.loc[:, 'cycle_duration_frames'] = np.nan
                df.loc[:, 'cycle_duration_hrs'] = np.nan
                df.loc[:, 'bud_serial'] = np.nan
                df.loc[:, 'bud_frame'] = df['position'] - 1
                df.loc[:, 'first_bud_frame'] = df.bud_frame.min()
                df.loc[:, 'bud_time_hours'] = (df.bud_frame*self.time_delta_mins)/60
                bud_frames_str_list = [str(int(val)) for val in df.bud_frame.values]
                df.loc[:, 'bud_serial'] = '-'.join(bud_frames_str_list)
                # Last "bud" is actually the frame before the cell dies, and
                # first frame is the "beginning" of the first cell cycle (the
                # frame at which bud appears) so cell cycle 0 needs to have the 
                # first duration
                # Last bud observed to start doesn't truly have an end time observed
                # because the last 'bud' frame is the last frame when the cell was
                # seen alive
                df.loc[df.index[0:-2], 'cycle_duration_frames'] = np.diff(df['position'])[0:-1]
                cycle_duration_hours = (df.cycle_duration_frames*self.time_delta_mins)/60
                df.loc[:, 'cycle_duration_hrs'] = cycle_duration_hours
                # Just in case the index is not zero to len(df)
                df.index = range(len(df))
                # Add median and mean filtered cycle duration columns
                constant_args = [
                    df,
                    'cycle_duration_hrs'
                ]
                name_with_kernel = True
                for kernelsize in [3, 5]:
                    trace_tools.mean_filter(*constant_args, kernelsize, name_with_kernel)
                    trace_tools.median_filter(*constant_args, kernelsize, name_with_kernel)
                # Create a column with the shape of each bud annotated
                pattern = "(round|long)"
                matches = df['name'].apply(lambda x: re.search(pattern, x))
                df.loc[0:, 'shape_matches'] = matches
                # This function will make we add nan for bud roi frames that
                # for whatever reason did not have a shape annotated
                def find_group(x):
                    if x is None:
                        return np.nan
                    else:
                        return x.group()
                df.loc[:, 'bud_shape'] = np.nan
                # Daughter shapes are annotated at bud appearance frame of the following
                # bud. I.e. the shape of bud i that appears at frame x is described
                # in the roi at frame y where bud i + 1 appears. The first value
                # in the shapes array is meaningless because it would describe a bud
                # that we never saw appear. The last ROI in the bud roi set annotates
                # the last frame before the cell dies. So that shape annotation stays
                df.loc[df.index[0:-1], 'bud_shape'] = df.shape_matches.apply(lambda x: find_group(x))[1:]
                n_round_buds = len(df[df.bud_shape=='round'])
                n_long_buds = len(df[df.bud_shape=='long'])
                df.loc[:, 'n_long_buds'] = n_long_buds
                df.loc[:, 'n_round_buds'] = n_round_buds
                df.loc[:, 'cycle_number'] = range(len(df))
                # Last index in the bud roi set is the pre-death frame so we'll give that
                # a cycle number of nan
                df.loc[df.index[-1], 'cycle_number'] = np.nan
                df.loc[:, 'rls_fraction'] = df['cycle_number']/df['cycle_number'].max()
                df.loc[:, 'dist_from_sen'] = df['cycle_number'].max() - df['cycle_number']
                # Last annotated "bud frame" is the frame at which the cells was last observed
                # alive, whether because it dies or disappears from view down the central
                #  hallway in the next (or almost next) frame.
                df.loc[:, 'rls'] = len(df) - 1
                # Annotate how the buds observation ended (either cell escaped or died during the
                # course of the experiment)
                # need to do this by reading in the annotated bud_rois_df.csv corresponding to 
                # this cell.
                bud_roi_path = df.bud_roi_set_path.iloc[0]
                bud_roi_fn = os.path.basename(bud_roi_path)
                xy = int(df.xy.iloc[0])
                rois_df_fn = f'{bud_roi_fn[0:12]}_xy{str(xy).zfill(2)}cell{str(int(cell_index)).zfill(3)}_bud_rois_df.csv'
                rois_df_path = bud_roi_path.replace(bud_roi_fn, rois_df_fn)
                bud_rois_annotated_df = pd.read_csv(rois_df_path)
                for col in bud_rois_annotated_df.columns:
                    df.loc[:, col] = bud_rois_annotated_df[col].iloc[0]
            else:
                df = pd.DataFrame({'cell_index': cell_index}, index=[0])
            df.loc[:, 'compartment_name'] = os.path.basename(os.path.dirname(rois_df_path))

        buddf = pd.concat(self.bud_roi_dfs, ignore_index=True)
        
        return buddf

def annotate_daughter_shapes(bud_rois_df):
    """
    Find the text "round" or "long" in the name column of
    the <bud_rois_df> and annotate the data found in place

    The <bud_rois_df> final frame or position row should
    be the death frame
    """    
    # Create a column with the shape of each bud annotated
    pattern = "(round|long)"
    matches = bud_rois_df['name'].apply(lambda x: re.search(pattern, x))
    bud_rois_df.loc[0:, 'shape_matches'] = matches
    # This function will add nan for bud roi frames that
    # for whatever reason did not have a shape annotated
    def find_group(x):
        if x is None:
            return np.nan
        else:
            return x.group()
    bud_rois_df.loc[:, 'bud_shape'] = np.nan
    # Daughter shapes are annotated at bud appearance frame of the following
    # bud. I.e. the shape of bud i that appears at frame x is described
    # in the roi at frame y where bud i + 1 appears. The first value
    # in the shapes array is meaningless because it would describe a bud
    # that we never saw appear. The last ROI in the bud roi set annotates
    # the last frame before the cell dies. So that shape annotation stays
    bud_rois_df.loc[:, 'bud_shape'] = bud_rois_df.shape_matches.apply(lambda x: find_group(x))
    # Offset the bud shapes because the shape annotated at the frame at which
    # bud i appears is the shape of bud i-1
    bud_rois_df.loc[0:bud_rois_df.index.max()-1, 'bud_shape'] = bud_rois_df.bud_shape[1:].values
    bud_rois_df.loc[bud_rois_df.index.max(), 'bud_shape'] = np.nan

def get_crop_roi_df_from_cell_index_compdir(cell_index, compartment_dir):

    crop_roi_df_fns = [fn for fn in os.listdir(compartment_dir) if 'crop_rois_df.csv' in fn]
    crop_roi_df_fn = [fn for fn in crop_roi_df_fns if str(cell_index).zfill(3) in fn][0]
    crop_roi_df_path = os.path.join(compartment_dir, crop_roi_df_fn)
    crop_roi_df = pd.read_csv(crop_roi_df_path)

    return crop_roi_df

if __name__ == '__main__':
    byc_database = DataBase()