import os
import pandas as pd

from byc import constants, utilities

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
        expt_paths = [f'{self.byc_data_dir}\\{folder}' for folder in expt_names]
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
        except:
            print(f"Couldn't find trace .csvs for expt {expt_name} in constants.byc_data_dir\nChecking in absolute path")
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

byc_database = DataBase()


