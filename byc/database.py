import os
import pandas as pd

class dataBase():
    _byc_data_dir = r'C:\Users\John Cooper\Projects\byc_data'
    # Set steady state data dir
    _ss_data_dir = r"C:\Users\John Cooper\Box Sync\Finkelstein-Matouschek\images"
    # Set flow cytometry data dir
    _fc_dat_dir = r"C:\Users\John Cooper\Box Sync\Finkelstein-Matouschek\images"
    _byc_trace_database_dir = r'C:\Users\John Cooper\Projects\byc\docs'

    def __init__(self):

        expt_names = os.listdir(self.byc_data_dir)
        expt_paths = [f'{self.byc_data_dir}\\{folder}' for folder in expt_names]
        byc_expt_dirs_dict = dict(zip(expt_names, expt_paths))  
        self._byc_expt_dirs_dict = byc_expt_dirs_dict
        self.trace_database_path = f'{self._byc_trace_database_dir}\\cell_trace_database.csv'
        self.trace_database_df = pd.read_csv(self.trace_database_path)

    @property
    def byc_data_dir(self):
        return self._byc_data_dir

    @byc_data_dir.setter
    def byc_data_dir(self, byc_data_dir):
        self._byc_data_dir = byc_data_dir

    @property
    def byc_expt_dirs_dict(self):
        return self._byc_expt_dirs_dict

    @property
    def trace_database_df(self):
        return self._trace_database_df

    @trace_database_df.setter
    def trace_database_df(self, trace_database_df):
        # Can set constraints here for making sure 
        # self.trace_database_df is set appropriately
        self._trace_database_df = trace_database_df

    def initialize_cell_trace_database(self, overwrite=True):
        """
        Completely clear the database .csv file and save a blank
        one with columns passed here. Reset self.trace_database_df 
        to the initialized df. 
        """
        # 
        cols = ['expt_name', 'trace_path', 'chase_index']
        trace_db_df = pd.DataFrame(columns = cols)
        if overwrite:
            writepath = f'{self._byc_trace_database_dir}\\cell_trace_database.csv'
        else:
            writepath = f'{self._byc_trace_database_dir}\\cell_trace_database_new.csv'

        trace_db_df.to_csv(write_path, index=False)
        self.trace_database_df = pd.read_csv(writepath)
        
    def add_expt_to_trace_database(self, cell_trace_df_paths, expt_name, chase_index):
        # Use the list of file names chosen by the user
        # to update a .csv with each row a different 
        # byc_expt_name (from byc_expt_dirs_dict)

        # This will broadcast the expt_name column (which starts as just one
        # row) to the length of the longest following column added
        trace_df = pd.DataFrame({'expt_name': expt_name,
                                 'trace_path': cell_trace_df_paths,
                                 'chase_index': chase_index})
        self.trace_database_df = pd.concat([self.trace_database_df, trace_df], ignore_index=True, sort=False)
        self.trace_database_df.to_csv(self.trace_database_path, index=False)

    def get_cell_trace_dfs(self, expt_name):
        # read the list of cell trace .csv paths in 
        # expt_dirs_dict into a list of dataframes and
        # return it.
        try:
            expt_df = self.trace_database_df[self.trace_database_df.expt_name == expt_name]
            traces_list = [pd.read_csv(trace_path) for trace_path in expt_df.trace_path]
        except:
            print(f'Could not find .csvs for expt {expt_name}')
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

byc_database = dataBase()