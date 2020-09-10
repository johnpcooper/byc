import re
import os

import numpy as np
import pandas as pd

from FlowCytometryTools import FCMeasurement

from byc import constants, utilities, files

class Experiments(object):
    """
    Look up and instantiate information about a flow
    experiment by date
    """
    def __init__(self):
        
        self.all_master_idx_paths = self.get_all_master_idx_paths()
        self.all_master_idx_dfs = self.get_all_master_idx_dfs()
        self.paths_dict = self.get_paths_dict()
        
    def get_all_master_idx_paths(self):
        """
        Return a list of directories containing
        master indices in constants.flow_data_dir
        """
        paths = utilities.get_all_master_index_paths(rootdir=constants.flow_data_dir)
        return paths
        
    def get_all_master_idx_dfs(self):
        """
        Return a list of all the master index
        dfs found in constants.flow_data_dir
        """
        dfs = utilities.get_all_master_index_dfs(master_index_paths=self.all_master_idx_paths)
        return dfs
    
    def get_paths_dict(self):
        """
        Make a dictionary where you can look up flow
        master_index path by experiment date
        """
        dates = []
        datepattern = constants.patterns.date
        
        for path in self.all_master_idx_paths:
            match = re.search(datepattern, path)
            if match:
                dates.append(match.group())
                
        if len(dates) == len(self.all_master_idx_paths):
            paths_dict = dict(zip(dates, self.all_master_idx_paths))
        else:
            paths_dict = None
            
        return paths_dict
    
    def master_idx_by_date(self, exptdate):
        """
        Create a master_index using information found
        in the exptdir
        """
        path = self.paths_dict[exptdate]
        datadir = os.path.join(os.path.dirname(path), 'data')
        os.path.exists(datadir)
        fns = os.listdir(datadir)
        # Create a master idx dataframe based on the files found
        # in this experiments datadir
        strains = []
        filepaths = []

        for fn in fns:
            if fn[-4:] == '.fcs':
                match = re.search(constants.patterns.strain_name, fn)
                if match:
                    strains.append(match.group())
                    filepath = os.path.join(datadir, fn)
                    filepaths.append(filepath)

        df = pd.DataFrame({'strain': strains,
                           'filepath': filepaths})
        # Add clone indices to the dataframe
        for strain in df.strain.unique():

            n_clones = len(df[df.strain == strain])
            df.loc[df.strain == strain, 'clone'] = [int(idx) for idx in range(1, n_clones+1, 1)]

        # Lookup each strain in constants.strains_dir/Strains_Database.csv
        # and add information found in the database
        strains_df = pd.read_csv(os.path.join(constants.strains_dir, 'Strains_Database.csv'))

        for idx in df.index:
            strain_name = df.loc[idx, 'strain']
            if strain_name in strains_df.name.values:
                for col in strains_df.columns:
                    df.loc[idx, col] = strains_df.loc[strains_df.name == strain_name, col].values[0]
                    
        return df

    def exptdf(self, exptdate):
        """
        Return a dataframe holding all flow observations
        found according to the master_idx
        """
        master_idx = self.master_idx_by_date(exptdate)
        sampledfs = []
        # Read in data and add identifying information
        for idx in master_idx.index:
            row = master_idx.loc[idx, :]    
            sampledf = FCMeasurement(ID=f'{row.strain}-{row.clone}', datafile=row.filepath).data
            
            for col in row.index:
                sampledf.loc[:, col] = row.loc[col]
                
            sampledfs.append(sampledf)
            
        exptdf = pd.concat(sampledfs, ignore_index=True)

        return exptdf