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
            print(f'Found unequal number of expt dates and master indices')
            print(dates)
            for path in self.all_master_idx_paths:
                print(path)
            paths_dict = None
            
        return paths_dict
    
    def master_idx_by_date(self, exptdate, timelapse=False):
        """
        Create a master_index using information found
        in the exptdir
        """
        path = self.paths_dict[exptdate]
        datadir = os.path.join(os.path.dirname(path), 'data')
        os.path.exists(datadir)
        if not timelapse:
            fns = os.listdir(datadir)
        else:
            dirs = os.listdir(datadir)
        # Create a master idx dataframe based on the files found
        # in this experiments datadir
        strains = []
        filepaths = []

        for fn in fns:
            print(fn)
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

    def exptdf(self, exptdate, **kwargs):
        """
        Return a dataframe holding all flow observations
        found according to the master_idx

        Optionally pass a master_index_df kwarg to avoid
        trying to automatically set a master_index_df
        """
        if 'master_index_df' in kwargs:
            master_idx = kwargs['master_index_df']
        else:
            master_idx = self.master_idx_by_date(exptdate)    

        sampledfs = []
        # Read in data and add identifying information
        # based on master index
        print(f'Found master index with {len(master_idx)}')
        for idx in master_idx.index:
            row = master_idx.loc[idx, :]
            print(f'Looking for data at {row.filepath}')

            if os.path.exists(row.filepath):
                print(f'Found data')
                sampledf = FCMeasurement(ID=f'{row.strain}-{row.clone}', datafile=row.filepath).data
                print(f'Found {len(sampledf)} measurements in this file')
                # Annotate sample df
                for col in row.index:
                    sampledf.loc[:, col] = row.loc[col]
                sampledfs.append(sampledf)
            else:
                print(f'No data found')

        if len(sampledfs) > 0:
            exptdf = pd.concat(sampledfs, ignore_index=True)
        else:
            exptdf = None
            print(f'No data found for exptdate {exptdate}')

        return exptdf

def label_from_mdf(mdfpath, alldf):
    mdf = pd.read_csv(mdfpath)
    # Add data labels to alldf based on the master index    
    mdf.set_index('sample_id', inplace=True)
    for col in mdf.columns:
        alldf.loc[:, col] = np.nan

    for sample_id in alldf.sample_id.unique():
        print(sample_id)
        if sample_id in list(mdf.index):
            print('Found entry in master index')
            # Add information from master index to
            # raw flow data
            for col in mdf.columns:
                alldf.loc[alldf.sample_id==sample_id, col] = mdf.loc[sample_id, col]
        else:
            print(f'sample_id {sample_id} not found in master index')
            print(f'sample_ids in master index:\n{mdf.index}')
    mdf.reset_index(inplace=True)
    return alldf

def label_from_substrate_index(alldf, substrate):
    
    substrates_index = pd.read_csv(constants.substrates_index_path)
    substrates_index.set_index('Substrate', inplace=True)

    for expr_method in substrates_index.columns:
        plasmid_number = substrates_index.loc[substrate, expr_method]
        if pd.notna(plasmid_number):
            plasmid = f'pJC{str(int(plasmid_number)).zfill(3)}'
            print(f'Found plasmid {plasmid}')
            if plasmid in alldf.plasmid.unique():
                print(f'Found {plasmid} in data')
                alldf.loc[alldf.plasmid==plasmid, 'substrate'] = substrate
                alldf.loc[alldf.plasmid==plasmid, 'expr_method'] = expr_method.replace('_URA_int', '')
            else:
                print(f'No data found for {plasmid}')
    substrates_index.reset_index(inplace=True)
    return alldf

def label_alldf(mdf, alldf, label_from_mdf_=True):

    substrates_index = pd.read_csv(constants.substrates_index_path)
    if label_from_mdf_:
        alldf = label_from_mdf(mdf, alldf)
    alldf.loc[:, 'substrate'] = np.nan
    alldf.loc[:, 'expr_method'] = np.nan
    for substrate in substrates_index.Substrate:
        alldf = label_from_substrate_index(alldf, substrate)
    
    return alldf

def background_subtract(alldf, **kwargs):
    """
    For each sample in the alldf, find its corresponding no-plasmid
    background sample and subtract that background sample's median
    fluorescence for each channel name in <channels>
    
    Return the alldf with background subtracted
    
    Chase method is necessary because of high background caused by 
    small molecules like tetracyline in DMSO or higher autofluorescence
    caused by integrated fluorescence across a larger cell size
    """
    classifier_list = ['background',
                      'Chase_method']
    classifers = kwargs.get('classifiers', classifier_list)
    channels_list = ['Alexa Fluor 488-A',
                'DsRed-A']
    channels = kwargs.get('channels', channels_list)
    # Filter out data that isn't properly labeled (i.e. from
    # plate wells that were collected but didn't actually contain
    # sample)
    for c in classifers:
        alldf = alldf[alldf[c].isna()==False]
    # All unique combinations of classifiers to use to
    # match samples to their background control
    no_plasmid_df = alldf.set_index(['plasmid']).loc['no-plasmid', :].reset_index()
    unique_sets = alldf.set_index(classifers).index.unique()
    print(f'Found unique sets of classfiers: {unique_sets}')
    # Set index to classifers columsn to make it quicker to
    # select data matching each background sample
    alldf.set_index(classifers, inplace=True)
    for channel in channels:
        newcolname = f'{channel}_bg_sub'
        normcolname = f'{channel}_bg_norm'
        alldf.loc[:, newcolname] = np.nan
        for unique_set in unique_sets:        
            bg_fluor_vals = no_plasmid_df.set_index(classifers).loc[unique_set, channel]
            bg_fluor_median = bg_fluor_vals.median()        
            bg_sub_vals = alldf.loc[unique_set, channel] - bg_fluor_median
            bg_norm_vals = alldf.loc[unique_set, channel]/bg_fluor_median
            alldf.loc[unique_set, newcolname] = bg_sub_vals
            alldf.loc[unique_set, normcolname] = bg_norm_vals
            print(f'Background for {unique_set} {channel} = {bg_fluor_median}')

    alldf.reset_index(inplace=True)
    
    return alldf