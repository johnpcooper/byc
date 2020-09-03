import re
import os
import pandas as pd
from snapgene_reader import snapgene_file_to_dict
from byc import files, utilities, constants

class Plasmids(object):    
    """
    A database of plasmids found in constants.plasmids_dir
    with the prefix 'pJC' + three digits
    """
    def __init__(self, choose_dir=False):
        # Set the minimum length required for
        # coding sequences found in the snapgene
        # file. Filters out leftover features
        # and linkers which is useful for extracting
        # constructs from plasmid name
        self.min_cds_length = 30        
        if choose_dir:
            # Functionality to manually select dir
            # here
            pass            
        else:
            self.plasmids_dir = constants.plasmids_dir
            self.plasmid_filenames = self.get_plasmid_filenames()
            self.plasmid_paths = self.get_plasmid_paths()
            self.snapgene_dicts = self.get_snapgene_dicts()
            self.plasmid_names = self.get_plasmid_names()
            self.plasmids_dict = self.get_plasmids_dict()
            self.features_dict = self.get_features_dict()
            
    def get_plasmid_filenames(self):
        """
        Return the list of snapgene .dna filenames
        in self.plasmids_ir
        """
        
        if os.path.exists(self.plasmids_dir):
            filetype = '.dna'
            filenames = os.listdir(self.plasmids_dir)
            filenames = [file for file in filenames if filetype in file]
            
        else:
            print(f"WARNING: path does not exist {self.plasmids_dir}")
            filenames = []
            
        return filenames
            
    def get_plasmid_paths(self):
        """
        Return the list of paths to snapgene files in
        constants.plasmids_dir
        """
        paths = [os.path.join(self.plasmids_dir, fn) for fn in self.plasmid_filenames]
        paths = [path for path in paths if os.path.exists(path)]
        return paths

    def get_snapgene_dicts(self):
        """
        Return the list of snapgene file dictionaries
        read from self.plasmid_paths
        """
        snapgene_dicts = [snapgene_file_to_dict(path) for path in self.plasmid_paths]
        return snapgene_dicts

    def get_plasmid_names(self):
        """
        Return the list of pJC000 etc. that
        currently exist in constants.plasmids_dir
        """
        plasmid_names = []
        for name in self.plasmid_filenames:
            match = re.search(constants.patterns.plasmid_name, name)
            if match:
                plasmid_names.append(match.group())

        return plasmid_names

    def get_plasmids_dict(self):
        """
        Return a dictionary where keys are plasmid names
        ('pJC001' etc.) and values are snapgene file dicts
        """
        d = dict(zip(self.plasmid_names, self.snapgene_dicts))
        return d
    
    def get_features_dict(self):
        """
        Return a dict where plasmid names ('pJC001' etc.) are keys
        referring to dataframes made from each plasmid's snapgene
        dict feature set
        """
        features_dict = {}
        for plasmid_name, sg_dict in self.plasmids_dict.items():
            features_df = pd.DataFrame(sg_dict['features'])
            features_dict[plasmid_name] = features_df.sort_values(by='end').reset_index(drop=True)
            
        return features_dict