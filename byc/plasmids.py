import re
import os
import pandas as pd
import numpy as np
from snapgene_reader import snapgene_file_to_dict
from byc import files, utilities, constants

class Plasmids(object):    
    """
    A database of plasmids found in constants.plasmids_dir
    """
    def __init__(self, choose_dir=False, read_files=False):
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
            self.plasmid_names = self.get_plasmid_names()
            self.paths_dict = dict(zip(self.plasmid_names, self.plasmid_paths))

            if read_files:
                self.snapgene_dicts = self.get_snapgene_dicts()                
                self.plasmids_dict = self.get_plasmids_dict()
                self.features_dict = self.get_features_dict()
            
    def get_plasmid_filenames(self):
        """
        Return the list of snapgene .dna filenames
        in self.plasmids_dir
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
        plasmids_subdirpaths = [os.path.join(self.plasmids_dir, subdir) for subdir in constants.plasmid_subdirs]

        snapgene_filepaths = []
        for dirpath in plasmids_subdirpaths:
            print(f'Checking {dirpath} for .dna files')
            if os.path.exists(dirpath):

                pattern = r'(.*)(.dna$)'
                filenames = os.listdir(dirpath)
                filenames = [fn for fn in filenames if re.search(pattern, fn) != None]
                paths = [os.path.join(dirpath, fn) for fn in filenames]
                plasmids_subdirpaths = [path for path in paths if os.path.exists(path)]
                
            else:
                print(f"WARNING: path does not exist {dirpath}")
                plasmids_subdirpaths = []

            snapgene_filepaths.extend(plasmids_subdirpaths)

        return snapgene_filepaths

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
        for name in self.plasmid_paths:
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


def ng_ul_to_fmol_ul(dsDNA_len_bp, dsDNA_concn_ng_ul):
    """
    Return DNA concentration in femtomoles per microliter. 
    Takes arrays or single values
    
    X = dsDNA_concn_ng_ul
    Y = dsDNA_len_bp
    
                          X ng | 1e-9 g |    1 mole   | 1e15 femtomole
    dsDNA_concn_fmol_ul = --------------------------------------------
                          1 ul |  1 ng  | 650g * Y bp |    1 mole
    """
    g_per_bp = 650
    g_per_base = g_per_bp/2

    dsDNA_concn_fmol_ul = (dsDNA_concn_ng_ul*np.power(10, 6))/ (g_per_bp*dsDNA_len_bp)
    return dsDNA_concn_fmol_ul


def extend_nanodrop_df(spec_df_path, target_fmol_per_ul=40, writeoutput=True):
    """
    
    """
    specdf =pd.read_csv(spec_df_path)
    plsmds = Plasmids()

    if 'Sample ID' in specdf.columns and 'Nucleic Acid' in specdf.columns:

        plasmid_names = specdf.loc[:, 'Sample ID']
        plasmid_paths = [plsmds.paths_dict[name] for name in plasmid_names]

        snapgene_dicts = [snapgene_file_to_dict(path) for path in plasmid_paths]
        plasmid_size_bp = [len(d['seq']) for d in snapgene_dicts]

        specdf.loc[:, 'plasmid_path'] = plasmid_paths
        specdf.loc[:, 'plasmid_size_bp'] = plasmid_size_bp

        fmol_ul = ng_ul_to_fmol_ul(specdf.plasmid_size_bp, specdf.loc[:, 'Nucleic Acid'])
        specdf.loc[:, 'fmol_per_ul'] = fmol_ul

        # Calcuate dilution factor needed to get to target_fmol_per_ul fmol/ul
        dil_factors = specdf.fmol_per_ul/target_fmol_per_ul
        specdf.loc[:, f'dil_factor_for_{target_fmol_per_ul}_fmol_per_ul'] = dil_factors

        if writeoutput:
            specdf.to_csv(spec_df_path)
        return specdf

    else:
        print("Sample ID and Nucleic Acid not found in dataframe columns")
        return None