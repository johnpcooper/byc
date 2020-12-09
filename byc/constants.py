import os
import pandas as pd

def package_path(**kwargs):
    """
    Return the path to the local installation
    of byc
    """
    import byc
    package_path = byc.__file__
    package_path = os.path.dirname(package_path)
    return package_path

def source_path(**kwargs):
    """
    Return the path to the local installation
    source of byc. (Includes 'data', 'notebooks',
    'bin', 'docs', etc.
    """
    import byc
    dirname = os.path.dirname(byc.__file__)
    source_path = os.path.dirname(dirname)
    return source_path
    
package_path = package_path()
source_path = source_path()

database_paths = {'cell_trace_db_path': os.path.join(package_path, 'cell_trace_database.csv'),
                  'cell_roi_db_path': os.path.join(package_path, 'cell_roi_database.json')}

byc_data_dir = os.path.join(source_path, 'data\\')
legacy_byc_data_dir = "C:\\Users\\John Cooper\\Projects\\byc_data\\"
steady_state_data_dir = "C:\\Users\\John Cooper\\Box Sync\\Finkelstein-Matouschek\\images\\"
flow_data_dir = "C:/Users/John Cooper/Box Sync/Finkelstein-Matouschek/flow_cytometry"

master_index_cols = ['date',
                     'expt_type',
                     'channels_collected',
                     'collection_interval'
                     'channel_exposures',
                     r'tetracycline[]',
                     'tetracycline_perfusion_frames',
                     'genetic_background',
                     'construct',
                     'expr_type',
                     'xy',
                     'cell_index',
                     'death',
                     'sen_start',
                     'sen_end',
                     'note',
                     'late_daughter_shape',
                     'Pos',
                     'dist_from_sen']

ss_master_index_cols = ['date',
                        'expt_type',
                        'fluor_channel_names'
                        'channel_names',
                        'raw_channel_names',
                        'channel_exposures',
                        r'tet_concn_uM',
                        'n_end',
                        'c_end',
                        'plasmid',
                        'clone',
                        'exptdir',
                        'genetic_background',
                        'construct',
                        'expr_type',
                        'substrate_ctl',
                        'note']


default_channels = ['yfp', 'dsred']

# Some constants for quickly accessing construct names and features
plasmids_dir = "C:/Users/John Cooper/Box Sync/Finkelstein-Matouschek/yeast_engineering/plasmids/"
plasmid_subdirs = ['YeastToolkit', 'JC', 'BLS']
ytk_dir = "C:/Users/John Cooper/Box Sync/Finkelstein-Matouschek/yeast_engineering/plasmids/YeastToolkit"            
strains_dir = "C:/Users/John Cooper/Box Sync/Finkelstein-Matouschek/yeast_engineering/strains"

class Patterns(object):

    def __init__(self):

        self.master_index_file = self.get_master_index_filename()
        self.date = self.get_date()
        self.plasmid_name = self.get_plasmid_name()
        self.strain_name = self.get_strain_name()
        self.crop_roi_df_file = self.get_crop_roi_df_file()
        self.measurement_roi = self.get_measurement_roi_zip()
        self.crop_roi = self.get_crop_roi_zip()
        self.bud_roi_df_file = self.get_bud_roi_df_file()

    def get_master_index_filename(self):
        """
        Return the pattern used in a regex search
        to identify master index files
        """
        # Has string 'master_index' anyhwere in string,
        # ends with '.csv'
        pattern = r"(.*)(master_index)(.*)(.csv$)"
        return pattern

    def get_date(self):
        """
        Return the pattern for dates. 
        Looks like: '20200723'
        """
        centuries = '|'.join([str(num) for num in range(19, 22)])
        decades_years = '|'.join([str(num).zfill(2) for num in range(0, 100)])
        months = '|'.join([str(num).zfill(2) for num in range(0, 12)])
        days = '|'.join([str(num).zfill(2) for num in range(0, 31)])

        pattern = f"({centuries})({decades_years})({months})({days})"
        return pattern

    def get_plasmid_name(self):
        """
        Return the plasmid name pattern:
        'pJC', 'pYTK', 'BLS', or 'pJBK' + 3 digits
        """
        indices = '|'.join([str(num).zfill(3) for num in range(0, 1000)])
        pattern = f'(pJC|pYTK|pBLS|pJBK)({indices})'
        return pattern

    def get_strain_name(self):
        """
        Return the strain name pattern
        `JPC` or `BLS` + 3 digits
        """
        indices = '|'.join([str(num).zfill(3) for num in range(0, 1000)])
        pattern = f'(BLS|JPC)({indices})'
        return pattern

    def get_crop_roi_df_file(self):
        """
        Return the crop roi df csv pattern
        These types of files are generated when
        Fiji runs python script byc.imagejpc.addcellroi.py
        in the plugin imagejpc/utilities/save_cell_roi_set.py
        and 'crop' roi type is selected
        """
        pattern = r'(.*)(crop_rois_df)(.*)(.csv$)'
        return pattern

    def get_crop_roi_zip(self):
        """
        Return the crop roi df csv pattern
        These types of files are generated when
        Fiji runs python script byc.imagejpc.addcellroi.py
        in the plugin imagejpc/utilities/save_cell_roi_set.py
        and 'crop' roi type is selected
        """
        pattern = r'(.*)(crop_rois)(.*)(.zip|.roi$)'
        return pattern

    def get_measurement_roi_df_file(self):
        """
        Return the measurement roi df csv pattern
        These types of files are generated when
        Fiji runs python script byc.imagejpc.addcellroi.py
        in the plugin imagejpc/utilities/save_cell_roi_set.py
        and 'measurement' roi type is selected
        """
        pattern = r'(.*)(measurement_rois_df)(.*)(.csv$)'
        return pattern

    def get_measurement_roi_zip(self):
        """
        Return the measurement roi .zip or .roi
        filename pattern
        These types of files are generated when
        Fiji runs python script byc.imagejpc.addcellroi.py
        in the plugin imagejpc/utilities/save_cell_roi_set.py

        These are stacks of single cell outlines over a byc
        experiment created in Fiji based on cropped cell stacks
        output from byc/bin/run_segmentation.py
        and 'measurement' roi type is selected
        """
        pattern = r'(.*)(measurement_rois)(.*)(.zip|.roi$)'
        return pattern

    def get_bud_roi_df_file(self):
        """
        Return the bud roi df csv pattern
        These types of files are generated when
        Fiji runs python script byc.imagejpc.addcellroi.py
        in the plugin imagejpc/utilities/save_cell_roi_set.py
        and 'bud' roi type is selected
        """
        pattern = r'(.*)(bud_rois_df)(.*)(.csv$)'
        return pattern

patterns = Patterns()