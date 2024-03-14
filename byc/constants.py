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
# Define data loacations for various types of experiments. Note the ones
# that need to be fully hardcoded - steady_state_data_dir, steady_state_rep_image_dir,
# and flow_data_dir
byc_data_dir = os.path.join(source_path, 'data/')
compartment_index_path = os.path.join(byc_data_dir, 'meta/compartments_index.csv')
steady_state_data_dir = "C:/Users/johnp/Box/Finkelstein-Matouschek/images/"
steady_state_rep_image_dir  = os.path.join(source_path, 'notebooks')
flow_data_dir = "C:/Users/johnp/Box/Finkelstein-Matouschek/flow_cytometry"
plots_dir = os.path.join(source_path, 'data/meta/plots')

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


default_fluor_channels = ['yfp', 'rfp']
default_channel_names = ['bf', 'yfp', 'rfp']
all_channel_names = ['bf', 'yfp', 'rfp', 'mko', 'bfp', 'gfp']
default_raw_channel_names = ['Brightfield', 'YFP', 'RFP']

# Some constants for quickly accessing construct names and features
plasmids_dir = "C:/Users/johnp\\Box/Finkelstein-Matouschek/yeast_engineering/plasmids/"
plasmid_subdirs = ['YeastToolkit', 'JC', 'BLS']
ytk_dir = "C:/Users/johnp\\Box/Finkelstein-Matouschek/yeast_engineering/plasmids/YeastToolkit"            
strains_dir = "C:/Users/johnp\\Box/Finkelstein-Matouschek/yeast_engineering/strains"
strains_db_path = r'C:\Users\johnp\\Box\Finkelstein-Matouschek\yeast_engineering\strains\JPC_Strains.xlsx'
substrates_index_path = r'C:\Users\johnp\\Box\Finkelstein-Matouschek\images\meta_analysis\Analysis\Substrates_Index.csv'

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
        self.tet_concn = self.get_tet_concn()
        self.estradiol_concn = self.get_estradiol_concn()
        self.genotype = self.get_genotype()
        self.culture_condition = self.get_culture_condition()
        self.timepoint_mins = self.get_timepoint_mins()
        self.clone_number = self.get_clone_number()
        self.cell_index = self.get_cell_index()
        self.channel_name = self.get_channel_name()

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
        months = '|'.join([str(num).zfill(2) for num in range(0, 13)])
        days = '|'.join([str(num).zfill(2) for num in range(0, 32)])

        pattern = f"({centuries})({decades_years})({months})({days})"
        return pattern

    def get_plasmid_name(self):
        """
        Return the plasmid name pattern:
        'pJC', 'pYTK', 'BLS', or 'pJBK' + 3 digits
        """
        indices = '|'.join([str(num).zfill(3) for num in range(0, 1000)])
        pattern = f'(pJC|pYTK|pBLS|pJBK|pJZ)({indices})'
        return pattern

    def get_strain_name(self):
        """
        Return the strain name pattern
        `JPC` or `BLS` + 3 digits
        """
        indices = '|'.join([str(num).zfill(3) for num in range(0, 1000)])
        pattern = f'(BLS|JPC)({indices})'
        return pattern

    def get_tet_concn(self):
        """
        Return the pattern used for naming data that
        has been treated with tetracycline typically
        to inhibit translation of proteasome substrates

        e.g. '250uM-Tet'
        """
        numbers = '|'.join([str(num).zfill(3) for num in range(0, 1000)])
        pattern = f"({numbers})(uM-Tet)"
        return pattern

    def get_genotype(self):
        """
        Return a pattern with a collection of typical
        genetic backgrounds found in experiments

        Should just make a .csv from which to draw these values
        """
        genotypes = ['BY4741',
                     'CMY3465',
                     'rpn4d',
                     'ubr2d',
                     'BY4742',
                     'pdr5d']
        genotypes = '|'.join(genotypes)
        pattern  = f'({genotypes})'
        return pattern

    def get_estradiol_concn(self):
        """
        Return the pattern used for naming data that
        has been treated with beta-estradiol, typcially
        to drive expression via pZ4EV

        e.g. '200nM-Estradiol'
        """
        numbers = '|'.join([str(num).zfill(3) for num in range(0, 1000)])
        pattern = f"({numbers})(nM-Estradiol)"
        return pattern

    def get_culture_condition(self):
        """
        Continuous OD or left at room temperature.
        Likely one-off from 20210421 experiment
        """
        conditions = '|'.join(['constantOD', 'RT'])
        pattern = f'({conditions})'
        return pattern

    def get_timepoint_mins(self):
        """
        Pattern for doing timecourse imaging of cells
        on coverslip

        e.g. "time010"
        """
        numbers = '|'.join([str(num).zfill(3) for num in range(0, 1000)])
        pattern = f"(time)({numbers})"
        return pattern

    def get_clone_number(self):
        """
        E.g. 'clone1'
        """
        numbers = '|'.join([str(num) for num in range(1,20)])
        pattern = f"(clone)({numbers})"
        pattern = f"(C)({numbers})"
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

    def get_cell_index(self):
        """
        Return pattern for finding cell index string
        like "cell001" in filenames and other kinds of
        text
        """
        numbers = '|'.join([str(num).zfill(3) for num in range(0, 1000)])
        pattern = f"(cell)({numbers})"
        return pattern

    def get_channel_name(self):
        """
        Return pattern for finding a channel name
        in a string. e.g. "rfp", "yfp" etc.
        """
        channels = [
            'bf',
            'gfp',
            'yfp',
            'rfp',
            'mko',
            'e2c'
        ]
        channels_str = '|'.join(channels)
        pattern = f"({channels_str})"
        return pattern

patterns = Patterns()
patterns_list = [patterns.genotype,
                 patterns.strain_name,
                 patterns.plasmid_name,
                 patterns.date,
                 patterns.tet_concn,
                 patterns.estradiol_concn,
                 patterns.clone_number]

col_names = ['background',
             'strain_name',
             'plasmid',
             'date',
             'tet_concn',
             'estradiol_concn',
             'clone_number']

patterns_dict = dict(zip(col_names, patterns_list))

# Documentation of final replicates used for manuscript etc.
compartments_dict = {
    'wt_rkk': [
        '20221019_byc_JPC228_UBL-YFP-Su9_full',
        '20221028_byc_JPC228_DsRed-P2A-UBL-YFP-Su9_3x-int',
        '20221104_byc_JPC220_R-YFP-Su9_Hsp104-mCherry',
        '20221118_byc_JPC220_R-YFP-Su9_Hsp104-mCherry',
        '20221201_byc_JPC220_R-YFP-Su9_HSP104-mCherry'
    ],
    'wt_ubl': [
        '20230126_byc_JPC228_DsRed-P2A-UBL-YFP-Su9_BY4741',
        '20230201_byc_JPC228-20230127-int_UBL-YFP-Su9x3_BY4741'
    ],
    'wt_odc': [
        '20220610_byc_JPC122_-YFP-ODC_two-chase',
        '20220616_byc_JPC122_YFP-ODC_BY4741_two-chase',
        '20220211_byc_JPC122_YFP-ODC(47)_BY4741',
        '20220825_byc_JPC122_two-chase_1000xetOH',
        '20221028_byc_JPC122_YFP-ODC_x3-int',
        '20221201_byc_JPC122_YFP-ODC'
    ],
    'mub1_rkk': [
        '20230822_byc_JPC274_R-YFP-Su9_in_mub1d',
        '20230518_byc_JPC274_R-YFP-Su9_mub1d',
        '20230810_byc_JPC274_R-YFP-Su9_mub1d'
    ],
    'ubr2_rkk': [
        '20230810_byc_JPC262_R-YFP-Su9_ubr2d',
        '20230714_byc_JPC262_R-YFP-Su9_ubr2d',
        '20230822_byc_JPC262_R-YFP-Su9_in_ubr2d'
    ],   
    'ubr2_ubl': [
        '20230518_byc_JPC259_DsRed-P2A-UBL-YFP-Su9_ubr2d'
    ],
    'rpn4_rkk': [
        '20230303_byc_JPC263_R-YFP-Su9_rpn4d',
        '20230609_byc_JPC277_R-YFP-Su9_in_rpn4d_cue5d'
    ],
    'rpn4_odc': [
        '20230224_byc_JPC261_YFP-ODC_rpn4d'
    ],
    'rpn4_ubl': [
        '20230201_byc_JPC226-20230127-int_UBL-YFP-Su9x3_rpn4d'
    ],
    'wt_dsred': [
        '20230201_byc_JPC228-20230127-int_UBL-YFP-Su9x3_BY4741_rfp'
    ],
    'hsp104OE_rkk': [
        '20220915_byc_JPC193_R-YFP-Su9_HSP104-OE',
        '20221215_byc_JPC193_R-YFP-Su9_HSP104-OE'
    ],
    'ssa1OE_rkk': [
        '20221104_byc_JPC199_R-YFP-Su9_SSA1-OE'
    ],
    'ssa2OE_rkk': [
        '20220610_byc_JPC196_SSA2-OE_in_JPC121_R-YFP-Su9_two-chase'
    ],
    'ubr1OE_rkk': [
        '20220118_byc_JPC136_pJC359_Ubr1_OE_w_pJC495',
        '20230324_byc_JPC136_R-YFP-Su9_Ubr1-OE'
    ]
}