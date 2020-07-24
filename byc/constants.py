import os

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

default_channels = ['yfp', 'dsred']                     

class Patterns(object):

    def __init__(self):

        self.master_index_file = self.get_master_index_filename()
        self.date = self.get_date()

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

patterns = Patterns()