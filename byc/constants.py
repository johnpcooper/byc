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
    source_path = os.path.abspath(os.path.join(dirname, '..'))
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
                     'path',
                     'xy',
                     'cell_index',
                     'death',
                     'sen_start',
                     'sen_end',
                     'note',
                     'late_daughter_shape',
                     'Pos',
                     'dist_from_sen']