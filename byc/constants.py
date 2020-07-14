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

byc_data_dir = os.path.join(source_path, 'data')