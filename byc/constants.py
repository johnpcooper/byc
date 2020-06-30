def package_path(**kwargs):
    """
    Return the path to the local installation
    of byc
    """
    import byc
    file = byc.__file__
    try: # this will work on linux based OS
        end = file.rindex('/')
    except ValueError: # this will work on windows
        end = file.rindex('\\')
    package_path = file[:end+1]
    
    return package_path
    
package_path = package_path()

data_paths = {'cell_trace_db_path': f'{package_path}cell_trace_database.csv'}