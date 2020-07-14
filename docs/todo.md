## byc/python

`byc.dataBase` updates are needed. Need to make function to replace existing database entry for existing

Need to make standard_analysis.set_processed_traces integrated with byc.database.

* Make it add new trace dfs to database if they're not already in there
* Name experiments using columns in the `master_index_df` (date, type, plasmid etc.)
* Populate database fields for chase index etc.
* Add column in database for path to `master_cells_index.csv`. 

## Imagej

Write a macro in imagej that installs all my scripts and macros. 

1. Record macro for installing a script 

2. Loop the macro over a list of script paths listed in imagejpc_reqs.txt or something