import os
import sys
import numpy as np
import pandas as pd

from byc import constants, database

def parse_args():
    if len(sys.argv) > 1:
        compartments_string = sys.argv[1]
        new_compartments = compartments_string.split('|')
        return new_compartments
    
new_compartments = parse_args()
print(f'Adding new compartments:')
for compartment in new_compartments:
    print(compartment)

# Read in existing database dataframes
extant_traces_df, extant_fits_df, extant_buds_df = database.read_in_trace_fits_buds_dfs()

if len(new_compartments) > 0:
    # Only read in, annotate, and compile raw data specified in the argv
    # for running this refresh_database.py script
    fits_df, traces_df = database.get_byc_fits_df(
        return_traces_df=True,
        compartment_names=new_compartments
        )
elif len(new_compartments) == 0:
    # Instantiate fits, traces, and buds databases for all 
    # compartments found in constants.byc_data_dir
    fits_df, traces_df = database.get_byc_fits_df(
        return_traces_df=True,
        compartment_names=[]
        )
# xy data are needed for reading in and annotating bud information
# in database.BudDataBase(fits_df)
fits_df = fits_df[~fits_df.xy.isna()]
traces_df = traces_df[~traces_df.xy.isna()]
fits_df = database.label_from_manual_byc_index(fits_df)
traces_df = database.label_from_manual_byc_index(traces_df)

for df in [fits_df, traces_df]:
    df.loc[:, 'bud_rois_path'] = df.loc[:, 'bud_roi_set_path']


BudDB = database.BudDataBase(fits_df[:])
buds_df = BudDB.buddf
# Not sure where 'esca pe' is coming from but whatever
buds_df.loc[buds_df.end_event_type=='esca pe', 'end_event_type'] = 'escape'
# Make sure that both bud roi sets are labelled as escaped. There are cases
# where I forget to mark cell as having escaped when I save the bud roi set
# again and generate the bud_rois_df.csv
for bud_serial in buds_df.bud_serial.unique():
    if 'escape' in buds_df.loc[buds_df.bud_serial==bud_serial, 'end_event_type'].values:
        buds_df.loc[buds_df.bud_serial==bud_serial, 'end_event_type'] = 'escape'
buds_df = database.label_from_manual_byc_index(buds_df)

# Screen out any data that was calculated erroneously
buds_df.loc[buds_df.dist_from_sen<1, 'cycle_duration_hrs'] = np.nan

# Add a "measdex" which will be unique for every measurement
# made regardless of whether it's one of multiple on the same cell
for df in [fits_df, traces_df, buds_df]:
    # Separates unique deacy rate measurements
    df.loc[:, 'measdex'] = df.compartment_name.str.cat(df.cell_index.astype(int).astype(str), sep='-cell')
# Label the buds df and trace and fits dfs with their budding information.
# This will allow us to tell which measurements came from the same cell
# Separates unique cells (some have multiple decay rate measurements)
buds_df.loc[:, 'celldex'] = buds_df.compartment_name.str.cat(buds_df.bud_serial, sep='-bud_hours:')

unique_measdices = fits_df.measdex.unique()
bud_dfs = []
# Annotate cycle_duration_at_chase etc.
for i, measdex in enumerate(unique_measdices):
    print(f'Annotating measurement {i+1} of {len(unique_measdices)}')
    cell_bud_df = database.annotate_bud_info_on_fits_df(measdex, fits_df, buds_df)
    # function will return None if it fails to annotate
    if cell_bud_df is None:
        continue
    else:
        bud_dfs.append(cell_bud_df)
# Annotate information that can be done across all rows of the fits_df
sep_border = 5
fits_df.loc[fits_df.dist_from_sen<sep_border, 'post_sep'] = True
fits_df.loc[fits_df.dist_from_sen>sep_border, 'post_sep'] = False
n_long_thresh = 1
fits_df.loc[fits_df.n_long_daughters>=n_long_thresh, 'produced_elongated_daughters'] = True
fits_df.loc[fits_df.n_long_daughters==0, 'produced_elongated_daughters'] = False
fits_df.loc[fits_df.dist_from_sen>5, 'post_sep'] = False
fits_df.loc[fits_df.dist_from_sen<5, 'post_sep'] = True

# Use data from above to annotate the buds_df
# Compartment_bud_id gets set on buds_df when we run
# database.annotate_bud_info_on_fits_df() above
for i, compartment_bud_id in enumerate(buds_df.compartment_bud_id.unique()):
    if compartment_bud_id in fits_df.compartment_bud_id.unique():
        pass
    else:
        print(f'Cell not found in fits_df\n{compartment_bud_id}')
        continue
    print(f'Annotating bud set {i+1} of {len(buds_df.compartment_bud_id.unique())} with fits_df information')
    fits_df_slice = fits_df[fits_df.compartment_bud_id==compartment_bud_id]
    # Annotate the buds_df
    # aging mode data here will be the same for both measurements
    # in the fits_df because the data are derived from the bud roi set
    # for the cell in both those measurements
    buds_df.loc[buds_df.compartment_bud_id==compartment_bud_id, 'produced_elongated_daughters'] = fits_df_slice.produced_elongated_daughters.values[0]
    buds_df.loc[buds_df.compartment_bud_id==compartment_bud_id, 'n_long_daughters'] = fits_df_slice.n_long_daughters.values[0]
    # Sometimes I forget to annotate end_event_type correctly when I 
    # re-save bud ROI sets. So if any of them have escape annotated, set
    # that as the end event type for all belonging to that cell
    if 'escape' in buds_df.loc[buds_df.compartment_bud_id==compartment_bud_id, 'end_event_type'].values:
        buds_df.loc[buds_df.compartment_bud_id==compartment_bud_id, 'end_event_type'] = 'escape'
    # Whether the cell contains an aggregate at chase
    # is only accurately annotated in crop roi sets, so if there were multiple
    # measurements for a single cell, consider it to contain an aggregate in 
    # the buds df if there was an aggregate in one of those measurements
    if len(fits_df_slice)==1:
        buds_df.loc[buds_df.compartment_bud_id==compartment_bud_id, 'contains_aggregate'] = fits_df_slice.contains_aggregate.values[0]
    else:
        # If there are multiple measurements from the same cell, then if
        # there was an aggregate present at one of those measurements, 
        # annotate that in buds_df
        contains_aggregate_list = fits_df_slice.contains_aggregate.values
        if 'yes' in contains_aggregate_list:
            contains_aggregate = 'yes'
        elif np.nan in list(contains_aggregate_list):
            contains_aggregate = np.nan
        else:
            contains_aggregate = 'no'
        buds_df.loc[buds_df.compartment_bud_id==compartment_bud_id, 'contains_aggregate'] = contains_aggregate

buds_df.loc[buds_df.first_bud_frame<=18, 'observed_since_start'] = True
buds_df.loc[buds_df.first_bud_frame>18, 'observed_since_start'] = False

# Avoid adding duplicate cells 
new_fits_df = fits_df.loc[~fits_df.measdex.isin(extant_fits_df.measdex), :]
new_traces_df = traces_df.loc[~traces_df.measdex.isin(extant_traces_df.measdex), :]
new_buds_df = buds_df.loc[~buds_df.measdex.isin(extant_buds_df.measdex), :]
new_fits_df = pd.concat([extant_fits_df, new_fits_df])
new_traces_df = pd.concat([extant_traces_df, new_traces_df])
new_buds_df = pd.concat([extant_buds_df, new_buds_df])
# Write newly created databases to harddrive
database.write_trace_fits_buds_dfs(new_traces_df, new_fits_df, new_buds_df, save_non_gzip=True)