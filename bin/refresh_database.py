import numpy as np
from byc import database

fits_df, traces_df = database.get_byc_fits_df(return_traces_df=True)
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
# Added binning for decay trace summaries
bin_borders_dict = {
    1: (0, 4),
    2: (6, 15),
    3: (15, np.Inf)
}

traces_df.loc[:, 'dist_from_sen_bin'] = np.nan

for key, borders in bin_borders_dict.items():
    binmask = traces_df.dist_from_sen.between(borders[0], borders[1], inclusive='right')
    traces_df.loc[binmask, 'dist_from_sen_bin'] = str(borders)

# Add indices to differentiate multiple measurements in same cell
for df in [fits_df, traces_df, buds_df]:
    # Separates unique deacy rate measurements
    measdex = df.compartment_name.str.cat(df.cell_index.map(str), sep='-cell')
    df.loc[:, 'measdex'] = measdex
# Label the buds df and trace and fits dfs with their budding information.
# This will allow us to tell which measurements came from the same cell
# Separates unique cells (some have multiple decay rate measurements)
buds_df.loc[:, 'celldex'] = buds_df.compartment_name.str.cat(buds_df.bud_serial, sep='-bud_hours:')

database.write_trace_fits_buds_dfs(traces_df, fits_df, buds_df, save_non_gzip=True)