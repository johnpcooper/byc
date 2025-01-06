import os
import numpy as np

from byc import constants, database

trace_df, fits_df, buds_df = database.read_in_trace_fits_buds_dfs()
unique_measdices = fits_df.measdex.unique()
bud_dfs = []
for i, measdex in enumerate(unique_measdices):
    print(f'Annotating measurement {i+1} of {len(unique_measdices)}')
    cell_bud_df = database.annotate_bud_info_on_fits_df(measdex, fits_df, buds_df)
    # function will return None if it fails to annotate
    if cell_bud_df is None:
        continue
    else:
        bud_dfs.append(cell_bud_df)
# Annotate information that can be done across the whole fits_df
sep_border = 5
fits_df.loc[fits_df.dist_from_sen<sep_border, 'post_sep'] = True
fits_df.loc[fits_df.dist_from_sen>sep_border, 'post_sep'] = False
n_long_thresh = 1
fits_df.loc[fits_df.n_long_daughters>=n_long_thresh, 'produced_elongated_daughters'] = True
fits_df.loc[fits_df.n_long_daughters==0, 'produced_elongated_daughters'] = False
fits_df.loc[fits_df.dist_from_sen>5, 'post_sep'] = False
fits_df.loc[fits_df.dist_from_sen<5, 'post_sep'] = True

savepath = os.path.join(constants.byc_data_dir, 'meta/fits_df_annotated.csv')
fits_df.to_csv(savepath, index=False)
print(f'Saved file at\n{savepath}')
savepath = os.path.join(constants.byc_data_dir, 'meta/buds_df_annotated.csv')
buds_df.to_csv(savepath, index=False)
print(f'Saved file at\n{savepath}')
# Use data from above to annotate the buds_df
print(f'Setting compartment_bud_id on buds_df')
# This is already done on fits_df via database.annotate_bud_info_on_fits_df
# database.set_bud_id(buds_df)
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

# Save annotated data. Ultimately will get rid of the "annotated" part
# in the filename once I know this work

savepath = os.path.join(constants.byc_data_dir, 'meta/buds_df_annotated.csv')
buds_df.to_csv(savepath, index=False)
print(f'Saved file at\n{savepath}')