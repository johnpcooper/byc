import os
import sys

import pandas as pd
import numpy as np

from byc import constants, segmentation
from byc import standard_analysis as sa

if __name__=="__main__":

    if len(sys.argv) == 2:
        compartmentname = sys.argv[1]
        # Extract experiment name from compartment name
        # that can be used to get absolute path to the 
        # compartment directory
        splitname = compartmentname.split('_')
        exptname = '_'.join(splitname[0:2])
        channels = ['bf', 'yfp']
        chase_frame_dict = {
            0: 14,
            148: 158
        }
    else:
        print("Please pass a compartment name argument variable")
        quit()
    args = [exptname,
            compartmentname]

    kwargs = {
        'chase_frame_dict': chase_frame_dict
    }
    print(f'Proceeding with compartment {compartmentname} from expt {exptname}')
    mdf = sa.create_and_annotate_mdf(*args, **kwargs)
    ds = sa.bycDataSet(mdf=mdf)
    # Now need to roll the stuff below into a single function
    # that reads in crop_roi_dfs and annotates fluorescence
    # 
    crop_roi_dfs = segmentation.make_cell_roi_dfs(mdf, bycdataset=ds)

    chname = channels[0]
    bfcellstacks, bfcellstacksdict = segmentation.get_cell_stacks(crop_roi_dfs, channel_name=chname)
    chname = channels[1]
    yfpcellstacks, yfpcellstacksdict = segmentation.get_cell_stacks(crop_roi_dfs, channel_name=chname)
    # segmentation.get_cell_stacks() runs segmentation.cropped_stack_from_cellroidf
    # which adds relative x, y center coordinates according to the x and y
    # buffer size used in the crop. Thus, I concatanate after that function
    # has been run
    crop_rois_df = pd.concat(crop_roi_dfs, ignore_index=True)

    # Find radial intensity peaks and segment cell areas
    cellstacksdict = bfcellstacksdict
    peak_idx_to_use=0
    from byc.segmentation import find_radial_avg_intensity_peaks
    crop_rois_df = find_radial_avg_intensity_peaks(crop_rois_df,
                                                cellstacksdict,
                                                peak_idx_to_use=peak_idx_to_use)

    allframesdfs = []
    for cell_index in crop_rois_df.cell_index.unique():
        args = [cell_index,
                crop_rois_df,
                cellstacksdict]
        allframesdf = segmentation.cell_stack_I_by_distance_and_theta(*args)
        allframesdfs.append(allframesdf)
        
    alldf = pd.concat(allframesdfs, ignore_index = True)
    writedir = os.path.join(constants.byc_data_dir, crop_rois_df.compartment_reldir.iloc[0])
    writepath = os.path.join(writedir, 'segmentation_df.csv')
    print(f'Writing segmentation dataframe to\n{writepath}')
    alldf.to_csv(writepath, index=False)
    print(f'Wrote segmentation dataframe to\n{writepath}')
    # Measure instensity inside segmented cell area
    allframesdfs_measured = []
    for cell_index in mdf.cell_index.unique():
        allframesdf = alldf[alldf.cell_index==cell_index]
        if type(list(yfpcellstacksdict.keys())[cell_index]) == str:
            measurement_stack = yfpcellstacksdict[str(cell_index)]
        else:
            measurement_stack = yfpcellstacksdict[cell_index]
        frame_cell_masks = []
        measurements = []
        for frame_idx in allframesdf.frame_rel.unique():
            mask = segmentation.get_frame_cell_mask(allframesdf,
                                    measurement_stack,
                                    frame_idx)
            frame_cell_masks.append(mask)
            measurement = np.mean(measurement_stack[int(frame_idx)][mask])
            measurements.append(measurement)
            allframesdf.loc[allframesdf.frame_rel==frame_idx, 'Mean_yfp_auto'] = measurement
        allframesdfs_measured.append(allframesdf)
        
    allmeasureddf = pd.concat(allframesdfs_measured, ignore_index=True)
    # Split combined cell dataframes back into individual ones
    celldfs = [allmeasureddf[allmeasureddf.cell_index==cidx] for cidx in allmeasureddf.cell_index.unique()]
    
    path = os.path.join(writedir, f"{compartmentname}_alldf_measured.csv")
    allmeasureddf.to_csv(path, index=False)
    print(f'Wrote all trace measurements df at \n{path}')