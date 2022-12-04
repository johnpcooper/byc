import os
import sys

import pandas as pd
import numpy as np

from byc import constants, segmentation
from byc import standard_analysis as sa
from byc import plotting

if __name__=="__main__":
    # Pass either one or two arg variables

    # Argv[1] is compartmentname or the directory name
    # holding analysis of this experiment in 
    # constants.source_dir/data/<compartmentname>

    # argv[2] is the channels collected and needing to be 
    # analyzed in this experiment e.g. 'bf|yfp|rfp'

    # segmentation df is redundant to allmeasureddf
    # which contains all the segmentation information
    # per frame as well as mean intensity within cell
    # ROI per frame
    write_segmentation_df = False
    plot_segmentation_results = True

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
    elif len(sys.argv) == 3:
        compartmentname = sys.argv[1]
        # Extract experiment name from compartment name
        # that can be used to get absolute path to the 
        # compartment directory
        splitname = compartmentname.split('_')
        exptname = '_'.join(splitname[0:2])
        channels = str(sys.argv[2]).split(' ')
        # Need to strip non-alpha numerics because somehow
        # quote symbols end up in string when passing
        # 'ch1 ch2 ch3' as argv
        channels = [''.join(filter(str.isalnum, name)) for name in channels]
        print(f'Using channels {channels}')
        chase_frame_dict = {
            0: 14,
            148: 158
        }
    else:
        sys.exit("Please pass a compartment name argument variable")

    args = [exptname,
            compartmentname]

    kwargs = {
        'chase_frame_dict': chase_frame_dict,
        'channels': channels
    }
    print(f'Proceeding with compartment {compartmentname} from expt {exptname}')
    mdf = sa.create_and_annotate_mdf(*args, **kwargs)
    mdf.loc[:, 'channels_collected'] = ' '.join(channels)
    ds = sa.bycDataSet(mdf=mdf)
    # Now need to roll the stuff below into a single function
    # that reads in crop_roi_dfs and annotates fluorescence
    # 
    crop_roi_dfs = segmentation.make_cell_roi_dfs(mdf, bycdataset=ds)
    cellstacksdicts_dict = {}
    for channel in channels:
        cellstacks, cellstacksdict = segmentation.get_cell_stacks(crop_roi_dfs, channel_name=channel)
        cellstacksdicts_dict[channel] = cellstacksdict
    bfcellstacksdict = cellstacksdicts_dict['bf']
    # segmentation.get_cell_stacks() runs segmentation.cropped_stack_from_cellroidf
    # which adds relative x, y center coordinates according to the x and y
    # buffer size used in the crop. Thus, I concatanate after that function
    # has been run
    crop_rois_df = pd.concat(crop_roi_dfs, ignore_index=True)

    # Find radial intensity peaks and segment cell areas
    peak_idx_to_use=0
    allframesdfs = []
    for cell_index in crop_rois_df.cell_index.unique():
        args = [cell_index,
                crop_rois_df,
                bfcellstacksdict]

        kwargs = {
            'use_img_inverse': False,
            'use_constant_circle_roi': True,
            'default_radius_px': 9
        }
        allframesdf = segmentation.cell_stack_I_by_distance_and_theta(*args, **kwargs)
        # Filter outliers from within theta groups
        allframesdf = segmentation.filter_outlier_peaks(allframesdf)
        if plot_segmentation_results:
            try:
                plotting.save_segmentation_visualization(allframesdf, 'bf', bfcellstacksdict)
                # plotting.save_segmentation_visualization(allframesdf, 'bf', bfcellstacksdict, draw_outline=False)
                # plotting.save_segmentation_visualization(allframesdf, 'yfp', yfpcellstacksdict)
            except Exception as E:
                print(f'Failed to save cell vides with following exception:\n{E}')
        allframesdfs.append(allframesdf)
        
    alldf = pd.concat(allframesdfs, ignore_index = True)
    writedir = os.path.join(constants.byc_data_dir, crop_rois_df.compartment_reldir.iloc[0])
    writepath = os.path.join(writedir, 'segmentation_df.csv')
    if write_segmentation_df:
        print(f'Writing segmentation dataframe to\n{writepath}')
        alldf.to_csv(writepath, index=False)
        print(f'Wrote segmentation dataframe to\n{writepath}')
    else:
        print(f'Not writing segmentation dataframe')
    # Measure instensity inside segmented cell area
    allframesdfs_measured = []
    for cell_index in mdf.cell_index.unique():
        allframesdf = alldf[alldf.cell_index==cell_index]
        for channel, cellstacksdict in cellstacksdicts_dict.items():
            if type(list(cellstacksdict.keys())[0]) == str:
                measurement_stack = cellstacksdict[str(cell_index)]
            else:
                measurement_stack = cellstacksdict[cell_index]
            frame_cell_masks = []
            measurements = []
            for frame_idx in allframesdf.frame_rel.unique():
                mask = segmentation.get_frame_cell_mask(allframesdf,
                                        measurement_stack,
                                        frame_idx)
                frame_cell_masks.append(mask)
                measurement = np.mean(measurement_stack[int(frame_idx)][mask])
                measurements.append(measurement)
                allframesdf.loc[allframesdf.frame_rel==frame_idx, f'Mean_{channel}_auto'] = measurement
        allframesdfs_measured.append(allframesdf)
        
    allmeasureddf = pd.concat(allframesdfs_measured, ignore_index=True)
    # Split combined cell dataframes back into individual ones
    celldfs = [allmeasureddf[allmeasureddf.cell_index==cidx] for cidx in allmeasureddf.cell_index.unique()]
    
    path = os.path.join(writedir, f"{compartmentname}_alldf_measured.csv.gzip")
    print(f'Writing measured segmentation dataframe to\n{path}')
    allmeasureddf.to_csv(path, index=False, compression='gzip')
    print(f'Wrote all trace measurements df at \n{path}')
    # Aggregate dataset into just one entry per 
    # cell per frame (get rid of intensity vs. theta,
    # radial distance per frame)
    idx = ['cell_index',
        'frame_rel']
    table = pd.pivot_table(allmeasureddf, index=idx).reset_index()
    filename = f'{exptname}_alldf.csv'
    savepath = os.path.join(writedir, filename)
    table.to_csv(savepath, index=False)
    print(f'Saved trace table at\n{savepath}')