import os

import tifffile as tf
import tkinter as tk
import tkinter.filedialog as dia
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from skimage.io import imsave, concatenate_images
from skimage.filters import threshold_otsu
from skimage import img_as_uint
import skimage
from scipy.signal import medfilt, find_peaks
from read_roi import read_roi_file, read_roi_zip

from byc import constants, utilities, files, standard_analysis

class Cell_Stack(object):
    """
    Upon instantiation, this class will ask the user to select a master index
    .csv file containing cell indices. The methods of the class then operate
    on rois etc. referenced by the master index. For now, the only public variable
    of this class is self.channel_names which is useful for saving stacks later.
    """
    
    def __init__(self, threshold_channel):
        
        # Set the master index df for this experiment
        self.master_cells_df = self.set_master_cells_df("Choose the .csv master index of this experiment")
        self.threshold_channel = threshold_channel
        
        # It may be useful at some point to add code here to run methods below upon instantiation.
    
    def set_fp(self, prompt):
        """
        Return the path to a file of interest. Call with the prompt you would like to 
        display in the dialog box.
        """

        # create the dialog box and set the fn
        root = tk.Tk()
        fp = dia.askopenfilename(parent=root, title=prompt)
        root.destroy() # very important to destroy the root object, otherwise python 
        # just keeps running in a dialog box

        return fp # return the path to the file you just selected
    
    def set_master_cells_df(self, prompt):
        """ 
        Return a dataframe read from the master cells index .csv
        """

        # define the path to the index .csv
        master_cells_fp = self.set_fp(prompt)
        # define the DataFrame for the master expt index
        master_cells_df = pd.read_csv(master_cells_fp)

        return master_cells_df

    def get_compartment_dir(self, master_cells_df, cell_index):
        # Works on older formatted master indexes
        if 'path' in master_cells_df.columns:
            compartment_dir = master_cells_df.path[cell_index]
        # Works for newer formatted master indexes where
        # compartment_dir and compartment_reldir are determined automatically
        # during addition of various cell roi types
        elif 'compartment_reldir' in master_cells_df.columns:
            compartment_reldir = master_cells_df.compartment_reldir[cell_index]
            compartment_dir = os.path.join(constants.byc_data_dir, compartment_reldir)
        else:
            print(f"No compartment directory in master index")
            compartment_dir = None

        return compartment_dir

    def set_cell_crop_roi_dfs(self, master_cells_df):
        """
        Return a list of DataFrames, one for each cell. The coordinates in each
        of these DataFrames will be used to crop from the image stacks in 
        set_cropped_cell_stack_list()
        """
    
        cell_crop_roi_dfs = []
        abs_i = 0
        for cell_index in master_cells_df.cell_index:
            # 
            compartment_dir = self.get_compartment_dir(master_cells_df, abs_i)
            expt_date = int(master_cells_df.date[abs_i])
            expt_type = master_cells_df.expt_type[abs_i]
            cell_rois_fp = f"{compartment_dir}\\{expt_date}_{expt_type}_cell{str(cell_index).zfill(3)}_crop_rois.zip"
            print(f'Looking for cell rois at {cell_rois_fp}')
            # cell_rois is an OrderedDict so I split the object up and put it in a dataframe 
            # to get relevant data out of it
            cell_rois = read_roi_zip(cell_rois_fp)
            roi_keys = list(cell_rois.keys())

            frames = [cell_rois[roi_key]['position'] - 1 for roi_key in roi_keys]

            x_lbs = [cell_rois[roi_key]['left'] for roi_key in roi_keys]
            x_ubs = [cell_rois[roi_key]['left'] + cell_rois[roi_key]['width'] for roi_key in roi_keys]
            
            # The y crop boundary settings is a bit weird to me because y increases as you move 
            # down the image, not up.
            y_ubs = [cell_rois[roi_key]['top'] + cell_rois[roi_key]['height'] for roi_key in roi_keys]            
            y_lbs = [cell_rois[roi_key]['top'] for roi_key in roi_keys]
            
            width = [cell_rois[roi_key]['width'] for roi_key in roi_keys]
            height = [cell_rois[roi_key]['height'] for roi_key in roi_keys]

            cell_rois_df = pd.DataFrame({'cell_index': cell_index,
                                         'frame' : frames,
                                         'x_lb' : x_lbs,
                                         'y_lb' : y_lbs,
                                         'x_ub' : x_ubs,
                                         'y_ub' : y_ubs,
                                         'width' : width,
                                         'height' : height})

            cell_crop_roi_dfs.append(cell_rois_df)

            abs_i += 1
        
        return cell_crop_roi_dfs
    
    def set_cell_channel_stacks(self, master_cells_df, cell_index):
        """
        Return a list of tf.imread() objects (one for each channel collected in expt)
        read according to data in the master_cells_df
        """ 
        compartment_dir = self.get_compartment_dir(master_cells_df, cell_index)
        expt_date = int(master_cells_df.date[cell_index])
        expt_type = master_cells_df.expt_type[cell_index]
        xy = str(int(master_cells_df.xy[cell_index]))
        # Should be the channel names used to create stacks as output
        # of process. Channel names are separated by spaces in master_cells_df
        self.channel_names = master_cells_df.channels_collected[cell_index].split()

        # Set a list of paths to the channel stacks
        channel_stacks_fps = [f"{compartment_dir}\\{expt_date}_{expt_type}_xy{xy.zfill(2)}_{channel_name}_stack.tif" for channel_name in self.channel_names]
        
        print(f"Paths to channel stacks for cell {cell_index} found based on master df")
        # Let the user know which stacks were read
        for fp in channel_stacks_fps:
            print(fp)
            
        # Set a list of the actual channel stack imread() objects
        cell_channel_stacks = [tf.imread(filepath) for filepath in channel_stacks_fps]
        
        return cell_channel_stacks
    
    def set_cell_cropped_stacks_dict(self, master_cells_df, cell_rois_dfs, cell_index):
        print(f"Setting cell cropped stacks dict for cell: {cell_index}")
        cell_rois_df = cell_rois_dfs[cell_index]
        channel_stacks = self.set_cell_channel_stacks(master_cells_df, cell_index)
        
        cell_cropped_channels_list = []        
        for channel_stack in channel_stacks:

            whole_fov_shape = channel_stack[0].shape
            whole_fov_half_width = whole_fov_shape[1] / 2

            cell_crop_stack = []

            # Iterate over each frame number that will be in the final cropped stack
            j = 0 # j is the index of the currently used crop roi
            for frame_index in range(cell_rois_df.frame.min(), cell_rois_df.frame.max()+1):
                # if the current frame_index has an roi assigned, use
                # that roi to crop for that frame
                if frame_index in cell_rois_df.frame.values:

                    print(f"Frame {frame_index} has an roi assigned")

                    y_lb = cell_rois_df.y_lb[j]
                    y_ub = cell_rois_df.y_ub[j]
                    x_lb = cell_rois_df.x_lb[j]
                    x_ub = cell_rois_df.x_ub[j]

                    print(f"Cropping cell {cell_index} frame {frame_index} with parameters x{x_lb},{x_ub} and y{y_lb},{y_ub}")

                    j += 1

                    source_image = channel_stack[frame_index]
                    cropped_image = source_image[y_lb: y_ub, x_lb: x_ub]

                    cell_crop_stack.append(cropped_image)

                elif frame_index not in cell_rois_df.frame.values:

                    y_lb = cell_rois_df.y_lb[j]
                    y_ub = cell_rois_df.y_ub[j]
                    x_lb = cell_rois_df.x_lb[j]
                    x_ub = cell_rois_df.x_ub[j]

                    source_image = channel_stack[frame_index]
                    cropped_image = source_image[y_lb: y_ub, x_lb: x_ub]

                    print(f"Cropping cell {cell_index} frame {frame_index} with parameters x: ({x_lb},{x_ub}) and y: ({y_lb},{y_ub})")

                    cell_crop_stack.append(cropped_image)

            cell_cropped_channels_list.append(cell_crop_stack)
            
        self.cell_cropped_channels_dict = dict(zip(self.channel_names, cell_cropped_channels_list))
        
        return self.cell_cropped_channels_dict
    
    def add_cell_otsu_thresholded_stack(self):
        
        """ Return nothing. Add the otsu thresholded cropped channel stack to the 
            cell_cropped_channels_ditc, this channel is specified on 
            instantiation of the Cell_Stacks class. """
        
        print("Applying otsu threshold")
        
        # Set cell_crop_stack to the stack specified on Cell_Stack instantiation
        # as the channel to use for thresholding
        cell_crop_stack = self.cell_cropped_channels_dict[self.threshold_channel]
        
        thresholded_stack = []

        for i in range(0, len(cell_crop_stack)):

            image = cell_crop_stack[i]
            threshold_value = threshold_otsu(image)
            binary = image > threshold_value

            # add fresh binary image to stack
            thresholded_stack.append(img_as_uint(binary))
            
        thresh_key = f'{self.threshold_channel}_otsu'
        self.cell_cropped_channels_dict[thresh_key] = thresholded_stack
    
    def set_resized_cell_cropped_channels_dict(self, cell_cropped_channels_dict):
        
        """ Return a dictionary of channel stacks resized so that each frame of 
            the stack is the same dimensions as the largest frame. This is necessary
            for saving the cell stack. """
        
        channel_keys = list(cell_cropped_channels_dict.keys())
        # Set the cropped stack that will be used as a reference for shape for this cell
        cropped_stack = cell_cropped_channels_dict[channel_keys[0]]

        # Set an array of shapes for each frame in the cell crop
        frame_shapes_list = [frame.shape for frame in cropped_stack]
        frame_shapes_array = np.array(frame_shapes_list, dtype='uint16')

        # Set max height in each dimension (y and x or height and width)
        max_height = frame_shapes_array[:, 0].max()
        print(f'max height = {max_height}')
        max_width = frame_shapes_array[:, 1].max()
        print(f'max width = {max_width}')

        # Set two lists of offsets telling you how each frame differs
        # from the shape of the max frame 
        frame_height_offsets = []
        for frame_height in frame_shapes_array[:, 0]:

            height_offset = max_height - frame_height
            # print(f"max height {max_height} - frame height {frame_height} = offset {height_offset}")
            frame_height_offsets.append(height_offset)

        frame_width_offsets = []
        for frame_width in frame_shapes_array[:, 1]:

            width_offset = max_width - frame_width
            frame_width_offsets.append(width_offset)

        resized_channels = []
        for channel_key in channel_keys:
            # iterate through the cropped stack for each channel
            cropped_stack = cell_cropped_channels_dict[channel_key]

            # Iterate through each frame of the cropped stack and resize it
            # to whatever max dimensions are for this cell stack
            resized_cropped_stack = []
            for i in range(0, len(cropped_stack)):
                h_offset = frame_height_offsets[i]
                w_offset = frame_width_offsets[i]
                frame = cropped_stack[i].copy()

                h_filler_value = frame[-1, :].min()
                h_filler_array = np.full((h_offset, frame.shape[1]), h_filler_value, dtype='uint16')
                #print(f"h_filler_array shape {(h_filler_array.shape)}")
                resized_frame = np.append(frame, h_filler_array, axis=0)
                # print(f"Shape of resized frame after adding h_filler_array {resized_frame.shape} and height offest {h_offset}")

                w_filler_value = frame[:, -1].min()
                w_filler_array = np.full((resized_frame.shape[0], w_offset), w_filler_value, dtype='uint16')
                resized_frame = np.append(resized_frame, w_filler_array, axis=1)

                resized_cropped_stack.append(resized_frame)

            resized_channels.append(resized_cropped_stack)
        print(channel_keys)
        cell_resized_channels_dict = dict(zip(channel_keys, resized_channels))
        return cell_resized_channels_dict
    
def save_cell_stacks():
    
    threshold_channel = input("Choose channel to use for thresholding: ")
    cs = Cell_Stack(threshold_channel=threshold_channel)
    try:
        cell_indices = cs.master_cells_df.sub_coord
    except:
        cell_indices = cs.master_cells_df.cell_index
    abs_i = 0
    for cell_index in cell_indices:
        print('Cell index = ', cell_index) 
        # Run all the cropping and processing methods of Cell_Stack on the cell
        # with cell_index
        cell_rois_dfs = cs.set_cell_crop_roi_dfs(cs.master_cells_df)
        cell_cropped_channels_dict = cs.set_cell_cropped_stacks_dict(cs.master_cells_df, cell_rois_dfs, abs_i)
        cs.add_cell_otsu_thresholded_stack()
        resized_channels_dict = cs.set_resized_cell_cropped_channels_dict(cs.cell_cropped_channels_dict)

        # Save each stack (Note: compartment_dir is the dir holding all
        # xy FOVs for that flow compartment and therefore isolated genetic
        # + environmental etc. condition)
        # Need to use absolute i so the right compartment_dir etc. is
        # selected when there are multiple compartments in a master_index
        # and thus cell_index restarts
        compartment_dir = cs.get_compartment_dir(cs.master_cells_df, abs_i)
        expt_date = int(cs.master_cells_df.date[abs_i])
        expt_type = cs.master_cells_df.expt_type[abs_i]
        xy = str(int(cs.master_cells_df.xy[abs_i]))

        for channel_name, stack in resized_channels_dict.items():
            filename = f'{expt_date}_{expt_type}_xy{str(xy).zfill(2)}_cell{str(cell_index).zfill(3)}_{channel_name}_stack.tif'
            save_path = f'{compartment_dir}//{filename}'
            try:
                imsave(save_path, concatenate_images(stack))
            except:
                print(f"Could not save stacks for cell {cell_index}, img dims must agree")
        abs_i += 1

def fill_crop_roi_df(cell_index, mdf, dataset, bf_stack=None):
    """
    Using the entry for cell <cell_index> in master index dataframe (<mdf>), 
    annotate the mdf with the x, y coordinates for the center of the crop
    roi for each frame, add fluorescence and brightfiled intensity measurements

    <dataset> Is an instantiation of standard_analysis.bycDataSet(mdf=<mdf>)

    Return the filled crop_roi_df
    """
    cell_row = mdf[mdf.cell_index==cell_index]
    # Read in various files
    rois_path = cell_row.crop_rois_path.iloc[0]
    bf_stack_path = cell_row.bf_stack_path.iloc[0]
    if bf_stack is None:
        bf_stack = skimage.io.imread(bf_stack_path)
    rois_df = files.read_rectangular_rois_as_df(rois_path)
    rois_df.loc[:, 'frame'] = rois_df.position - 1
    first_frame = rois_df.frame.min()
    final_frame = rois_df.frame.max()

    recent_frame = first_frame
    filled_rois_df = pd.DataFrame()

    stack_frames = range(first_frame, final_frame+1)

    i=0
    for frame in stack_frames:
        if frame in list(rois_df.frame):
            recent_frame = frame
            for col in rois_df.columns:
                filled_rois_df.loc[i, col] = rois_df.loc[rois_df.frame==frame, col].iloc[0]
        else:
            for col in rois_df.columns:
                filled_rois_df.loc[i, col] = rois_df.loc[rois_df.frame==recent_frame, col].iloc[0]
        # Change position and frame from parent values to actualy
        # values
        filled_rois_df.loc[i, 'frame'] = frame
        filled_rois_df.loc[i, 'position'] = frame + 1
        i += 1

    # Get center of crop roi and add it to the filled roi df
    x_center = filled_rois_df.left + filled_rois_df.width/2
    y_center = filled_rois_df.top + filled_rois_df.height/2
    filled_rois_df.loc[:, 'x_center'] = x_center
    filled_rois_df.loc[:, 'y_center'] = y_center
    filled_rois_df.loc[:, 'bf_stack_path'] = bf_stack_path
    filled_rois_df.loc[:, 'frame_rel'] = filled_rois_df.frame - first_frame
    # Add fluorescence and brightfild intensity measurements etc.
    alldf = pd.concat(dataset.cell_trace_dfs, ignore_index=True)
    cell_trace_df = alldf[alldf.cell_index==cell_index].reset_index()
    for col in list(cell_trace_df.columns):
        # Don't overwrite existing values in the
        # crop_roi_df
        if col not in list(filled_rois_df.columns):
            filled_rois_df.loc[:, col] = np.nan
            filled_rois_df.loc[:, col] = cell_trace_df.loc[:, col]
    # Add all information from this cell's entry in the
    # master index dataframe 
    for col in cell_row.columns:
        filled_rois_df.loc[:, col] = cell_row[col].iloc[0]
    
    return filled_rois_df

### Following two functions need to be rolled into
### standard_analysis.bycDataSet.make_cell_trace_dfs()

def find_channel_name(meas_rois_path, crop_rois_df):
    """
    Find the light channel name in the <meas_rois_path>.
    Typically 'bf', 'yfp', or 'rfp' etc. <meas_rois_path>
    doesn't have to necessarily be the path to measurement
    (cell outline) rois
    
    Return the name as a string
    """
    channels_collected = crop_rois_df.channels_collected.iloc[0]
    channels_collected = channels_collected.split()
    found = False
    for channel in channels_collected:
        query = f'_{channel}_stack'
        if query in meas_rois_path:
            found = True
            return channel
    if found == False:
        print(f'No channel name found in <meas_rois_path>:/n{meas_rois_path}')

def get_channel_dfs(crop_rois_df):
    """
    Find the channel intensity data measured using
    the measurement (cell outline) rois curated by
    the user
    
    This function does basically the same thing as:
    
        dataset = byc.DataSet(mdf=<mdf>)
        dfs = dataset.cell_trace_dfs()
    
    But is simpler. I should replace the
    byc.DataSet.make_cell_trace_dfs() function with it
    and also use it in 
    
    Return found channel_dfs in a list
    """
    meas_rois_path = crop_rois_df.outline_rois_path.iloc[0]
    to_replace = '_measurement_rois.zip'
    replacement = '.tif'
    channel_df_path_temp = meas_rois_path.replace(to_replace, replacement)
    channel_dfs = []
    # The measurement ROIs .zip was named using the path to the image
    # that was active when the measurement (cell outline) ROIs were
    # ceated using save_cell_roi_set.py (fiji plugin from imagejpc)
    source_channel_name = find_channel_name(meas_rois_path, crop_rois_df)
    for channel in crop_rois_df.channels_collected.iloc[0].split():
        if channel == source_channel_name:
            channel_df_path = channel_df_path_temp.replace('.tif', '.csv')
        else:
            channel_df_path = channel_df_path_temp.replace(source_channel_name, channel)
            channel_df_path = channel_df_path.replace('.tif', '.csv')
        channel_df = pd.read_csv(channel_df_path)
        # Annotate which channel the channelescence (light) trace
        # data came from
        newcols = [f'{channel}_{name}' for name in channel_df.columns]
        channel_df.columns = newcols
        channel_dfs.append(channel_df)
        print(f'Added channel df to {channel} to <channel_dfs> from:/n{channel_df_path}')

def make_cell_roi_dfs(mdf, bycdataset=None, return_crop_stacks=False):
    """
    For each cell in the <mdf>, read its crop rois into a dataframe,
    annotate the dataframe with x, y coordinates of crop ROI center,
    annotate with light intensity measurements 
    """
    if bycdataset is None:
        bycdataset = standard_analysis.bycDataSet(mdf=mdf)
        print(f'Read in data set using provided master index')
    crop_roi_dfs = []
    
    for i, cell_index in enumerate(list(mdf.cell_index.unique())):
        # Only read in a new brightfield stack if it's different from
        # the brightfield stack for the previous cell
        print(f'Making roi dataframe for cell {cell_index}')
        if i == 0:
            bf_stack_path = mdf.loc[mdf.cell_index==cell_index, 'bf_stack_path'].iloc[0]
            prev_bf_stack_path = bf_stack_path
            bf_stack = skimage.io.imread(bf_stack_path)
            print(f'Read in brightfield image with shape {bf_stack.shape} from/{bf_stack_path}')
        else:
            new_bf_stack_path = mdf.loc[mdf.cell_index==cell_index, 'bf_stack_path'].iloc[0]
            if prev_bf_stack_path == new_bf_stack_path:
                print('No new bf stack needed')
            else:
                bf_stack = skimage.io.imread(new_bf_stack_path)
                prev_bf_stack_path = new_bf_stack_path
                print(f'Read in brightfield image with shape {bf_stack.shape} from/{new_bf_stack_path}')
        # segmentation.fill_crop_roi_df annotates center of crop rois on each
        # frame and will soon include original measured fluorescence
        crop_rois_df = fill_crop_roi_df(cell_index, mdf, bycdataset, bf_stack=bf_stack)
        crop_roi_dfs.append(crop_rois_df)

    return crop_roi_dfs


def cropped_stack_from_cellroidf(cellroidf, source_stack=None, **kwargs):
    """
    Return cropped_frames, cellroidf where cropped_frames is cropped
    using center coordinates found in cellroidf, and cellroidf has
    been annotated with new x and y center coordinates (x_center_rel,
    y_center_rel) relative to their newly cropped frame
    """
    y_buffer = kwargs.get('y_buffer', 20)
    x_buffer = kwargs.get('x_buffer', 50)
    # Not using the timestamp and framestamp kwargs for now,
    # but would like to include options to timestamp
    timestamp = kwargs.get('timestamp', False)
    framestamp = kwargs.get('framestamp', False)

    cellroidf.loc[:, 'x_center_rel'] = np.nan
    cellroidf.loc[:, 'y_center_rel'] = np.nan
    # Read in the brightfield stack this cell was cropped
    # from if no source stack provided during function call
    if source_stack is None:
        source_stack_path = cellroidf.bf_stack_path.iloc[0]
        source_stack = skimage.io.imread(source_stack_path)
    cropped_frames = []
    for frame_idx in cellroidf.frame:
        frame = source_stack[int(frame_idx)]
        # 
        x_center = int(cellroidf.loc[cellroidf.frame==frame_idx, 'x_center'])
        y_center = int(cellroidf.loc[cellroidf.frame==frame_idx, 'y_center'])
        x_upper = x_center + x_buffer
        x_lower = x_center - x_buffer
        y_upper = y_center + y_buffer
        y_lower = y_center - y_buffer
        # Annotate relative center x and y base on limits
        # defined above
        cellroidf.loc[cellroidf.frame==frame_idx, 'x_center_rel'] = x_center - x_lower
        cellroidf.loc[cellroidf.frame==frame_idx, 'y_center_rel'] = y_center - y_lower
        # Make iure the crop boundaries aren't outside of the 
        # bounds of the frame
        if x_lower < 0:
            x_lower = 0
        if y_lower < 0:
            y_lower = 0
        if x_upper > frame.shape[1]:
            x_upper = frame.shape[1]
        if y_upper > frame.shape[0]:
            y_upper = frame.shape[1]
        cropped_frame = frame[y_lower:y_upper, x_lower:x_upper]
        cropped_frames.append(cropped_frame)

    return cropped_frames, cellroidf

def get_cell_crop_stack(cell_roi_df,
                        source_stack,
                        x_buffer=None,
                        y_buffer=None,
                        save_cell_stacks=True):
    """
    Crop the <source_stack> (typically brightfield or a fluorescent channel)
    using the x_center and y_center columns found in <cell_roi_df> with 
    <x_center>, <y_center> +/- <x_buffer>, <y_buffer>
    
    If x_buffer, y_buffer defaults to 50, 20 in segmentation.cropped_stack_from_cellroidf
    which is used below
    
    Return the cropped cell stack
    """
    crop_args = [cell_roi_df]
    crop_kwargs = {'source_stack': source_stack}
    if x_buffer is None:
        pass
    else:
        crop_kwargs['x_buffer'] = x_buffer
        pass
    if y_buffer is None:
        pass
    else:
        crop_kwargs['y_buffer'] = y_buffer
    
    cellstack, cell_roi_df = cropped_stack_from_cellroidf(*crop_args, **crop_kwargs)
    cell_index = cell_roi_df.cell_index.iloc[0]
    crop_rois_path = cell_roi_df.crop_roi_set_path.iloc[0]
    compartment_name = cell_roi_df.compartment_name.iloc[0]
    auto_compartment_name = f'{compartment_name}_auto'
    crop_rois_path_mod = crop_rois_path.replace(compartment_name, auto_compartment_name)
    # Make a directory in which to save this next gen-generated cell
    # crop stack
    if not os.path.exists(crop_rois_path_mod):
        os.mkdir(crop_rois_path_mod)
    stack_path = crop_rois_path_mod.replace('crop_rois.zip', 'crop.tif')
    cellstackarr = np.array(cellstack)
    if save_cell_stacks:
        skimage.io.imsave(stack_path, cellstackarr)
        print(f'Saved cell {cell_index} crop tif at\n{stack_path}')
        
    return cellstackarr

def get_cell_stacks(crop_roi_dfs,
                    channel_name='bf',
                    x_buffer=None,
                    y_buffer=None,
                    save_cell_stacks=True,
                    read_stack=False):

    """
    For each crop_roi_df in <crop_roi_dfs>, crop the cell tracking
    ROI from the source stack found for <channel_name>

    Return the list of cell_stacks and save each cell stack by default
    """
    
    cell_stacks = []
    cell_stacks_dict = {}
    for i, cellroidf in enumerate(crop_roi_dfs):
        cell_index = cellroidf.cell_index.iloc[0]
        # Only read in a new stack if the new cell has a different
        # Source xy than the previous cell
        path_colname = f'{channel_name}_stack_path'
        if i == 0:
            source_stack_path = cellroidf.loc[:, path_colname].iloc[0]
            prev_source_stack_path = source_stack_path
            source_stack = skimage.io.imread(source_stack_path)
            print(f'Read in {channel_name} image with shape {source_stack.shape} from/{source_stack_path}')
        else:
            new_source_stack_path = cellroidf.loc[:, path_colname].iloc[0]
            if prev_source_stack_path == new_source_stack_path:
                print('No new source stack needed')
            else:
                source_stack = skimage.io.imread(new_source_stack_path)
                prev_source_stack_path = new_source_stack_path
                print(f'Read in {channel_name} image with shape {source_stack.shape} from/{new_source_stack_path}')
        
        args = [cellroidf,
                source_stack]
        kwargs = {'x_buffer': x_buffer,
                  'y_buffer': y_buffer,
                  'save_cell_stacks': save_cell_stacks}

        cell_stack = get_cell_crop_stack(*args, **kwargs)
        cell_stacks_dict[str(cell_index)] = cell_stack
        cell_stacks.append(cell_stack)
    return cell_stacks, cell_stacks_dict


def radial_avg_intensity(img, x_center, y_center, **kwargs):
    """
    Return two 1D arrays: (rad, intensity)
    rad is a set of radial distances from the x and y
    center coordinates and intensity is average intensity
    in rings of each of those distances away from center
    """
    # Bin size is in pixels
    bin_size = kwargs.get('bin_size', 1)
    # Get image parameters
    a = img.shape[0]
    b = img.shape[1]
    # X is a matrix that gives the X coordinate at each pixel 
    # Y is a matrix that gives the Y coordinate at each pixel
    # X = 0 and Y = 0 is at the center point
    [X, Y] = np.meshgrid(np.arange(b) - x_center, np.arange(a) - y_center)
    # R is then the radial distance from the center to each pixel 
    # (a squared + b squared = c squared)
    R = np.sqrt(np.square(X) + np.square(Y))
    # rad is a 1d array that will be used as the
    # x axis later. It's a range that covers all
    # possible radial distances from our center pixel
    # to whatever the maximum distance was from center
    # found in R
    rad = np.arange(1, np.max(R), 1)
    # intensity will be filled with average
    # intensity at each radial distance from
    # center found in rad
    intensity = np.zeros(len(rad))
    index = 0
    for i in rad:
        # For each radial distance found in rad (the x axis), select
        # all the values in R (radial distances from center at each
        # coordinate) and store those coordinates in a boolean mask.
        # Then use that mask to select all the intensity values from
        # the original image and take the mean
        mask = (np.greater(R, i - bin_size) & np.less(R, i + bin_size))
        values = img[mask]
        intensity[index] = np.mean(values)
        index += 1
    return rad, intensity

def find_radial_avg_intensity_peaks(crop_rois_df, cellstacksdict):
    """
    After reading in a <crop_rois_df> and <cellstacksdict> as follows,
    pass them this find_radial_avg_intensity_peaks
    
        crop_roi_dfs = segmentation.make_cell_roi_dfs(mdf, bycdataset=ds)
        crop_rois_df = pd.concat(crop_roi_dfs, ignore_index=True)
        chname = 'bf'
        bfcellstacks, bfcellstacksdict = segmentation.get_cell_stacks(crop_roi_dfs, channel_name=chname)
        chname = 'yfp'
        yfpcellstacks, yfpcellstacksdict = segmentation.get_cell_stacks(crop_roi_dfs, channel_name=chname)
        
    Return crop_rois_df with radial avg intensity peaks etc. annotated

    This 'peak_avg_intensity_radius' can be used to draw a circle with center
    at center of cropped cell stack (in the cellstacksdict). The circle can then
    be used as an approximate cell outline ROI for that frame of that cell 
    crop stack
    """
    crop_rois_df.loc[:, 'circle_slice_radius'] = np.nan
    crop_rois_df.loc[:, 'intensity_avg'] = np.nan
    crop_rois_df.loc[:, 'peak_avg_intensity_radius'] = np.nan

    for cell_index in crop_rois_df.cell_index.unique():
        cell_bool = crop_rois_df.cell_index==cell_index
        cellroidf = crop_rois_df[cell_bool]
        cellstack = cellstacksdict[str(cell_index)]

        radii, intensities, peak_radii = [], [], []
        for frame_idx in cellroidf.frame_rel:
            frame_idx = int(frame_idx)
            x_center_rel = cellroidf.x_center_rel.iloc[frame_idx]
            y_center_rel = cellroidf.y_center_rel.iloc[frame_idx]
            rad, intensity = radial_avg_intensity(cellstack[frame_idx],
                                                  x_center_rel,
                                                  y_center_rel)


            peaks = find_peaks(intensity, height=np.array([300, None]))
            peak_radius = peaks[0][0]
            # Need to condense radius vs. intensity
            # into single row entries.
            joined_rad = '|'.join([str(rad_val) for rad_val in rad])
            radii.append(joined_rad)
            joined_intensity = '|'.join([str(intensity_val) for intensity_val in intensity])
            intensities.append(joined_intensity)
            peak_radii.append(peak_radius)

        crop_rois_df.loc[cell_bool, 'circle_slice_radius'] = radii
        crop_rois_df.loc[cell_bool, 'intensity_avg'] = intensities
        crop_rois_df.loc[cell_bool, 'peak_avg_intensity_radius'] = peak_radii
        
    return crop_rois_df

def get_polar_coordinates(img, x_center, y_center):
    """
    Return (R, theta): two arrays of the same shape as img.
    R is radial distance from x_center, y_center.
    Theta is angle from horizontal in radians.
    """
    # Get image parameters
    a = img.shape[0]
    b = img.shape[1]
    # X is a matrix that gives the X coordinate at each pixel 
    # Y is a matrix that gives the Y coordinate at each pixel
    # X = 0 and Y = 0 is now the center point
    [X, Y] = np.meshgrid(np.arange(b) - x_center, np.arange(a) - y_center)
    # R is then the radial distance from the center to each pixel 
    # (a squared + b squared = c squared)
    R = np.sqrt(np.square(X) + np.square(Y))
    theta = np.arctan2(Y, X)
    
    return R, theta

def intensity_by_distance_and_theta(img, x_center, y_center, **kwargs):
    """
    Return a polar intensity dataframe with calculated
    intensity at a range of theta, radial distance polar
    coordinates
    """
    distance_max = kwargs.get('distance_max', 20)
    distance_min = kwargs.get('distance_min', 0)
    theta_bin_size = kwargs.get('theta_bin_size', np.pi/16)
    dist_bin_size = kwargs.get('dist_bin_size', 1)
    thetas = np.arange(-np.pi, np.pi, theta_bin_size)
    dists = np.arange(distance_min, distance_max, dist_bin_size)
    # Initialize dataframe to be filled with intensities
    # as a function of (theta_bin, radial_distance_bin)
    df = pd.DataFrame(index=np.arange(0, len(thetas)*len(dists), 1))
    for i, t in enumerate(thetas):
        df.loc[i*len(dists):i*len(dists) + len(dists), 'theta'] = t
        for dist_i, dist in enumerate(dists):
            df.loc[i*len(dists) + dist_i, 'radial_distance'] = dist
    # Fill the dataframe above with calculated average intensity
    # for all pixels that fall within each theta, radial_distance bin
    R, theta = get_polar_coordinates(img, x_center, y_center)
    for index in df.index:
        t = df.loc[index, 'theta']
        dist = df.loc[index, 'radial_distance']
        # print(f'Looking for pixels at {t} radians, {dist} pixels')
        theta_mask = (np.greater(theta, t - theta_bin_size) & np.less(theta, t + theta_bin_size))
        dist_mask = (np.greater(R, dist - dist_bin_size) & np.less(R, dist + dist_bin_size))
        values = img[(theta_mask & dist_mask)]
        df.loc[index, 'intensity'] = np.mean(values)

    df.loc[:, 'x_center'] = x_center
    df.loc[:, 'y_center'] = y_center

    return df

def polar_intensity_peaks(polar_intensity_df, **kwargs):
    """
    For each theta value in polar_intensity_df, find peaks
    in intensity at varying radial_distance from center. 
    Add a column to the polar_intensity_df that records 
    distance at which largest peak was found, etc.

    polar_intensity_df should be created using
    byc.segmentation.intensity_by_distance_and_theta
    """
    # Number of pixels to add to distance at which
    # peak was found
    peak_dist_offset = kwargs.get('peak_offset', -4)
    peak_height_threshold = kwargs.get('peak_height_threshold', 100)
    medfilt_kern_size = kwargs.get('medfilt_kern_size', 7)
    df = polar_intensity_df
    df.loc[:, 'peak_X'] = np.nan
    df.loc[:, 'peak_Y'] = np.nan
    x_center = df.x_center.iloc[0]
    y_center = df.y_center.iloc[0]

    for t in df.theta.unique():
        intensity = df.loc[df.theta == t, 'intensity']
        intensity_medfilt = medfilt(intensity, kernel_size=medfilt_kern_size)
        peaks = find_peaks(intensity,
                           height=np.array([peak_height_threshold, None]))
        # Make sure peaks were found before trying to add them to the dataframe
        if len(peaks[0]) > 0:
            highest_peak_index = np.argmax(peaks[1]['peak_heights'])
            highest_peak_dist = peaks[0][highest_peak_index] + 1
            first_peak = peaks[0][0] + 1
            df.loc[df.theta==t, 'highest_peak_dist'] = highest_peak_dist + peak_dist_offset
            df.loc[df.theta==t, 'first_peak_dist'] = first_peak
            # Add cartesian coordinates for max amplitude peak found

            df.loc[df.theta==t, 'peak_X'] = highest_peak_dist*np.cos(t) + x_center
            df.loc[df.theta==t, 'peak_Y'] = highest_peak_dist*np.sin(t) + y_center
        else:
            print(f'No peaks found at {t} radians')

    return df
