import tifffile as tf
import tkinter as tk
import tkinter.filedialog as dia
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from skimage.io import imsave, concatenate_images
from skimage.filters import threshold_otsu
from skimage import img_as_uint
from read_roi import read_roi_file, read_roi_zip

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
    
    def set_cell_crop_roi_dfs(self, master_cells_df):
        """
        Return a list of DataFrames, one for each cell. The coordinates in each
        of these DataFrames will be used to crop from the image stacks in 
        set_cropped_cell_stack_list()
        """
    
        cell_crop_roi_dfs = []

        for cell_index in master_cells_df.index:

            expt_path = master_cells_df.path[cell_index]
            expt_date = int(master_cells_df.date[cell_index])
            expt_type = master_cells_df.expt_type[cell_index]
            cell_rois_fp = f"{expt_path}\\{expt_date}_{expt_type}_cell{str(cell_index).zfill(2)}_crop_rois.zip"

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
        #cell_rois_df.to_csv(f"{expt_path}\\{expt_date}_{expt_type}_cell{str(cell_index).zfill(2)}_crop_rois.csv")
        return cell_crop_roi_dfs
    
    def set_cell_channel_stacks(self, master_cells_df, cell_index):
        """
        Return a list of tf.imread() objects (one for each channel collected in expt)
        read according to data in the master_cells_df
        """ 

        expt_path = master_cells_df.path[cell_index]
        expt_date = int(master_cells_df.date[cell_index])
        expt_type = master_cells_df.expt_type[cell_index]
        xy = str(int(master_cells_df.xy[cell_index]))
        # Should be the channel names used to create stacks as output
        # of alignment. Channel names are separated by spaces in master_cells_df
        self.channel_names = master_cells_df.channels_collected[cell_index].split()

        # Set a list of paths to the channel stacks
        channel_stacks_fps = [f"{expt_path}\\{expt_date}_{expt_type}_xy{xy.zfill(2)}_{channel_name}_stack.tif" for channel_name in self.channel_names]
        
        print(f"Paths to channel stacks for cell {cell_index} found based on master df")
        # Let the user know which stacks were read
        for fp in channel_stacks_fps:
            print(fp)
            
        # Set a list of the actual channel stack imread() objects
        cell_channel_stacks = [tf.imread(filepath) for filepath in channel_stacks_fps]
        
        return cell_channel_stacks
    
    def set_cell_cropped_stacks_dict(self, master_cells_df, cell_rois_dfs, cell_index):

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
            print(f"max height {max_height} - frame height {frame_height} = offset {height_offset}")
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
                print(f"Shape of resized frame after adding h_filler_array {resized_frame.shape} and height offest {h_offset}")

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
    
    for cell_index in cs.master_cells_df.sub_coord:    
    
        # Run all the cropping and processing method of Cell_Stack on the cell
        # with cell_index
        cell_rois_dfs = cs.set_cell_crop_roi_dfs(cs.master_cells_df)
        cell_cropped_channels_dict = cs.set_cell_cropped_stacks_dict(cs.master_cells_df, cell_rois_dfs, cell_index)
        cs.add_cell_otsu_thresholded_stack()
        resized_channels_dict = cs.set_resized_cell_cropped_channels_dict(cs.cell_cropped_channels_dict)

        # Save each stack 
        expt_path = cs.master_cells_df.path[cell_index]
        expt_date = int(cs.master_cells_df.date[cell_index])
        expt_type = cs.master_cells_df.expt_type[cell_index]
        xy = str(int(cs.master_cells_df.xy[cell_index]))

        for channel_name, stack in resized_channels_dict.items():
            filename = f'{expt_date}_{expt_type}_xy{xy}_cell{cell_index}_{channel_name}_stack.tif'
            save_path = f'{expt_path}//{filename}'
            try:

                imsave(save_path, concatenate_images(stack))
            except:
                print(f"Could not save stacks for cell {cell_index}, img dims must agree")