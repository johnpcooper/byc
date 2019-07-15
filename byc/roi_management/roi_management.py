import tifffile as tf
import tkinter as tk
import tkinter.filedialog as dia
import skimage as ski
import skimage.io as io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skimage.filters import try_all_threshold
from skimage.filters import threshold_minimum
from skimage.filters import threshold_otsu
from skimage import img_as_uint

class Cell_Stacks(object):
    
    """ Doc strings """
    
    def __init__(self):
        
        # set the files that will be used for every cell analyzed
        self.master_cells_df = self.set_master_cells_df("Choose the .csv for the master index of this experiment")
        self.cells_dfs = self.set_cells_dfs(self.master_cells_df)

        self.image_stack_dsred = self.set_image_stack("Choose the dsred .tif stack ")
        self.image_stack_yfp = self.set_image_stack("Choose the yfp .tif stack ")
        self.image_stack_bf = self.set_image_stack("Choose the brightfield .tif stack ")
    
    def set_fp(self, prompt):
    
        """ Return the path to a file of interest. Call with the prompt you would like to 
            display in the dialog box."""

        # create the dialog box and set the fn
        root = tk.Tk()
        fp = dia.askopenfilename(parent=root, title=prompt)
        root.destroy() # very important to destroy the root object, otherwise python 
        # just keeps running in a dialog box

        return fp # return the path to the file you just selected
    
    def set_image_stack(self, prompt):

        """ Return a tf.imread() object referring to the image stack"""

        fp = self.set_fp(prompt)
        print("Path to the stack .tif selected:\n%s" % fp)
        image_stack = tf.imread(fp)
        # let the user see the tif stack they just instantiated
        #print("Displaying 0th image in this stack:")
        #io.imshow(image_stack[0])

        return image_stack
    
    def set_master_cells_df(self, prompt):

        """ Return a dataframe read from the master cells index .csv"""

        # define the path to the index .csv
        master_cells_fp = self.set_fp(prompt)
        # define the filename for the master expt index
        master_cells_df = pd.read_csv(master_cells_fp)

        return master_cells_df

    def set_cells_dfs(self, master_cells_df):

        """Return a list of sub-dataframes of the master_cells_df. The sub-dataframes
           contain data for only one cell according to its sub_coord"""

        # create the cells_dfs list
        cells_dfs = []

        # add to the cells_dfs a dataframe at index i for every
        # unique value in the master_cells_df['sub_coord']
        for i in master_cells_df['sub_coord'].unique():
            # set the logic mask to a sub-dataframe of master_cells_df containing
            # only values whose 'sub_coord' value is value
            print("value is ", i)
            logic_mask = (master_cells_df['sub_coord'] == i)
            cells_dfs.append(master_cells_df[logic_mask])
            print("cells_dfs is now %d elements long" % len(cells_dfs))

        return cells_dfs
    
    def set_cell_crop_params(self, image_stack, cells_dfs, cell_sub_coord):

        """ Return and save to the cwd a cropped image of the cell corresponding to
            cell_sub_coord by making a sub-matrix of the image_stack """

        # set the upper bounds as Y_ub and X_ub. 'Y' and 'X' are currently the lower bounds    
        cells_dfs[cell_sub_coord]['Y_ub'] = cells_dfs[cell_sub_coord]['Y'] + cells_dfs[cell_sub_coord]['Height']
        cells_dfs[cell_sub_coord]['X_ub'] = cells_dfs[cell_sub_coord]['X'] + cells_dfs[cell_sub_coord]['Width']
    
    def set_cropped_cell_stack(self, image_stack, cells_dfs, cell_sub_coord):
    
        """ Return a cropped stack made of the ROIs in cells_dfs """

        y_lb = cells_dfs[cell_sub_coord]['Y']
        y_ub = cells_dfs[cell_sub_coord]['Y_ub']
        x_lb = cells_dfs[cell_sub_coord]['X']
        x_ub = cells_dfs[cell_sub_coord]['X_ub']
        
        whole_fov_shape = self.image_stack_dsred[0].shape
        whole_fov_half_width = whole_fov_shape[1] / 2
        
        cropped_image_stack = []

        # iterate over ROI number
        for i in range(cells_dfs[cell_sub_coord]['Pos'].index.min(), cells_dfs[cell_sub_coord]['Pos'].index.max() + 1):

            print("Cropping slices with ROI %d of cell %d: " % (i, cell_sub_coord))
            

            if cells_dfs[cell_sub_coord]['Pos'].index.max() > i >= cells_dfs[cell_sub_coord]['Pos'].index.min():

                first_stack_index = cells_dfs[cell_sub_coord]['Pos'][i] - 1
                last_stack_index = cells_dfs[cell_sub_coord]['Pos'][i + 1] - 1

                # if true, then this is a single position ROI
                if (first_stack_index - last_stack_index) == 1:

                    stack_index = first_stack_index

                    print("Cropping stack index %d (single index ROI)" % stack_index)

                    source_image = image_stack[stack_index]
                    cropped_image = source_image[y_lb[i]: y_ub[i], x_lb[i]: x_ub[i]]
                    
                    # if the image is on the left side of the fov, then it should be flipped so that
                    # when it's resized extra pixels will be added on to the side that it's not dividing
                    # toward (ie pixels are always added to the right side for now, and this flipping loop
                    # makes it so that the left side is always the side pointing toward the trench)
                    if x_lb[i] < whole_fov_half_width:
                        cropped_image = np.flip(cropped_image)                        
                    else:
                        pass
                        
                    cropped_image_stack.append(cropped_image)                

                elif (first_stack_index - last_stack_index) != 1:                
                    # iterate over the slices of image_stack between the first two positions
                    for stack_index in range(first_stack_index, last_stack_index):
                        # this loop stops just short of last_stack_index, which will be the first index
                        # of the next ROI loop
                        print("Cropping stack index %d (multi-index ROI)" % stack_index)

                        source_image = image_stack[stack_index]
                        cropped_image = source_image[y_lb[i]: y_ub[i], x_lb[i]: x_ub[i]]
                        
                        if x_lb[i] < whole_fov_half_width:
                            cropped_image = np.flip(cropped_image)                        
                        else:
                            pass
                        
                        cropped_image_stack.append(cropped_image)

                else:
                    print("Error: Stack index %d may contain multiple ROIs" % first_stack_index)

            # add a cropped stack index at the last position of the image_stack. This may
            # potentially cause problems if the final stack index is part of a mult-index ROI
            elif i == cells_dfs[cell_sub_coord]['Pos'].index.max():
                first_stack_index = cells_dfs[cell_sub_coord]['Pos'][i] - 1

                stack_index = first_stack_index

                print("Cropping stack index %d (final stack index)" % stack_index)
                source_image = image_stack[stack_index]
                cropped_image = source_image[y_lb[i]: y_ub[i], x_lb[i]: x_ub[i]]
                
                if x_lb[i] < whole_fov_half_width:
                    cropped_image = np.flip(cropped_image)                        
                else:
                    pass
                
                cropped_image_stack.append(cropped_image)

            else: 
                pass
        
        return cropped_image_stack
    
    def set_slice_shapes_array(self, cropped_stack):

        """ Return a np.array() containing the shape of each slice in
            cropped_stack. For example, slice_shapes_array[:, 0] is a 
            list of the heights of each slice. """

        cropped_stack = cropped_stack.copy()
        slice_shapes_list = []

        for i in range(0, len(cropped_stack)):
            # Create a list of tuples, one tuple for each slice, where
            # the first element of the tuple is the height (# rows)
            # and the second element is width (# columns) of that slice

            slice_shapes_list.append(cropped_stack[i].shape)
        # convert the list of tuples into a 2d np.array() object
        slice_shapes_array = np.array(slice_shapes_list, dtype='uint16')

        return slice_shapes_array


    def set_max_dims(self, slice_shapes_array):

        """ Return a tuple which is (max height in slice_shapes_array, max width in slice_shapes_array)"""

        # maximum value in the heights column
        height_max = slice_shapes_array[:, 0].max()
        # maximum value in the widths column
        width_max = slice_shapes_array[:, 1].max()

        print("Largest slice in stack is %d rows by %d columns" % (height_max, width_max))

        return (height_max, width_max)
    
    def set_slice_offsets(self, cropped_stack, slice_shapes_array, height_max, width_max):
    
        """ Return a tuple which is (list of height_max - slice_shapes_array[:, 0],
            list of width_max - slice_shapes_array[:, 1]) """

        # define the differences in height and width between each slice and the maxima

        slice_height_offsets = []
        slice_width_offsets = []

        for i in range(0, len(cropped_stack)):

            slice_height_offsets.append(height_max - slice_shapes_array[i, 0])
            slice_width_offsets.append(width_max - slice_shapes_array[i, 1])

        return (slice_height_offsets, slice_width_offsets)
    
    def set_resized_stack(self, cropped_stack, slice_height_offsets, slice_width_offsets):
    
        """ Return the cropped_stack with added zeros to make each slice the same 
            size as the largest slice in the stack """

        resized_stack = []

        for i in range(0, len(cropped_stack)):

            image = cropped_stack[i]
            resized_image = image.copy()
            image_h_offset = slice_height_offsets[i]
            image_w_offset = slice_width_offsets[i]

            h_filler_value = resized_image[-1, :].min()
            h_filler = np.full((image_h_offset, resized_image.shape[1]), h_filler_value, dtype='uint16')
            resized_image = np.append(resized_image, h_filler, axis=0)

            w_filler_value = resized_image[:, -1].min()
            w_filler = np.full((resized_image.shape[0], image_w_offset), w_filler_value, dtype='uint16')
            resized_image = np.append(resized_image, w_filler, axis=1)

            resized_stack.append(resized_image)

        return resized_stack
    
    def set_min_thresholded_stack(self, resized_stack, plot_results):
    
        """ Return a minimum thresholded version of resized_stack. Note that the returned 
            image is a matrix of booleans, not integers. The stack can be converted to 
            integers using skimage.io.img_as_uint. plot_results is a boolean for whether to 
            plot the results of thresholding for quality control. I typically set this to False.
            Cell area in the returned stack is True, everywhere else is False. """

        thresholded_stack = []

        for i in range(0, len(resized_stack)):

            print("Applying minimum threshold to image %s of %s" % (i, len(resized_stack)-1))

            image = resized_stack[i]
            thresh = threshold_minimum(image)
            binary = image > thresh

            # add fresh binary image to stack
            thresholded_stack.append(binary)

            if plot_results == True:    

                fig, axes = plt.subplots(ncols=2, figsize=(8, 3))
                ax = axes.ravel()

                ax[0].imshow(image, cmap=plt.cm.gray)
                ax[0].set_title('Original image')

                ax[1].imshow(binary, cmap=plt.cm.gray)
                ax[1].set_title('Minimum Threshold')

                for a in ax:
                    a.axis('off')

                plt.show()

            else: 
                pass

        return thresholded_stack
    
    def set_otsu_thresholded_stack(self, resized_stack, plot_results):
    
        """ Return a otsu thresholded version of resized_stack. Note that the returned 
            image is a matrix of booleans, not integers. The stack can be converted to 
            integers using skimage.io.img_as_uint. plot_results is a boolean for whether to 
            plot the results of thresholding for quality control. I typically set this to False.
            Cell area in the returned stack is True, everywhere else is False. """

        thresholded_stack = []

        for i in range(0, len(resized_stack)):

            print("Applying otsu threshold to image %s of %s" % (i, len(resized_stack)-1))

            image = resized_stack[i]
            thresh = threshold_otsu(image)
            binary = image > thresh

            # add fresh binary image to stack
            thresholded_stack.append(binary)

            if plot_results == True:    

                fig, axes = plt.subplots(ncols=2, figsize=(8, 3))
                ax = axes.ravel()

                ax[0].imshow(image, cmap=plt.cm.gray)
                ax[0].set_title('Original image')

                ax[1].imshow(binary, cmap=plt.cm.gray)
                ax[1].set_title('Otsu Threshold')

                for a in ax:
                    a.axis('off')

                plt.show()

            else: 
                pass

        return thresholded_stack

    def set_cropped_stacks(self, cell_index, cells_dfs):

        """ Return a stack cropped from the source dsred, yfp, and brightfield stacks
            according to the cell index parameters found in master index .csv (this gets 
            read into a pd.DataFrame at instantiation of Cell_Stacks() """  

        self.set_cell_crop_params(self.image_stack_dsred, cells_dfs, cell_index)
        self.cropped_stack_dsred = self.set_cropped_cell_stack(self.image_stack_dsred, cells_dfs, cell_index)
        self.cropped_stack_yfp = self.set_cropped_cell_stack(self.image_stack_yfp, cells_dfs, cell_index)
        self.cropped_stack_bf = self.set_cropped_cell_stack(self.image_stack_bf, cells_dfs, cell_index)

        return self.cropped_stack_dsred, self.cropped_stack_yfp, self.cropped_stack_bf


    def set_cropped_cells_lists(self, cells_dfs):
        
        """ Return three lists - dsred_stacks, yfp_stacks_, and bf_stacks,
            with each element being the cropped stack of that cell """
        
        cropped_dsred_stacks = []
        cropped_yfp_stacks = []
        cropped_bf_stacks = []
        
        # remember that cells_dfs is a list of dfs, one for each
        # unique cell in the master_cells_dfs
        for cell_index in range(0, len(cells_dfs)):
            
            print("Setting cropped stacks for each channel (dsred, yfp, bf) of cell %s" % cell_index)
            
            cropped_dsred_stack, cropped_yfp_stack, cropped_bf_stack = self.set_cropped_stacks(cell_index, cells_dfs)
            
            cropped_dsred_stacks.append(cropped_dsred_stack)
            cropped_yfp_stacks.append(cropped_yfp_stack)
            cropped_bf_stacks.append(cropped_bf_stack)
            
        return cropped_dsred_stacks, cropped_yfp_stacks, cropped_bf_stacks
            
    def set_resized_cells_lists(self, cropped_dsred_stacks, cropped_yfp_stacks, cropped_bf_stacks, cells_dfs):
        
        """ Return three lists - dsred_stacks, yfp_stacks_, and bf_stacks,
            wit each element being the resized stack of that cell """
        
        resized_dsred_stacks = []
        resized_yfp_stacks = []
        resized_bf_stacks = []
        
        # remember that cells_dfs is a list of dfs, one for each
        # unique cell in the master_cells_dfs
        for cell_index in range(0, len(cells_dfs)):
            
            print("Setting resized stacks for each channel (dsred, yfp, bf) of cell %s" % cell_index)
            
            # use the cropped_dsred_stack of current cell to define slice shapes
            slice_shapes_array = self.set_slice_shapes_array(cropped_dsred_stacks[cell_index])
            height_max, width_max = self.set_max_dims(slice_shapes_array)
            
            # the slice shapes offsets
            slice_height_offsets, slice_width_offsets = self.set_slice_offsets(cropped_dsred_stacks[cell_index],
                                                                             slice_shapes_array,
                                                                             height_max,
                                                                             width_max)
            
            resized_dsred_stack = self.set_resized_stack(cropped_dsred_stacks[cell_index], slice_height_offsets, slice_width_offsets)
            resized_yfp_stack = self.set_resized_stack(cropped_yfp_stacks[cell_index], slice_height_offsets, slice_width_offsets)
            resized_bf_stack = self.set_resized_stack(cropped_bf_stacks[cell_index], slice_height_offsets, slice_width_offsets)
            
            resized_dsred_stack = io.concatenate_images(resized_dsred_stack)
            resized_yfp_stack = io.concatenate_images(resized_yfp_stack)
            resized_bf_stack = io.concatenate_images(resized_bf_stack)
            
            resized_dsred_stacks.append(resized_dsred_stack)
            resized_yfp_stacks.append(resized_yfp_stack)
            resized_bf_stacks.append(resized_bf_stack)
            
        return resized_dsred_stacks, resized_yfp_stacks, resized_bf_stacks

    def set_min_thresholded_cells_lists(self, resized_dsred_stacks, cells_dfs):
        
        """ Return a list of min thresholded stacks made from the resized_stack 
            list passed to this function """
        
        min_thresholded_dsred_stacks = []
        
        # remember that cells_dfs is a list of dfs, one for each
        # unique cell in the master_cells_dfs
        for cell_index in range(0, len(cells_dfs)):
            print("Determining min threshold for resized dsred stack of cell %s" %  cell_index)
            
            resized_stack = resized_dsred_stacks[cell_index]
            thresholded_stack = self.set_min_thresholded_stack(resized_stack, plot_results=False)
            final_thresholded_stack = io.concatenate_images(img_as_uint(thresholded_stack))
            
            min_thresholded_dsred_stacks.append(final_thresholded_stack)
            
            #now save the stacks
            cell_df = cells_dfs[cell_index]
            
            path = str(cell_df['path'][cell_df.index.min()] + "\\")
            expt_date = cell_df['date'][cell_df.index.min()]
            expt_type = cell_df['expt_type'][cell_df.index.min()]
            xy = cell_df['xy'][cell_df.index.min()]
            sub_coord = cell_df['sub_coord'][cell_df.index.min()]
            
            stack_title = str("%s_%s_xy%s_cell%s" % (expt_date, expt_type, xy, sub_coord))
            io.imsave(path + stack_title + "dsred_min_thresh.tif", final_thresholded_stack)
            
        return min_thresholded_dsred_stacks

    def set_otsu_thresholded_cells_lists(self, resized_dsred_stacks, cells_dfs):
        
        """ Return a list of otsu thresholded stacks made from the resized_stack 
            list passed to this function """
        
        otsu_thresholded_dsred_stacks = []
        
        # remember that cells_dfs is a list of dfs, one for each
        # unique cell in the master_cells_dfs
        for cell_index in range(0, len(cells_dfs)):
            print("Determining otsu threshold for resized dsred stack of cell %s" %  cell_index)
            
            resized_stack = resized_dsred_stacks[cell_index]
            thresholded_stack = self.set_otsu_thresholded_stack(resized_stack, plot_results=False)
            final_thresholded_stack = io.concatenate_images(img_as_uint(thresholded_stack))
            
            otsu_thresholded_dsred_stacks.append(final_thresholded_stack)
            
            #now save the stacks
            cell_df = cells_dfs[cell_index]
            
            path = str(cell_df['path'][cell_df.index.min()] + "\\")
            expt_date = cell_df['date'][cell_df.index.min()]
            expt_type = cell_df['expt_type'][cell_df.index.min()]
            xy = cell_df['xy'][cell_df.index.min()]
            sub_coord = cell_df['sub_coord'][cell_df.index.min()]
            
            stack_title = str("%s_%s_xy%s_cell%s" % (expt_date, expt_type, xy, sub_coord))
            io.imsave(path + stack_title + "dsred_otsu_thresh.tif", final_thresholded_stack)
            
        return otsu_thresholded_dsred_stacks
        

    def save_resized_stacks(self, resized_dsred_stacks, resized_yfp_stacks, resized_bf_stacks, cells_dfs):
        
        """ Return nothing, save all cropped and resized stacks in the directory
            specified in the master cells index """
        
        for cell_index in range(0, len(cells_dfs)):
            
            dsred_stack = resized_dsred_stacks[cell_index]
            yfp_stack = resized_yfp_stacks[cell_index]
            bf_stack  = resized_bf_stacks[cell_index]
            
            cell_df = cells_dfs[cell_index]
            
            path = str(cell_df['path'][cell_df.index.min()] + "\\")
            expt_date = cell_df['date'][cell_df.index.min()]
            expt_type = cell_df['expt_type'][cell_df.index.min()]
            xy = cell_df['xy'][cell_df.index.min()]
            sub_coord = cell_df['sub_coord'][cell_df.index.min()]
        
            
            stack_title = str("%s_%s_xy%s_cell%s" % (expt_date, expt_type, xy, sub_coord))
            
            print("Saving resized stacks for cell %s in %s" % (cell_index, path))
            io.imsave(path + stack_title + "dsred.tif", dsred_stack)
            io.imsave(path + stack_title + "yfp.tif", yfp_stack)
            io.imsave(path + stack_title + "bf_stack.tif", bf_stack)
            
    def save_min_thresh_stacks(self, resized_dsred_stacks, resized_yfp_stacks, resized_bf_stacks):
        
        """ Return nothing, save all cropped, resized, and thresholded stacks in the directory 
            specified in the master cells index """
    
    