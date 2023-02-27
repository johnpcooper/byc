import os
import re

import tifffile
import skimage.io
import skimage.filters
import skimage.transform
import skimage.util
import skimage.morphology
import numpy as np

from byc import rotation, utilities, constants, files, registration

class bycImageSet(object):
    """
    On instantation, prompt the user to select the input
    directory holding their raw tifs. Then create some
    information about that directory etc.
    """
    def __init__(self, **kwargs):
        
        self.input_dir_path = kwargs.get('input_dir_path', None)
        if self.input_dir_path == None:
            self.input_dir_path = files.select_directory('Choose the directory holding the expt you want to align')
        
        self.output_dir_path = self.output_dir_path(self.input_dir_path)
        self.fov_dir_paths_dict = self.fov_dir_paths_dict(self.input_dir_path)
        self.display_and_comments = files.get_byc_display_and_comments_dict(self.input_dir_path)
        self.channel_names = [channel['Name'] for channel in self.display_and_comments['Channels']]
        
    def output_dir_path(self, input_dir_path):
        """
        Create and return the path to an output dir
        in the source directory for input_dir_path
        """
        output_dir_path = os.path.join(os.path.dirname(input_dir_path), 'output')
        # Only make an output directory if none exists
        if not os.path.exists(output_dir_path):
            print(f'Making output dir at:\n{output_dir_path}')
            os.mkdir(output_dir_path)
            
        return output_dir_path
    
    def fov_dir_paths_dict(self, input_dir_path):
        """
        Find all the directories in input_dir_path
        that look like FOV directories and return a
        dictionary with keys = FOV index found in 
        filename referring to the path to that FOV
        """
        fov_dir_paths = []
        fov_indices = []
        for dirpath, dirnames, filenames in os.walk(input_dir_path):
            
            if 'Pos' in dirpath:
                fov_index = int(dirpath[dirpath.rindex('Pos') + 3:])
                fov_indices.append(fov_index)
                fov_dir_paths.append(dirpath)
    
        fov_dir_paths_dict = dict(zip(fov_indices, fov_dir_paths))
        return fov_dir_paths_dict
    
    def fov_channels_dict(self, fov_index, fov_dir_paths_dict, crop_frames=True, start_pixel_x=None, fraction=0.5):
        """
        Return fov_channels_dict: keys are channel
        names found in display_and_comments and elements
        are lists of images (one for each timepoint) 
        corresponding to that channel
        """
        fov_dir_path = fov_dir_paths_dict[fov_index]
        allfiles = os.listdir(fov_dir_path)
        tiffiles = [fn for fn in allfiles if '.tif' in fn]
        
        fov_channels_dict = {}
        for channel in self.channel_names:
            channel_paths = [os.path.join(fov_dir_path, fn) for fn in tiffiles if channel in fn]
            # Files don't get listed in order on unix based OS so 
            # we need to sort so that frames get read in in order
            channel_paths.sort()
            fov_channels_dict[channel] = [skimage.io.imread(path) for path in channel_paths]
        # Crop out data outside the middle half of each frame because
        # there are no cells there and we want as lean a dataset as possible
        if crop_frames == True:
            for key in fov_channels_dict.keys():
                stack = fov_channels_dict[key]
                frame0 = stack[0]
                total_width = frame0.shape[1]
                crop_width = int(np.round(total_width*fraction))
                
                if start_pixel_x == None:
                    start = int(np.round((total_width - crop_width)/2))
                else:
                    start = int(start_pixel_x)
                end = start + crop_width
                print(f'Cropping x axis to keep pixels {start} through {end} of total width {total_width}')
                cropped_stack = [frame[:, start:end] for frame in stack]
                fov_channels_dict[key] = cropped_stack

        
        return fov_channels_dict

def get_median_rotational_offset(bf_stack):
    """
    Choose three random images from the bf_stack,
    calculate their rotational offset from vertical,
    and return the median of those three offsets in
    degrees
    """

    image_rotations = []

    if len(bf_stack) > 3:
        rotation_sample_inds = np.random.choice(range(len(bf_stack)), 3, replace=False)
        rotation_sample_frames = [bf_stack[i] for i in rotation_sample_inds]
    else:
        rotation_sample_frames = bf_stack

    for index, image in enumerate(rotation_sample_frames):
        
        print(f'Determining rotation offset for sample frame {index+1} of {len(rotation_sample_frames)}')
        image_rotation = rotation.ImageRotation(image, compute=True)
        image_rotations.append(image_rotation)
        print(f'Found rotational offset of {str(image_rotation.offset)[0:6]} degrees')

    rotational_offsets_sample = [ir.offset for ir in image_rotations]
    
    median_offset = np.median(np.asarray(rotational_offsets_sample))
    return median_offset

def rotate_stack(stack, rotational_offset):
    """
    Take the stack and rotated it by
    rotational_offset degrees and return 
    the rotated stack as a list
    """
    rotated_images = [skimage.transform.rotate(img, rotational_offset) for img in stack]
    return rotated_images

def translate_stack(stack, translational_offsets):
    """
    Take the stack and translate each frame according
    to the offset with the same index in translational_offsets
    and return the stack of translated_images
    """
    assert len(stack) == len(translational_offsets), "Length of stack and offsets not the same!"
    
    translated_images = []
    for offset, image in zip(translational_offsets, stack):
        x = offset[1]
        y = offset[0]
        new_image = skimage.transform.warp(image, skimage.transform.AffineTransform(translation=(-x, -y)))
        translated_images.append(new_image)
        
    return translated_images

def rotate_channels(channels_dict, offset):
    """
    Return the dictionary of rotated
    channel stacks
    """
    rotated_channels = {}
    for channel_name, channel_stack in channels_dict.items():
        
        rotated_channel = rotate_stack(channel_stack, offset)
        rotated_channels[channel_name] = rotated_channel

    return rotated_channels

def translate_channels(channels_dict, offsets):
    """
    Return the dictionary of translated 
    channel stacks
    """
    translated_channels = {}
    for channel_name, channel_stack in channels_dict.items():

        translated_channel = translate_stack(channel_stack, offsets)
        translated_channels[channel_name] = translated_channel

    return translated_channels

def align_fov(
        fov_index,
        byc_image_set,
        write_output=True,
        rotate=False,
        **kwargs):
    """
    Pass this function an process.bycImageSet() instance which
    it will use to get data for the fov at fov_index, align
    each channel collected based on rotational and translational
    process of the 'Brightfield' channel, then if write_output
    save the translated channel_stacks in byc_image_set_output_dir_path.

    If write_output==False, return the translated channels_dict
    """
    save_unaligned = kwargs.get('save_unaligned', False)
    channels = byc_image_set.fov_channels_dict(fov_index,
                                               byc_image_set.fov_dir_paths_dict)
    if rotate:
        # Find a good rotational offset and rotate images in each channel
        # stack by that rotational offset
        median_offset = get_median_rotational_offset(channels['Brightfield'])
        rotated_channels = rotate_channels(channels, median_offset)
    else:
        rotated_channels = channels
    # Find translational offsets
    # Registration typicall works better when I 
    # use an image later in the stack for
    base_image_index = int(len(rotated_channels['Brightfield'])/2)
    # base_image_index = 75
    base_image = rotated_channels['Brightfield'][base_image_index]
    # Calculate the tranlsational registration offsets
    offsets = registration.determine_offsets(base_image, rotated_channels['Brightfield'])
    translated_channels = translate_channels(rotated_channels, offsets)
    # Save each translated channel stack
    fov_dir_path = byc_image_set.fov_dir_paths_dict[fov_index]
    # Extract the experiment date from the filepath
    date_pattern = constants.patterns.date
    match = re.search(date_pattern, fov_dir_path)
    if match:
        expt_date = match.group()
    else:
        print(f'No date found in FOV dir. path:\n{fov_dir_path}')
        expt_date = '00000000'
        print(f'Using {expt_date} as date')
    if write_output:

        fov_str = str(fov_index).zfill(2)
        basefilename = f'{expt_date}_byc_xy{fov_str}'
        base_writepath = os.path.join(byc_image_set.output_dir_path, basefilename)

        for channel_name, stack in translated_channels.items():
            print(f'Saving {channel_name} stack...')
            stack = skimage.util.img_as_uint(skimage.io.concatenate_images(stack))
            if channel_name == 'Brightfield':
                writepath = f'{base_writepath}_bf_stack.tif'
            else:
                writepath = f'{base_writepath}_{channel_name.lower()}_stack.tif'

            tifffile.imsave(writepath, stack)
            if save_unaligned:
                unaligned_stack = skimage.util.img_as_uint(skimage.io.concatenate_images(rotated_channels[channel_name]))
                tifffile.imsave(writepath.replace('stack.tif', 'unaligned_stack.tif'), unaligned_stack)

    else:
        return translated_channels


def align_byc_expt(**kwargs):
    """
    Find data using bycImageSet() and align all channels and
    fovs in that data, then save the data
    """
    input_path = kwargs.get('input_path', None)
    write_output = kwargs.get('write_output', True)
    save_unaligned = kwargs.get('save_unaligned', False)

    if input_path != None and os.path.exists(input_path):
        byc_image_set = bycImageSet(input_path)
    else:
        # Lack of input path arg to bycImageSet()
        # will prompt the user to select one
        byc_image_set = bycImageSet()

    for fov_index in byc_image_set.fov_dir_paths_dict.keys():
        print(f'Aligning FOV {fov_index+1} of {len(byc_image_set.fov_dir_paths_dict.keys())}')
        align_fov(fov_index, byc_image_set, write_output=write_output, save_unaligned=save_unaligned)
    print('Finished aligning')