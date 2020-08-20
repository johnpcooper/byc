from skimage.filters import rank, threshold_otsu, sobel_v
from skimage.morphology import disk, remove_small_objects
from skimage import transform
from skimage.morphology import skeletonize
from skimage.feature import register_translation
from skimage.registration import phase_cross_correlation
from skimage.transform import AffineTransform
import skimage.io as io
from skimage.util import img_as_uint
import random

from scipy import ndimage
import numpy as np

import time
import math
import logging
import os
import logging
import tifffile as tf
import tkinter as tk
import tkinter.filedialog as dia


log = logging.getLogger(__name__)

class Constants(object):
    FIFTEEN_DEGREES_IN_RADIANS = 0.262
    ACCEPTABLE_SKEW_THRESHOLD = 5.0
    NUM_CATCH_CHANNELS = 28
    
def create_vertical_segments(image_data):
    """
    Creates a binary image with blobs surrounding areas that have a lot of vertical edges
    """
    # find edges that have a strong vertical direction
    vertical_edges = sobel_v(image_data)
    # Sepearate out the areas where there is a large amount of vertically
    # oriented stuff
    return _segment_edge_areas(vertical_edges)
    
def _segment_edge_areas(edges, disk_size=9, mean_threshold=200, min_object_size=500):
    """
    
    Takes a greyscale image (with brighter colors corresponding to edges) and returns
    a binary image with high edge density and black indicates low density
    
    param image_data: a 2D numpy array
    
    """
    
    # convert the greyscale edge information into black and white image
    threshold = threshold_otsu(edges)
    # Filter out the edge data below the threshold, effectively removing some noise
    raw_channel_areas = edges <= threshold
    return ndimage.binary_fill_holes(raw_channel_areas)


# from fylm/service/rotation.py
def _determine_rotation_offset(image):
    """
    Finds rotational skew so that the sides of the central trench are (nearly) perfectly vertical.
    :param image:   raw image data in a 2D (i.e. grayscale) numpy array
    :type image:    np.array()
    """
    segmentation = create_vertical_segments(image)
    # Draw a line that follows the center of the segments at each point, which should be roughly vertical
    # We should expect this to give us four approximately-vertical lines, possibly with many gaps in each line
    skeletons = skeletonize(segmentation)
    # Use the Hough transform to get the closest lines that approximate those four lines
    hough = transform.hough_line(skeletons, np.arange(-Constants.FIFTEEN_DEGREES_IN_RADIANS,
                                                      Constants.FIFTEEN_DEGREES_IN_RADIANS,
                                                      0.0001))
    # Create a list of the angles (in radians) of all of the lines the Hough transform produced, with 0.0 being
    # completely vertical
    # These angles correspond to the angles of the four sides of the channels, which we need to correct for
    angles = [angle for _, angle, dist in zip(*transform.hough_line_peaks(*hough))]
    if not angles:
        log.warn("Image skew could not be calculated. The image is probably invalid.")
        return 0.0
    else:
        # Get the average angle and convert it to degrees
        offset = sum(angles) / len(angles) * 180.0 / math.pi
        if offset > Constants.ACCEPTABLE_SKEW_THRESHOLD:
            log.warn("Image is heavily skewed. Check that the images are valid.")
        return offset
    
def rotate_image(image, offset):
    """
    
    Return an image (np.array()) rotated by the number of degrees
    returned by _determine_rotation_offset(image)
    
    """
    
    return transform.rotate(image, offset)

# from fylm/service/registration.py
def _determine_registration_offset(base_image, uncorrected_image):
    """
    
    Finds the translational offset required to align this image with all others in the stack.
    Returns dx, dy adjustments in pixels *but does not change the image!*
    
    :param base_image:   a 2D numpy array that the other image should be aligned to
    :param uncorrected_image:   a 2D numpy array
    :returns:   float, float
    
    """

    # Get the dimensions of the images that we're aligning
    base_height, base_width = base_image.shape
    uncorrected_height, uncorrected_width = uncorrected_image.shape

    # We take the area that roughly corresponds to the catch channels. This has two benefits: one, it
    # speeds up the registration significantly (as it scales linearly with image size), and two, if
    # a large amount of debris/yeast/bacteria/whatever shows up in the central trench, the registration
    # algorithm goes bonkers if it's considering that portion of the image.
    # Thus we separately find the registration for the left side and right side, and average them.
    left_base_section = base_image[:, int(base_width * 0.1): int(base_width * 0.3)]
    left_uncorrected = uncorrected_image[:, int(uncorrected_width * 0.1): int(uncorrected_width * 0.3)]
    right_base_section = base_image[:, int(base_width * 0.7): int(base_width * 0.9)]
    right_uncorrected = uncorrected_image[:, int(uncorrected_width * 0.7): int(uncorrected_width * 0.9)]

    # 
    left_dy, left_dx = phase_cross_correlation(left_base_section, left_uncorrected, upsample_factor=20)[0]
    right_dy, right_dx = phase_cross_correlation(right_base_section, right_uncorrected, upsample_factor=20)[0]

    return (left_dy + right_dy) / 2.0, (left_dx + right_dx) / 2.0

def translate_image(uncorrected_image, translational_offset):
    x = translational_offset[1]
    y = translational_offset[0]
    new_image = transform.warp(uncorrected_image, transform.AffineTransform(translation=(-x, -y)))
    return new_image

def get_channel_names(fov_path):

    fov_slices_filenames = os.listdir(fov_path)
    rep_image_filename = fov_slices_filenames[0]
    rep_image = tf.imread(fov_path + '/%s' % rep_image_filename)
    len_channels = len(rep_image)
    # if it's a single channel image, len_channels will return number of pixels, not 
    # the actual number of channels so I'm filtering for that here
    if len_channels > 10:
        len_channels = 1
    else: 
        pass

    print("Detected {} channels".format(len_channels))
    print("Enter names below:")

    channel_names = []
    for i in range(0, len_channels):
        channel_name = input("Name for channel {}: ".format(i))
        channel_names.append(channel_name)

    return len_channels, channel_names

def align_images(fov_path, channel_names):

    fov_slice_filenames = [file for file in os.listdir(fov_path) if file.endswith('.tif')]

    images = []
    for filename in fov_slice_filenames:
        image = tf.imread(fov_path + '/%s' % filename)
        images.append(image)
        print(filename)

    # ascertain the number of channels in the image, assume that 
    # channel 0 is brightfield

    # Calculate rotational offset for three random frames
    # in FOV stack
    if len(images) >= 3:
        frames_to_align = np.random.choice(range(len(images)), 3, replace=False)
    elif len(image) < 3:
        frames_to_align = [0]

    rotational_offsets = []
    for frame in frames_to_align:
        image = images[frame]    
        print(f"Determining rotation offset for frame {frame}")
        if image.ndim < 3: # if image only has one channel
        # assume that that one channel is the brightfield channel
        # and align based on that channel
            vis_channel = image
        elif image.ndim >= 3: # if image has more than one channel,
        # align based on the first channel which should be brightfied (bf)
            vis_channel = image[0]
        rotational_offset = _determine_rotation_offset(vis_channel)
        print(rotational_offset)
        rotational_offsets.append(rotational_offset)

    rotational_offset = np.median(np.array(rotational_offsets))
    final_rt_offests_arr = np.full(shape=(len(images)), fill_value=(rotational_offset))

    # create a list of rotationally aligned images rotated according to the final_rt_offsets_arr array
    rotated_images = []
    index = 0 # again, progress bar index
    for i, frame in enumerate(images):
        print("Rotating image %d of %d" % (i, len(images)))

        if frame.ndim < 3:
            rotated_channels_i = rotate_image(frame, final_rt_offests_arr[i])
        elif frame.ndim >= 3:
            channels_i = []
            for j in range(0, len(frame)):
                channels_i.append(frame[j])

            index += 1

            rotated_channels_i = []
            for j in range(0, len(channels_i)):
                rotated_channels_i.append(rotate_image(channels_i[j], final_rt_offests_arr[i]))

        rotated_images.append(rotated_channels_i)
    
    # calculate translational offsets
    index = 0
    translational_offsets = []
    for i in range(0, len(rotated_images)):
        image = rotated_images[i][0]
        index += 1
        print("Determining registration offset %d of %d" % (index, len(images)))
        try:
            # align based on feautres of image 0, I think I should change this to align
            # to image i-1
            translational_offset = _determine_registration_offset(rotated_images[0][1], image)
            print(f'Translational offset: {translational_offset}')
            translational_offsets.append(translational_offset)
        except:
            translational_offset = _determine_registration_offset(rotated_images[0], rotated_images[i])
            print(f'Translational offset: {translational_offset}')
            translational_offsets.append(translational_offset)
        
    # now translate the images
    translated_images = []
    index = 0
    for i in range(0, len(translational_offsets)):
        index += 1
        print("Translating image %d of %d" % (index, len(images)))
        
        translated_channels_i = []

        try:
            for j in range(0, len(images[0])):
                translated_channels_i.append(translate_image(rotated_images[i][j], translational_offsets[i]))
        except:
            translated_channels_i = (translate_image(rotated_images[i], translational_offsets[i]))

        translated_images.append(translated_channels_i)
    
    # make a dictionary to hold images for each channel
    keys = channel_names
    values = []

    if len(channel_names) > 1:
        for j in range(0, len(images[0])):        
            values.append([translated_images[i][j] for i in range(0, len(translated_images))])
    else:
        values.append([translated_images[i] for i in range(0, len(translated_images))])

    translated_images_dict = dict(zip(keys, values))

    return translated_images_dict

def save_stacks(translated_images_dict, save_path, fov_name):
    print("Saving stacks...")
    for channel, stack in translated_images_dict.items():
        print("Saving {} stack".format(channel))
        concat_stack = io.concatenate_images(img_as_uint(stack))
        tf.imsave(save_path + '/{}_{}_stack.tif'.format(fov_name, channel), concat_stack)