import os
import sys

import skimage
from skimage import io
import numpy as np

import matplotlib.pyplot as plt

from byc import files

def main(fraction=0.5, start_pixel_x=None):

    image_paths = files.select_files('choose images')
    writedir = os.path.join(os.path.dirname(image_paths[0]), 'cropped')
    if not os.path.exists(writedir):
        os.mkdir(writedir)
    print(f'Writing files to directory {writedir}')

    for path in image_paths:

        image = skimage.io.imread(path)
        if len(image.shape) == 3:        
            print(f'Opened image with dimension {image.shape} at {path}')

            total_width = image.shape[2]
            crop_width = int(np.round(total_width*fraction))
            
            if start_pixel_x == None:
                start = int(np.round((total_width - crop_width)/2))
            else:
                start = int(start_pixel_x)
            end = start + crop_width

            cropped_stack = image[:, :, start:end]
            writepath = path.replace(os.path.dirname(path), writedir)
            io.imsave(writepath, cropped_stack)

            print(f'Saved cropped image at {writepath}')

        elif len(image.shape) != 3:        
            print(f'Opened image with dimension {image.shape} at {path} ')
            print('Image array should have 3 dimensions')

if __name__ == '__main__':
    if sys.argv:
        print(sys.argv)
        if len(sys.argv) == 3:
            start_pixel_x = sys.argv[1]
            fraction = float(sys.argv[2])
        elif len(sys.argv) ==2:
            start_pixel_x = sys.argv[1]
            fraction = 0.5
            print(f'Defaulting to middle {fraction} of image starting at user input pixel {start_pixel_x}')
        else:
            start_pixel_x = None
            fraction= 0.5
            print(f'Defaulting to middle {fraction} of image starting 1/4 width into x axis')

        main(start_pixel_x=start_pixel_x, fraction=fraction)
    else:
        main()