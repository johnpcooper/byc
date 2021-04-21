import os

import skimage
from skimage import io
import numpy as np
from scipy import ndimage

import matplotlib.pyplot as plt

from byc import utilities, files, process, constants, rotation

def main():

    fraction=0.5
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

            start = int(np.round((total_width - crop_width)/2))
            end = start + crop_width

            cropped_stack = image[:, :, start:end]
            writepath = path.replace(os.path.dirname(path), writedir)
            io.imsave(writepath, cropped_stack)

            print(f'Saved cropped image at {writepath}')

        elif len(image.shape) != 3:        
            print(f'Opened image with dimension {image.shape} at {path} ')
            print('Image array should have 3 dimensions')

if __name__ == '__main__':
    main()