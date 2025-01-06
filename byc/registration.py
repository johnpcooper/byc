import os
import pandas as pd
import skimage.registration

def determine_offsets(base_image, rotated_images, save_offset_data=True):
    """
    Use skimage.registration.phase_cross_correlation() to 
    find the translational offset in (Y, X) pixels of each 
    image in raw_images from base_image. 
    
    Return the list of (Y, X) offset tuples
    """    
    offsets = []
    for index, image in enumerate(rotated_images):
        print(f'Determining registation offset for frame {index+1} of {len(rotated_images)}')
        pcc = skimage.registration.phase_cross_correlation(base_image,
                                                              image,
                                                              upsample_factor=20)
        # Just take the (Y, X) offset portion of the
        # phase_correlation output
        offset = pcc[0]
        print(f'Offset Y={offset[0]} and X={offset[1]}')
        offsets.append(offset)
    if save_offset_data:
        ys = [offset[0] for offset in offsets]
        xs = [offset[1] for offset in offsets]
        frames = range(len(offsets))
        offsets_df = pd.DataFrame(
            {
                'x_offset': xs,
                'y_offset': ys,
                'frame': frames
            }

        )
        offsets_df.to_csv('offsets.csv', index=False)
    return offsets

def filter_offsets(offsets, n=10):
    """
    Take an arry of (Y, X) offsets and replace outlier
    values with the average of the surrounding n offsets
    """
    pass