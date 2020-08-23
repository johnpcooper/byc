import skimage.registration

def determine_offsets(base_image, rotated_images):
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
        
    return offsets


