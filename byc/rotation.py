import skimage.filters
import skimage.transform
import skimage.util
import skimage.morphology
import numpy as np
from scipy import ndimage

class ImageRotation(object):
    """
    Methods for taking 2D numpy array image, segmenting 
    vertical edges, detecting the angle of those vertical
    edges, and rotating the image by the offset of that
    ange from vertical
    """
    def __init__(self, image, **kwargs):
        self.image = image
        self.angle_window = 15
        self.rank_mean_disk_size = 9
        self.rank_mean_threshold = 200
        self.min_small_object_size = 500
        # Run vertical rotational offset detection
        compute = kwargs.get('compute', False)
        if compute:
            self.vertical_edges_img = self.vertical_edges(image)
            self.otsu_vertical_edges_img = self.otsu_vertical_edges(self.vertical_edges_img)
            self.rank_mean_vertical_edges_img = self.rank_mean_vertical_edges(self.otsu_vertical_edges_img)
            self.clean_vertical_edges_img = self.clean_vertical_edges(self.rank_mean_vertical_edges_img)
            self.offset = self.detect_rotational_offset(self.clean_vertical_edges_img)
            self.rotated_image = self.rotate_image(self.image, self.offset)


    def vertical_edges(self, image):
        """
        Take the image and return it sobel_v filtered
        (Find and threshold for vertical edges)
        """
        vertical_edges = skimage.filters.sobel_v(image)

        return vertical_edges
     
    def otsu_vertical_edges(self, vertical_edges_img):
        """
        Take an image that has ideally been sobel vertical
        filtered and threshold the image using an otsu algorithm
        detected threshold value
        """
        threshold = skimage.filters.threshold_otsu(vertical_edges_img)
        otsu_vertical_edges_img = vertical_edges_img <= threshold
        otsu_vertical_edges_img = skimage.util.img_as_ubyte(otsu_vertical_edges_img)

        return otsu_vertical_edges_img
    
    def rank_mean_vertical_edges(self, otsu_vertical_edges_img, **kwargs):
        """
        Take the otsu_vertical_edges image and rank mean filter
        with a default threshold mean_threshold=200 and default
        disk_size=9
        """
        disk_size = kwargs.get('disk_size', self.rank_mean_disk_size)
        mean_threshold = kwargs.get('mean_threshold', self.rank_mean_threshold)
        disk = skimage.morphology.disk(disk_size)

        rank_mean_vertical_edges_img = skimage.filters.rank.mean(otsu_vertical_edges_img, disk) < mean_threshold
        return rank_mean_vertical_edges_img

    def clean_vertical_edges(self, rank_mean_vertical_edges_img, **kwargs):
        """
        Take a rank mean smooothed vertical edge image,
        remove small objects with a default 
        min_object_size=500, run ndimage.binary_fill_holes
        on that cleaned image, and return the cleaned image
        """
        min_object_size = kwargs.get('min_object_size', self.min_small_object_size)
        clean_vertical_edges_img = skimage.morphology.remove_small_objects(rank_mean_vertical_edges_img, min_size=min_object_size)
        clean_vertical_edges_img = ndimage.binary_fill_holes(clean_vertical_edges_img)

        return clean_vertical_edges_img

    def detect_rotational_offset(self, clean_vertical_edges_img):
        """
        Take a cleaned, vertical edge segmented image,
        skeletonize it, then run a hough transform to
        find the offsets from vertical most commonly
        found in the vertical edge segmented image.
        Return the median offset found in hough transform
        in degrees
        """
        self.skeletons = skimage.morphology.skeletonize(clean_vertical_edges_img)
        angles_range = np.arange(-np.deg2rad(self.angle_window), np.deg2rad(self.angle_window), 0.0001)
        
        hough = skimage.transform.hough_line(self.skeletons, angles_range)
        # Hough line peaks identifies the most prominent lines separateed by 
        # a certain distance and angle in a Hough transform

        # compute the hough_line_peaks using hspace, angles, and dists that
        # were output from hough
        hspace = hough[0]
        angles = hough[1]
        dists = hough[2]
        args = [hspace, angles, dists]
        hough_peaks = skimage.transform.hough_line_peaks(*args)
        # Starting from the inside: 
        # transform.hough_line_peaks(*hough) returns the output of hough_line_peaks when passed the 
        # three variables that are output from transform.hough above
        # zip(*transform.hough_line_peaks(*hough)) zips the three lists that are output from hough_line_peaks.
        angles = [angle for h, angle, dist in zip(*skimage.transform.hough_line_peaks(*hough))]
        # Defining self.dists (for each angle found, what distance is that angle along x
        # axis of image?) and self.angles (array of the angle found for each vertical line 
        # found in the hough transform)
        self.angles = np.asarray(angles)
        self.dists = [dist for h, angle, dist in zip(*skimage.transform.hough_line_peaks(*hough))]

        offset = np.median(self.angles)
        return np.rad2deg(offset)

    def rotate_image(self, image, offset):
        """            
        Return an image (np.array()) rotated by the number of degrees
        returned by _determine_rotation_offset(image)
        """
        return skimage.transform.rotate(image, offset)


