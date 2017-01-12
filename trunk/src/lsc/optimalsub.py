import subutil
from astropy.io import fits
import numpy as np


class ImageClass:
    """Contains the image and relevant parameters"""

    def __init__(self, image_filename, psf_filename):
        self.image_filename = image_filename
        self.psf_filename = psf_filename
        self.image_data = fits.getdata(image_filename)
        self.psf_data = subutil.read_psf_file(psf_filename)
        self.zero_point = 1.
        self.background_counts = 0.
        self.background_std = subutil.fit_noise(self.image_data)
        self.saturation_count = 10e6


def calculate_difference_image(science, reference, normalization='reference', output='output.fits'):
    """Calculate the difference image using the Zackey algorithm"""

    # match the gains
    science.zero_point, science.background_counts = subutil.solve_iteratively(science, reference)
    reference.zero_point = 1.
    zero_point_ratio = science.zero_point / reference.zero_point

    # create required arrays
    science_image = science.image_data - science.background_counts
    reference_image = reference.image_data - reference.background_counts
    science_psf = subutil.center_psf(subutil.resize_psf(science.psf_data, science_image.shape))
    reference_psf = subutil.center_psf(subutil.resize_psf(reference.psf_data, reference_image.shape))

    # do fourier transforms (fft)
    science_image_fft = np.fft.fft2(science_image)
    reference_image_fft = np.fft.fft2(reference_image)
    science_psf_fft = np.fft.fft2(science_psf)
    reference_psf_fft = np.fft.fft2(reference_psf)

    # calculate difference image
    denominator = science.background_std ** 2 * zero_point_ratio ** 2 * abs(science_psf_fft) ** 2
    denominator += reference.background_std ** 2 * abs(reference_psf_fft) ** 2
    difference_image_fft = reference_psf_fft * science_image_fft
    difference_image_fft -= zero_point_ratio * science_psf_fft * reference_image_fft
    difference_image_fft /= np.sqrt(denominator)
    difference_image = np.fft.ifft2(difference_image_fft)
    difference_image = normalize_difference_image(difference_image, science, reference, normalization=normalization)
    
    # save difference image to file
    hdu = fits.PrimaryHDU(np.real(difference_image))
    hdu.header = fits.getheader(science.image_filename)
    hdu.header['PHOTNORM'] = normalization
    hdu.header['CONVOL00'] = normalization
    hdu.writeto(output, clobber=True)

    return difference_image


def calculate_difference_image_zero_point(science, reference):
    """calculate the flux based zero point of the difference image"""

    zero_point_ratio = science.zero_point / reference.zero_point
    denominator = science.background_std ** 2 + reference.background_std ** 2 * zero_point_ratio ** 2
    difference_image_zero_point = zero_point_ratio / np.sqrt(denominator)
    return difference_image_zero_point


def normalize_difference_image(difference, science, reference, normalization='reference'):
    """normalize to user's choice of image"""

    difference_image_zero_point = calculate_difference_image_zero_point(science, reference)
    if normalization == 'reference' or normalization == 't':
        difference_image = difference * reference.zero_point / difference_image_zero_point
    elif normalization == 'science' or normalization == 'i':
        difference_image = difference * science.zero_point / difference_image_zero_point
    else:
        difference_image = difference
    return difference_image
