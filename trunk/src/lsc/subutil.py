from pyraf import iraf
from astropy.io import fits
import os
import numpy as np
import scipy
import statsmodels.api as stats


def center_psf(psf):
    """Center psf at (0,0)"""

    psf = np.roll(psf, psf.shape[0] / 2, 0)
    psf = np.roll(psf, psf.shape[1] / 2, 1)
    return psf


def convert_to_background_fft(gain, background_fft, science_background_std, reference_background_std):
    """Convert params in fourier space to image space"""

    return background_fft * np.sqrt(science_background_std ** 2 + gain ** 2 * reference_background_std ** 2)


def fit_noise(data):
    """Find the standard deviation of the image background by fitting the background to a gaussian"""

    edge = 50
    trimmed_data = data[edge:-edge, edge:-edge]
    trimmed_data = trimmed_data[trimmed_data < np.percentile(trimmed_data, 90)]
    histogram_data = np.histogram(trimmed_data, bins=100)
    x = histogram_data[1][:-1]
    y = histogram_data[0]
    parameters, covariance = scipy.optimize.curve_fit(gauss, x, y, p0=[np.max(y), np.median(trimmed_data), np.std(trimmed_data)], maxfev=1600)
    return parameters[2]


def gauss(x, a, b, c):
    """Return a gaussian function"""

    return a * np.exp(-(x-b)**2/(2*c**2))


def read_psf_file(psf_filename):
    """Extract normalized psf array from iraf psf file"""

    iraf.noao()
    iraf.digiphot()
    iraf.daophot(_doprint=0)
    iraf.seepsf(psf_filename, 'temp.psf.fits')
    psf = fits.open('temp.psf.fits')[0].data
    psf /= np.sum(psf)
    os.system('rm temp.psf.fits')
    return psf


def resize_psf(psf, shape):
    """Resize centered (0,0) psf to given shape"""

    psf_extended = np.zeros(shape)
    center = [psf_extended.shape[0] / 2, psf_extended.shape[1] / 2]
    stamp = psf.shape
    vert_offsets = [center[0] - stamp[0] / 2, center[0] + stamp[0] / 2]
    horiz_offsets = [center[1] - stamp[1] / 2, center[1] + stamp[1] / 2]
    psf_extended[vert_offsets[0]: vert_offsets[1] + 1, horiz_offsets[0]: horiz_offsets[1] + 1] = psf
    return psf_extended


def solve_iteratively(science, reference):
    """Solve for linear fit iteratively"""

    gain_tolerance = 1e-8
    background_fft_tolerance = 1e-8
    gain = 1.
    background_fft = 0.
    gain0 = 10e5
    background_fft0 = 10e5
    i = 0
    max_iterations = 10

    # trim edge of image to reduce edge effects
    center, d = science.image_data.shape[0] / 2, science.image_data.shape[0] / 16
    coords = [center - d, center + d, center - d, center + d]

    science_image = science.image_data[coords[0]:coords[1], coords[2]:coords[3]]
    reference_image = reference.image_data[coords[0]:coords[1], coords[2]:coords[3]]

    science_image_fft = np.fft.fft2(science_image)
    reference_image_fft = np.fft.fft2(reference_image)
    science_psf_fft = np.fft.fft2(resize_psf(science.psf_data, science_image.shape))
    reference_psf_fft = np.fft.fft2(resize_psf(reference.psf_data, reference_image.shape))

    while abs(gain - gain0) > gain_tolerance or abs(background_fft - background_fft0) > background_fft_tolerance:

        denominator = science.background_std ** 2 * abs(reference_psf_fft) ** 2
        denominator += gain ** 2 * reference.background_std ** 2 * abs(science_psf_fft) ** 2
        science_convolved_image_fft = reference_psf_fft * science_image_fft / np.sqrt(denominator)
        reference_convolved_image_fft = science_psf_fft * reference_image_fft / np.sqrt(denominator)
        science_convolved_image = np.real(np.fft.ifft2(science_convolved_image_fft))
        reference_convolved_image = np.real(np.fft.ifft2(reference_convolved_image_fft))
        science_convolved_image_flatten = science_convolved_image.flatten()
        reference_convolved_image_flatten = reference_convolved_image.flatten()

        # remove bad pixels
        science_good_pix = remove_bad_pix(science_convolved_image_flatten, saturation=science.saturation_count)
        reference_good_pix = remove_bad_pix(reference_convolved_image_flatten, saturation=reference.saturation_count)
        good_pix_in_common = np.intersect1d(science_good_pix, reference_good_pix)
        science_convolved_image_flatten = science_convolved_image_flatten[good_pix_in_common]
        reference_convolved_image_flatten = reference_convolved_image_flatten[good_pix_in_common]

        gain0, background_fft0 = gain, background_fft

        x = stats.add_constant(reference_convolved_image_flatten)
        y = science_convolved_image_flatten
        robust_fit = stats.RLM(y, x).fit()
        parameters = robust_fit.params
        gain = parameters[1]
        background_fft = parameters[0]

        if i == max_iterations:
            break
        i += 1
        print('Iteration {}:'.format(i))
        background = convert_to_background_fft(gain, background_fft, science.background_std, reference.background_std)
        print('Beta = {0}, background = {1}'.format(gain, background))

    print('Fit done in {} iterations'.format(i))

    covariance = robust_fit.bcov_scaled[0, 0]
    background = convert_to_background_fft(gain, background_fft, science.background_std, reference.background_std)

    print('Beta = ' + str(gain))
    print('Gamma = ' + str(background))
    print('Beta Variance = ' + str(covariance))
    return gain, background


def remove_bad_pix(data, saturation=None, remove_background_pix=True, significance=3.):
    """Remove saturated and background pixels from dataset"""

    if saturation is not None:
        not_saturated_pix = np.where(data < saturation)
    else:
        not_saturated_pix = np.arange(data.size)

    if remove_background_pix:
        threshold = np.median(data) + significance * np.std(data)
        signal_pix = np.where(data > threshold)
    else:
        signal_pix = np.arange(data.size)

    good_pix_in_common = np.intersect1d(not_saturated_pix, signal_pix)
    return good_pix_in_common
