import numpy as np
import scipy
from astropy.io import fits
from pyraf import iraf
import time
from os import system
import os
import argparse
from lsc.lscastrodef import sextractor
import statsmodels.api as sm
import matplotlib.pyplot as plt
import logging


def AlignImages(NewFile, RefFile):
    """Align reference image to new image when both images have WCS"""

    TempFile = NewFile.replace('.fits', '.ref.fits')
    iraf.wregister(RefFile, NewFile, TempFile)
    #system('cp {0} {1}'.format(TempFile, NewFile.replace('.fits', '.ref.fits')))
    AlignedRef = fits.getdata(TempFile)
    return AlignedRef


def CleanDataIndex(Data, SatCount = None, RmBackPix = True, Significance = 1.):
    """Remove saturated and background pixels from dataset"""

    if SatCount is not None:
        UnSatPixIndex = np.where(Data < SatCount)
    else:
        UnSatPixIndex = np.arange(Data.size)

    if RmBackPix:
        Thresh = np.median(Data) +  Significance * np.std(Data)
        NotBackPixIndex = np.where(Data > Thresh)
    else:
        NotBackPixIndex = np.arange(Data.size)
    
    GoodPixInCommon = np.intersect1d(UnSatPixIndex, NotBackPixIndex)
    return GoodPixInCommon


def CropImages(N, R, Pn, Pr):
    """Open and crop images and PSFs to the same size"""

    sN = N.shape
    #sR = R.shape
    #s = [np.min([sN[0], sR[0]]), np.min([sN[1], sR[1]])]
    #N = N[:s[0], :s[1]]
    #R = R[:s[0], :s[1]]

    Pn_extended = PadPSF(Pn, sN)
    Pr_extended = PadPSF(Pr, sN)

    return N, R, Normalize(Pn_extended), Normalize(Pr_extended)


def ExtractPSF(Pf):
    """Extract PSF array from iraf PSF file"""

    iraf.noao()
    iraf.digiphot()
    iraf.daophot(_doprint = 0)
    iraf.seepsf(Pf, 'temp.psf.fits')
    PSF = fits.open('temp.psf.fits')[0].data
    system('rm temp.psf.fits')

    return PSF


def FitB(NBackground, RBackground, Pn, Pr, sn, sr, NSaturation = None, RSaturation = None, Interactive = False):
    """Fit gain matching (beta) and background (gamma) parameters"""

    center, d = NBackground.shape[0] / 2, NBackground.shape[0] / 16
    a, b = center - d, center + d
    coords = [[a, b, a, b]]

    for star in coords:
        N = NBackground[star[0]:star[1], star[2]:star[3]]
        R = RBackground[star[0]:star[1], star[2]:star[3]]

        # trim PSFs to image size
        Pn = np.roll(Pn, Pn.shape[0] / 2, 0)
        Pn = np.roll(Pn, Pn.shape[1] / 2, 1)
        Pn = Pn[star[0]:star[1], star[2]:star[3]]
        Pn = np.roll(Pn, Pn.shape[0] / 2, 0)
        Pn = np.roll(Pn, Pn.shape[1] / 2, 1)

        Pr = np.roll(Pr, Pr.shape[0] / 2, 0)
        Pr = np.roll(Pr, Pr.shape[1] / 2, 1)
        Pr = Pr[star[0]:star[1], star[2]:star[3]]
        Pr = np.roll(Pr, Pr.shape[0] / 2, 0)
        Pr = np.roll(Pr, Pr.shape[1] / 2, 1)

        print 'x: {0}, y: {1}'.format(star[0] + d, star[2] + d)
        # iteratively solve for linear fit of beta and gamma
        #if Interactive:
        #    beta, gamma, betaError, gammaError = InteractiveFit(N, R, Pn, Pr, sn, sr, NSaturation = NSaturation, RSaturation = RSaturation)
        #else:
        beta, gamma, betaError, gammaError = IterativeSolve(N, R, Pn, Pr, sn, sr, NSaturation = NSaturation, RSaturation = RSaturation, Interactive = Interactive)

    return beta, gamma, betaError, gammaError


def FitNoise(Data):
    """Find the noise of the image by fitting the background to a gaussian"""

    Edge = 50
    xMid, yMid = Data.shape[0] / 2, Data.shape[1] / 2
    TrimmedData = Data[xMid-Edge: xMid+Edge, yMid-Edge: yMid+Edge]
    TrimmedData = TrimmedData[np.where(TrimmedData < np.percentile(TrimmedData, 90))]
    Hist = np.histogram(TrimmedData, bins = 100)
    x = Hist[1][:-1]
    y = Hist[0]
    Popt, Cov = scipy.optimize.curve_fit(Gauss, x, y, p0 = [1., np.median(Data), 1.], maxfev = 1600)
    return Popt


def FitPSF(ImageFile, FWHM = 5., Noise = 30., Verbose = True, Show = True, MaxCount = 15000.):
    """Fit the PSF given an image file name"""

    if Verbose:
        verb = 'yes'
    else:
        verb = 'no'

    PSFisGood = False

    CoordsFile = ImageFile + '.coo'
    MagFile = ImageFile + '.mag'
    PSTFile = ImageFile + '.pst'
    PSFFile = ImageFile + '.psf'
    OpstFile = ImageFile + '.opst'
    GroupFile = ImageFile + '.group'
    SeeFile = ImageFile + '.see'


    while not PSFisGood:
        DeleteList = [PSFFile + '.fits', CoordsFile, MagFile, PSTFile, OpstFile, GroupFile, SeeFile]
        for item in DeleteList:
            try:
                os.remove(item)
            except:
                pass
        #system('rm {} {} {} {} {} {} {}'.format(PSFFile + '.fits', CoordsFile, MagFile, PSTFile, OpstFile, GroupFile, SeeFile))

        try:
            # generate star catalog using daofind
            iraf.noao()
            iraf.digiphot()
            iraf.daophot(_doprint = 0)
            iraf.datapars.datamax = MaxCount
            iraf.datapars.sigma = Noise
            iraf.findpars.threshold = 5.
            iraf.datapars.datamin = 0
            iraf.datapars.datamax = MaxCount
            iraf.datapars.fwhm = FWHM
            iraf.daofind(ImageFile, output = CoordsFile, verify = 'no', display = 'no', verbose = verb)

            # uncomment the raw_input line if daofind adds stars that do not exist in catalog
            # this gives you time to manually remove nonexistent stars that cause a bad psf fit
            # this is temporary until daofind works better with images coadded with swarp
            #raw_input('Manually edit .coo file now if necessary; Press enter to continue ')

            # do aperture photometry
            a1, a2, a3, a4, = int(FWHM + 0.5), int(FWHM * 2 + 0.5), int(FWHM * 3 + 0.5), int(FWHM * 4 + 0.5)
            iraf.photpars.apertures = '{0},{1},{2}'.format(a2, a3, a4)
            iraf.centerpars.calgori = 'centroid'
            iraf.fitskypars.salgori = 'mode'
            iraf.fitskypars.annulus = 10
            iraf.fitskypars.dannulu = 10
            iraf.phot(ImageFile, CoordsFile, MagFile, verify = 'no', verbose = verb)

            # select PSF stars
            iraf.daopars.fitrad = a1
            iraf.daopars.nclean = 4
            iraf.daopars.varorder = 0
            iraf.daopars.recenter = 'yes'
            iraf.pstselect(ImageFile, MagFile, PSTFile, maxnpsf = 50, verify = 'no', verbose = verb)

            # make PSF
            iraf.psf(ImageFile, MagFile, PSTFile, PSFFile, OpstFile, GroupFile, verify = 'no', verbose = verb, interactive = 'no')

            # show psf to user for approval
            if Show:
                system ('rm {}'.format(SeeFile + '.fits'))
                iraf.seepsf(PSFFile, SeeFile)
                iraf.surface(SeeFile)
                PSFisGoodyn = raw_input('GoodPSF? y/n: ')
                if PSFisGoodyn == 'y':
                    PSFisGood = True
                else:
                    FWHMguess = raw_input('New FWHM: [{}] '.format(FWHM))
                    Noiseguess = raw_input('New Noise: [{}] '.format(Noise))
                    if FWHMguess != '':
                        FWHM = float(FWHMguess)
                    if Noiseguess != '':
                        Noise = float(Noiseguess)

            else:
                break

        except:
            if Show:
                print 'PSF fitting failed; try again with different parameters'
                FWHM = float(raw_input('New FWHM: '))
                Noise = float(raw_input('New Noise: '))
            else:
                print 'Unable to fit with given parameters'
                logging.info('PSF fitting failed for {}'.format(ImageFile))
                break

    #system('rm {} {} {} {} {} {} {}'.format(CoordsFile, MagFile, PSTFile, OpstFile, GroupFile, SeeFile))

    return PSFFile


def Gamma(beta, gammaPrime, sn, sr):
    """Convert params in fourier space to image space"""

    return gammaPrime * np.sqrt(sn ** 2 + beta ** 2 * sr ** 2)


def Gauss(x, a, b, c):
    """Return a gaussian function"""

    return a * np.exp(-(x-b)**2/(2*c**2))


def InteractiveFit(N, R, Pn, Pr, sn, sr, NSaturation = None, RSaturation = None):
    N_hat = np.fft.fft2(N)    
    R_hat = np.fft.fft2(R)
    Pr_hat = np.fft.fft2(Pr)
    Pn_hat = np.fft.fft2(Pn)

    GoodFit = False
    beta = 1.
    gamma = 0.

    while not GoodFit:

        # remove saturated pixels
        DnGoodPix = CleanDataIndex(N.flatten(), SatCount = NSaturation)
        DrGoodPix = CleanDataIndex(R.flatten(), SatCount = RSaturation)
        GoodPixInCommon = np.intersect1d(DnGoodPix, DrGoodPix)
        DnFlatten = N.flatten()[GoodPixInCommon]
        DrFlatten = R.flatten()[GoodPixInCommon]

        plt.plot(DrFlatten, DnFlatten, 'bo', DrFlatten, np.polyval([beta, gamma], DrFlatten), 'r-')
        plt.show(block = False)

        GoodFityn = raw_input('Good Fit? y/[n]/b: ')
        if GoodFityn == 'y':
            GoodFit = True
        elif GoodFityn == 'b':
            print 'Bad dataset: check your images for cosmic rays, saturated stars, registration errors and bad PSFs'
            quit()
        else:
            beta = float(input('New beta: [{}] '.format(beta)))
            gamma = float(input('New gamma: [{}] '.format(gamma)))
        plt.show()
    return beta, gamma, 0, 0 # functions expect error output


def IterativeSolve(N, R, Pn, Pr, sn, sr, NSaturation = None, RSaturation = None, Interactive = False):
    """Solve for linear fit iteratively"""

    BetaTolerance = 1e-8
    GammaPrimeTolerance = 1e-8
    beta = 1.
    gammaPrime = 0.
    beta0 = 10e5
    gammaPrime0 = 10e5
    i = 0
    MaxIteration = 10

    N_hat = np.fft.fft2(N)    
    R_hat = np.fft.fft2(R)
    Pr_hat = np.fft.fft2(Pr)
    Pn_hat = np.fft.fft2(Pn)

    while abs(beta - beta0) > BetaTolerance or abs(gammaPrime - gammaPrime0) > GammaPrimeTolerance:

        SqrtDen = np.sqrt(sn ** 2 * abs(Pr_hat) ** 2 + beta ** 2 * sr ** 2 * abs(Pn_hat) ** 2)
        Dn_hat = Pr_hat * N_hat / SqrtDen
        Dr_hat = Pn_hat * R_hat / SqrtDen
        Dn = np.real(np.fft.ifft2(Dn_hat))
        Dr = np.real(np.fft.ifft2(Dr_hat))
        DnFlatten = Dn.flatten()
        DrFlatten = Dr.flatten()

        # remove saturated pixels
        DnGoodPix = CleanDataIndex(DnFlatten, SatCount = NSaturation, Significance = 3.)
        DrGoodPix = CleanDataIndex(DrFlatten, SatCount = RSaturation, Significance = 3.)
        GoodPixInCommon = np.intersect1d(DnGoodPix, DrGoodPix)
        DnFlatten = DnFlatten[GoodPixInCommon]
        DrFlatten = DrFlatten[GoodPixInCommon]

        beta0, gammaPrime0 = beta, gammaPrime

        x = sm.add_constant(DrFlatten)
        y = DnFlatten


        RobustFit = sm.RLM(y, x).fit()
        Parameters = RobustFit.params
        Errors = RobustFit.bse
        beta = Parameters[1]
        gammaPrime = Parameters[0]
        betaError = Errors[1]
        gammaPrimeError = Errors[0]

        if i == MaxIteration: break
        i += 1
        print 'Iteration {}:'.format(i)
        print 'Beta = {0}, gamma = {1}'.format(beta, Gamma(beta, gammaPrime, sn, sr))

    print 'Fit done in {} iterations'.format(i)

    Cov = RobustFit.bcov_scaled[0,0]

    if Interactive:
        plt.plot(DrFlatten, DnFlatten, 'bo', DrFlatten, RobustFit.fittedvalues, 'r-')
        plt.show()
    #gammaPrime = FitNoise(y)[1] - FitNoise(x)[1]
    gamma = Gamma(beta, gammaPrime, sn, sr)


    if Cov > 1000.:
        print 'The is a poor fit. Try fitting the PSFs manually.'
        beta = 1.
        gamma = FitNoise(N)[1] - FitNoise(R)[1]
        betaError = 1.
        gammaError = 1.
        logging.info('Parameter fitting for images failed')
        raise ValueError

    #plt.plot(DrFlatten, DnFlatten, 'bo', DrFlatten, RobustFit.fittedvalues, 'r-')
    #plt.show()
    print 'Beta = ' + str(beta)
    print 'Gamma = ' + str(gamma)
    print 'Beta Variance = ' +str(Cov)
    return beta, gamma, betaError, gammaPrimeError


def Normalize(Image):
    """Normalize to sum = 1"""

    return Image / np.sum(Image)


def PadPSF(PSF, shape):
    """Pad PSF and center at (0,0)"""

    p = PSF.shape
    s = shape
    PSF_extended = np.zeros(s)
    PSF_extended[s[0] / 2 - p[0] / 2 - 1:s[0] / 2 + p[0] / 2, s[1] / 2 - p[1] / 2 - 1:s[1] / 2 + p[1] / 2] = PSF
    PSF_extended = np.roll(PSF_extended, s[0] / 2, 0)
    PSF_extended = np.roll(PSF_extended, s[1] / 2, 1)

    return PSF_extended

