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
    '''Align reference image to new image when both images have WCS'''

    TempFile = NewFile.replace('.fits', '.ref.fits')
    iraf.wregister(RefFile, NewFile, TempFile)
    #system('cp {0} {1}'.format(TempFile, NewFile.replace('.fits', '.ref.fits')))
    AlignedRef = fits.getdata(TempFile)
    return AlignedRef


def CleanDataIndex(Data, SatCount = None, RmBackPix = True, Significance = 1.):
    '''Clean saturated and background pixels from dataset'''

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


def ConvertWeightToVariance(image):
    '''Convert swarp weight image to variance; taken from externaldata.py'''

    # the output of swarp give the normalized weight 
    # we want to store the variance image
    # we invert the normalized weight and we "de-normalized" this image
    hd = fits.getheader(image)
    ar = fits.getdata(image)
    hd2 = fits.getheader(image.replace('.fits', '.weight.fits'))
    ar2 = fits.getdata(image.replace('.fits', '.weight.fits'))
    variance = 1 / ar2
    #  this is to take in account that the weight is normalized
    variance *= (np.median(np.abs(ar - np.median(ar)))*1.48)**2/np.median(variance)
    varimg = image.replace('.fits', '.var.fits')
    _saturate = GetSaturationCount(image)
    ar = np.where(ar2 == 0, _saturate, ar)

    fits.writeto(varimg, variance, hd2, clobber=True)

    # put the saturation all values where the weight is zero 

def CropImages(N, R, Pn, Pr):
    '''Open and crop images and PSFs to the same size'''

    sN = N.shape
    #sR = R.shape
    #s = [np.min([sN[0], sR[0]]), np.min([sN[1], sR[1]])]
    #N = N[:s[0], :s[1]]
    #R = R[:s[0], :s[1]]

    Pn_extended = PadPSF(Pn, sN)
    Pr_extended = PadPSF(Pr, sN)

    return N, R, Normalize(Pn_extended), Normalize(Pr_extended)


def ExtractPSF(Pf):
    '''Extract PSF array from iraf PSF file'''

    iraf.noao()
    iraf.digiphot()
    iraf.daophot(_doprint = 0)
    iraf.seepsf(Pf, 'temp.psf.fits')
    PSF = fits.open('temp.psf.fits')[0].data
    system('rm temp.psf.fits')

    return PSF


def FitB(NBackground, RBackground, Pn, Pr, sn, sr, NSaturation = None, RSaturation = None, Interactive = False):
    '''Fit gain matching (beta) and background (gamma) parameters'''

    center, d = NBackground.shape[0] / 2, NBackground.shape[0]
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
    '''Find the noise of the image by fitting the background to a gaussian'''

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
    '''Fit the PSF given an image file name'''

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
    '''Convert params in fourier space to image space'''

    return gammaPrime * np.sqrt(sn ** 2 + beta ** 2 * sr ** 2)


def Gauss(x, a, b, c):
    '''Return a gaussian function'''

    return a * np.exp(-(x-b)**2/(2*c**2))


def GetSaturationCount(filename):
    '''Get pixel saturation count for a given image'''

    sat = []
    Header = fits.getheader(filename)
    try:
        MaxLin = Header['MAXLIN']
        sat.append(MaxLin)
    except KeyError:
        pass
    try:
        Saturate = Header['SATURATE']
        sat.append(Saturate)
    except KeyError:
        pass
    if len(sat) < 1:
        print 'Saturation count not found in fits header'
        return None
    else:
        Saturation = min(sat)
        return Saturation


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
    '''Solve for linear fit iteratively'''

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

    gamma = Gamma(beta, gammaPrime, sn, sr)


    if Cov > 1000.:
        print 'The is a poor fit. Try fitting the PSFs manually. For now, Beta = 1'
        beta = 1.
        gamma = FitNoise(N)[1] - FitNoise(R)[1]
        betaError = 1.
        gammaError = 1.
        logging.info('Parameter fitting for images failed')
        raise ValueError

    print 'Beta = ' + str(beta)
    print 'Gamma = ' + str(gamma)
    print 'Beta Variance = ' +str(Cov)
    return beta, gamma, betaError, gammaPrimeError


def Normalize(Image):
    '''Normalize to sum = 1'''

    return Image / np.sum(Image)


def PadPSF(PSF, shape):
    '''Pad PSF and center at (0,0)'''

    p = PSF.shape
    s = shape
    PSF_extended = np.zeros(s)
    PSF_extended[s[0] / 2 - p[0] / 2 - 1:s[0] / 2 + p[0] / 2, s[1] / 2 - p[1] / 2 - 1:s[1] / 2 + p[1] / 2] = PSF
    PSF_extended = np.roll(PSF_extended, s[0] / 2, 0)
    PSF_extended = np.roll(PSF_extended, s[1] / 2, 1)

    return PSF_extended


def RemoveSaturatedFromCatalog(catalog):
    '''Remove stars with saturated pixels from catalog'''

    SafeIndex = np.where(catalog[:7].flatten() == 0)
    SafeSource = catalog[SafeIndex]
    return SafeSource


def SubtractBackground(data, NumStamps = 30, Border = 200):
    '''Subtract spatially varying background'''

    d = (data.shape[0] - 2 * Border) * (data.shape[1] - 2 * Border) / NumStamps ** 2
    LeftEdge = np.linspace(Border, data.shape[0] - Border, NumStamps)
    LeftEdge = [int(i) for i in LeftEdge]
    TopEdge = np.linspace(Border, data.shape[1] - Border, NumStamps)
    TopEdge = [int(i) for i in TopEdge]
    SmallImg = np.zeros(2 * [NumStamps])
    for iIndex, i in enumerate(TopEdge):
        for jIndex, j in enumerate(LeftEdge):
            #print i, i+d, j, j+d
            Stamp = data[i: i + d, j: j + d]
            Hist, BinEdge = np.histogram(Stamp, bins = 50)
            x = [(BinEdge[k] + BinEdge[k+1]) / 2 for k in range(len(BinEdge) - 1)]
            popt, cov = scipy.optimize.curve_fit(Gauss, x, Hist, p0 = [np.max(data), np.median(data), np.std(data)])
            #print popt
            SmallImg[jIndex, iIndex] = popt[1]
            #plt.imshow(Stamp)
            #plt.show()
            #plt.plot(x, Hist, x, Gauss(x, *popt))
            #plt.show()
    scale = [(data.shape[0] - 2 * d) / SmallImg.shape[0], (data.shape[1] - 2 * d) / SmallImg.shape[1]]
    print scale
    BackImg = scipy.ndimage.zoom(SmallImg, scale)
    s = [data.shape[0]/2 - BackImg.shape[0]/2, data.shape[1]/2 - BackImg.shape[1]/2]
    BackImgFull = np.rot90(np.fliplr(np.pad(BackImg, ((s[0], s[0]),(s[1], s[1])), 'constant')), 1)
#    BackImgFull = np.zeros(data.shape)
#    BackImgFull[Border:data.shape[0]-Border, Border:data.shape[1]-Border] = BackImg
    dataCorr = data - BackImgFull
    fits.writeto('back.fits', BackImgFull, clobber = True)
    fits.writeto('dataCorr.fits', dataCorr, clobber = True)

    plt.imshow(BackImg, interpolation = 'none')
    plt.show()
    plt.imshow(dataCorr, interpolation = 'none')
    plt.show()
    return dataCorr
