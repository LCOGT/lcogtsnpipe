#!/usr/bin/env python

import numpy as np
import scipy
from astropy.io import fits
from pyraf import iraf
import time
from os import system
import argparse
from lsc.lscastrodef import sextractor
import statsmodels.api as sm
import logging

from SubtractionUtil import *


import matplotlib.pyplot as plt

class OptimalSubtraction:
    '''Implementation of proper image subtraction from Zackey, Ofek, Gal-Yam 2016
       see https://arxiv.org/abs/1601.02655 for details'''


    def __init__(self, Nf, Rf, ArgsDict = {}):
        '''Take filenames and turn them into arrays and unpack arguments in ArgsDict'''

        logging.basicConfig(filename = Nf.replace('.fits', '.log'))

        self.ArgsDict = ArgsDict
        self.Nf, self.Rf = Nf, Rf
        self.MakeDefaults()
        

        #self.Nf, self.Rf = Nf, Rf

        self.CheckPSFsExist()

        self.CheckPSFfromIRAF()
        self.CheckAlign()

        self.N, self.R, self.Pn, self.Pr = CropImages(fits.getdata(self.Nf), fits.getdata(self.Rf), self.Pn, self.Pr)
        self.CheckBackground()

        # work in progress
        #self.CheckFlatten()

        self.sn = self.ArgsDict['NewNoise']
        self.sr = self.ArgsDict['RefNoise']

        self.CheckGainMatched()

    def MakeDefaults(self):
        '''Set up default optional arguments'''

        default = {
            'Align': False, 
            'Beta': 99999., 
            'Interactive': False, 
            'Flatten': False, 
            'Gamma': 99999., 
            'RefPSF': '', 
            'NewPSF': '', 
            'RefBackground': self.Rf, 
            'RefFWHM': 5., 
            'RefNoise': FitNoise(fits.getdata(self.Rf))[2], 
            'RefReadNoise': 0., 
            'RefSaturation': 10e6, #GetSaturationCount(self.Rf), 
            'RefVariance': '', 
            'NewBackground': self.Nf, 
            'NewFWHM': 5., 
            'NewNoise': FitNoise(fits.getdata(self.Nf))[2],
            'NewReadNoise': 0., 
            'NewSaturation': 10e6, #GetSaturationCount(self.Nf), 
            'NewVariance': '', 
            'PSFfromIRAF': True, 
            'SigmaX': 0.5, 
            'SigmaY': 0.5, 
            'Verbose': False
        }

        default.update(self.ArgsDict)
        self.ArgsDict = default


    def CheckAlign(self):
        '''Check if images need alignment'''

        if self.ArgsDict['Align']:
            self.R = AlignImages(self.Nf, self.Rf)
            self.Rf = self.Nf.replace('.fits', '.ref.fits')


    def CheckBackground(self):
        try:
            self.NBackground = fits.getdata(self.ArgsDict['NewBackground'])
        except ValueError:
            self.NBackground = self.N
        try:
            self.RBackground = fits.getdata(self.ArgsDict['RefBackground'])
        except ValueError:
            self.RBackground = self.R


    def CheckFlatten(self):
        if self.ArgsDict['Flatten']:
            self.N = SubtractBackground(self.N)
            self.R = SubtractBackground(self.R)


    def CheckGainMatched(self):
        '''Check if the gain matching parameters have been found'''

        beta = self.ArgsDict['Beta']
        gamma = self.ArgsDict['Gamma']
        # if user hasn't supplied beta and gamma, fit them
        if beta == 99999 or gamma == 99999:
            print 'Gain matching not done; fitting manually'

            self.beta, self.gamma = self.MatchGain()

        else:
            self.beta = beta
            self.gamma = gamma


    def CheckPSFsExist(self):
        '''Check if user provided a psf and make one if necessary'''

        Prf = self.ArgsDict['RefPSF']
        Pnf = self.ArgsDict['NewPSF']
        if Pnf != '':
            self.Pnf = Pnf
        else:
            FWHM = self.ArgsDict['NewFWHM']
            Noise = self.ArgsDict['NewNoise']
            self.Pnf = FitPSF(self.Nf, FWHM, Noise, Verbose = self.ArgsDict['Verbose'], Show = self.ArgsDict['Interactive'])
            self.ArgsDict.update({'PSFfromIRAF': True})

        if Prf != '':
            self.Prf = Prf
        else:
            FWHM = self.ArgsDict['RefFWHM']
            Noise = self.ArgsDict['RefNoise']
            self.Prf = FitPSF(self.Rf, FWHM, Noise, Verbose = self.ArgsDict['Verbose'], Show = self.ArgsDict['Interactive'])
            self.ArgsDict.update({'PSFfromIRAF': True})
            

    def CheckPSFfromIRAF(self):
        '''Check if user specified IRAF psf or actual psf'''

        try:
            PSFfromIRAF = self.ArgsDict['PSFfromIRAF']
            if PSFfromIRAF:
                self.Pn = ExtractPSF(self.Pnf)
                self.Pr = ExtractPSF(self.Prf)
            else:
                self.Pn = fits.getdata(self.Pnf)
                self.Pr = fits.getdata(self.Prf)
        except KeyError:
            self.Pn = fits.getdata(self.Pnf)
            self.Pr = fits.getdata(self.Prf)


    def D(self, normalize = '', diagnostic = False):
        '''Calculate proper subtraction image and normalize to zero point of reference or target'''

        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr

        #save convolved not gain matched images that will be subtracted
        if diagnostic:
            print 'Saving Deconvolved images'
            N_hat = np.fft.fft2(self.N)    
            R_hat = np.fft.fft2(self.R)
            Pr_hat = np.fft.fft2(self.Pr)
            Pn_hat = np.fft.fft2(self.Pn)
            Den = self.sr ** 2 * abs(Pn_hat) ** 2 + self.sn ** 2 * abs(Pr_hat) ** 2
            SqrtDen = np.sqrt(Den)
            DConvN = np.real(np.fft.ifft2(N_hat * Pr_hat / SqrtDen))
            hdu = fits.PrimaryHDU(DConvN)
            hdu.writeto('NDConvolve.fits', clobber = True)
            DConvR = np.real(np.fft.ifft2(R_hat * Pn_hat / SqrtDen))
            hdu = fits.PrimaryHDU(DConvR)
            hdu.writeto('RDConvolve.fits', clobber = True)


        # calculate D
        N = np.subtract(N, self.gamma)

        beta = self.beta

        N_hat = np.fft.fft2(N)
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2
        SqrtDen = np.sqrt(Den)

        D = np.fft.ifft2((Pr_hat * N_hat - beta * Pn_hat * R_hat) / SqrtDen)

        # apply user's normalization choice
        if normalize == 'i':
            DNormalize = D * beta / self.Fd()
        elif normalize == 't':
            DNormalize = D / self.Fd()
        else:
            DNormalize = D
        self.D_ = np.real(DNormalize)

        return np.real(DNormalize)


    def DNoise(self):
        a = self.N.shape
        deltaD = deltaD_hat * abs(np.divide(D_hat[a,a] * np.exp(2*np.pi*1j*np.dot(a,a)) - D_hat[0, 0], 2*np.pi*1j*D_hat))


    def Fd(self):
        '''Calculate the flux based zero point of D'''

        Fd = self.beta / np.sqrt(self.sn**2 + self.sr**2*self.beta**2)
        self.Fd_ = np.real(Fd)
        return np.real(Fd)


    def FindTransient(self, Threshold = 3., filename = 'transients.txt'):
        '''Write transient detections to file'''

        try:
            self.Scorr_
        except AttributeError:
            self.Scorr()
        edge = 0
        Scorr_ = self.Scorr_[edge: self.Scorr_.shape[0] - edge, edge: self.Scorr_.shape[1] - 100]
        maxima = (Scorr_ == scipy.ndimage.maximum_filter(Scorr_, 4))
        PeaksIndex = np.where(maxima)
        Peaks = Scorr_[PeaksIndex]
        Transients = np.array([Peaks, PeaksIndex[0], PeaksIndex[1]]).transpose()
        TransientsSortIndex = np.argsort(Transients[:,0])
        Transients = Transients[TransientsSortIndex[::-1]]
        Transients = Transients[np.where(Transients[:,0] >= Threshold)]
        np.savetxt(filename , Transients)


    def Flux(self, normalize = ''):
        '''Calculate transient Flux'''

        try:
            Flux = self.S_ / self.Fs_
        except AttributeError:
            Flux = self.S(normalize) / self.Fs()
        self.Flux_ = Flux
        return Flux


    def Fs(self):
        '''Calculate flux based zeropoint of S'''

        beta = self.beta
        Pr_hat = np.fft.fft2(self.Pr)
        Pn_hat = np.fft.fft2(self.Pn)
        sn, sr = self.sn, self.sr
        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2
        Fs = np.sum(beta ** 2 * abs(Pn_hat) * abs(Pr_hat) / Den)
        self.Fs_ = Fs
        return Fs

    def Kernels(self):
        '''Calculate convolution kernels'''

        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr
        N = np.subtract(N, self.gamma)

        N_hat = np.fft.fft2(N)    
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        beta = self.beta
        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2

        Kr_hat = (beta ** 2 * np.conj(Pr_hat) * abs(Pn_hat) ** 2) / Den
        Kr = np.real(np.fft.ifft2(Kr_hat))

        Kn_hat = (beta * np.conj(Pn_hat) * abs(Pr_hat) ** 2) / Den
        Kn = np.real(np.fft.ifft2(Kn_hat))
        return Kn, Kr


    def MakeCatalog(self, SortBy = 'magnitude'):
        '''Check for source catalog and make one if necessary'''

        try:
            cat = self.Catalog
            return cat

        except AttributeError:
            if SortBy == 'x':
                SortIndex = 0
            elif SortBy == 'y':
                SortIndex = 1
            elif SortBy == 'fwhm':
                SortIndex = 2
            elif SortBy == 'flux':
                SortIndex = 3
            elif SortBy == 'magnitude':
                SortIndex = 4
            elif SortBy == 'flag':
                SortIndex = 7
            else:
                SortIndex = 4

            sources = np.array(sextractor(self.Nf)).transpose()
            source_index = np.argsort(sources[:, SourceIndex])
            source_sorted = sources[source_index]

            self.Catalog = source_sorted
            return source_sorted


    def MatchGain(self):
        '''Call gain matching function'''

        NSat = self.ArgsDict['RefSaturation']
        RSat = self.ArgsDict['NewSaturation']
        beta, gamma, betaError, gammaError = FitB(self.N, self.R, self.Pn, self.Pr, self.sn, self.sr, 
                                                  NSaturation = NSat, RSaturation = RSat, Interactive = self.ArgsDict['Interactive'])

        return beta, gamma


    def Pd(self):
        '''Calculate PSF of D'''

        Pn_hat = np.fft.fft2(self.Pn)
        Pr_hat = np.fft.fft2(self.Pr)
        Den = self.sr ** 2 * self.beta ** 2 * abs(Pn_hat) ** 2 + self.sn ** 2 * abs(Pr_hat) ** 2
        SqrtDen = np.sqrt(Den)
        Pd_hat = self.beta * Pr_hat * Pn_hat / (self.Fd() * SqrtDen)
        Pd = np.fft.ifft2(Pd_hat)
        self.Pd_ = np.real(Pd)

        return np.real(Pd)


    def S(self):
        '''Calculate matched filter image S'''

        try:
            S_hat = self.Fd_ * np.fft.fft2(self.D_) * np.conj(np.fft.fft2(self.Pd_))
            S = np.fft.ifft2(S_hat)
        except AttributeError:
            S_hat = self.Fd() * np.fft.fft2(self.D()) * np.conj(np.fft.fft2(self.Pd()))
            S = np.fft.ifft2(S_hat)
        self.S_ = np.real(S)
        return np.real(S)

    def SaveD(self, filename, normalize = ''):
        '''Calculate and save proper subtraction image to database'''

        self.Df = filename
        self.D(normalize)
        hdu = fits.PrimaryHDU(np.real(self.D_))
        #hdu.header = fits.getheader(self.Nf)
        #print normalize, self.beta, self.gamma
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.header['CONVOL00'] = normalize
        hdu.writeto(self.Df, clobber = True)

    def SaveS(self, Sf, normalize = ''):
        '''Calculate and save S to database'''

        self.Sf = Sf
        self.S()
        hdu = fits.PrimaryHDU(np.real(self.S_))
        #hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.Sf, clobber = True)

    def SaveScorr(self, Scorrf, normalize = ''):
        '''Calculate and save S to database'''

        self.Scorrf = Scorrf
        self.Scorr()
        hdu = fits.PrimaryHDU(np.real(self.Scorr_))
        #hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.Scorrf, clobber = True)

    def SaveScorrThreshold(self, ScorrThreshf, Thresh = 3., normalize = ''):
        '''Save a Scorr such that pixel = 1 if > Thresh else pix = 0'''

        self.ScorrThreshf = ScorrThreshf
        self.Scorr()
        BrightIndex = np.where(self.Scorr_ >= Thresh)
        DarkIndex = np.where(self.Scorr_ < Thresh)
        ScorrThresh = np.copy(self.Scorr_)
        ScorrThresh[BrightIndex] = 1
        ScorrThresh[DarkIndex] = 0
        self.ScorrThresh = ScorrThresh        

        hdu = fits.PrimaryHDU(np.real(self.ScorrThresh))
        #hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.ScorrThreshf, clobber = True)


    def SaveSnoise(self, Snoisef, normalize = ''):
        '''Calculate and save S to database'''

        self.Snoisef = Snoisef
        self.Snoise()
        hdu = fits.PrimaryHDU(np.real(self.Snoise_))
        #hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.Snoisef, clobber = True)

    def SaveImageToWD(self):
        '''Save various images to working directory (testing only)'''

        Images = {'S.fits': self.S, 'Snoise.fits': self.Snoise, 'Scorr.fits': self.Scorr, 'D.fits': self.D, 'Flux.fits': self.Flux, 'Pd.fits': self.Pd}

        for element in Images:
            hdu = fits.PrimaryHDU(np.real(Images[element]))
            hdu.header=fits.getheader(self.Nf)
            hdu.writeto(element, clobber = True)


    def Scorr(self):
        '''Calculate Scorr'''

        try:
            Scorr = self.S_ / self.Snoise_
        except AttributeError:
            Scorr = self.S() / np.sqrt(self.Snoise())
        self.Scorr_ = np.real(Scorr)
        return np.real(Scorr)        


    def Snoise(self):
        '''Calculate the noise image for Scorr'''

        # this whole function needs optimization

        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr
        N = np.subtract(N, self.gamma)

        N_hat = np.fft.fft2(N)    
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        beta = self.beta
        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2

        Kn, Kr = self.Kernels()

        Kn2_hat = np.fft.fft2(Kn ** 2)
        Kr2_hat = np.fft.fft2(Kr ** 2)

        if self.ArgsDict['NewVariance'] != '':
            NVarFilename = self.ArgsDict['NewVariance']
            EpsN = fits.getdata(NVarFilename)
            s = N.shape
            EpsN = EpsN[:s[0], :s[1]]
            V_Sn = np.fft.ifft2(Kn2_hat * np.fft.fft2(EpsN))

        else:
            NReadNoise = self.ArgsDict['NewReadNoise']
            EpsN = np.add(abs(self.NBackground), NReadNoise ** 2)
            V_Sn = np.fft.ifft2(Kn2_hat * np.fft.fft2(EpsN))

        if self.ArgsDict['RefVariance'] != '':
            RVarFilename = self.ArgsDict['RefVariance']
            EpsR = fits.getdata(RVarFilename)
            s = N.shape
            EpsR = EpsR[:s[0], :s[1]]
            V_Sr = np.fft.ifft2(Kr2_hat * np.fft.fft2(EpsR))

        else:
            RReadNoise = self.ArgsDict['RefReadNoise']
            EpsR = np.add(abs(self.RBackground), RReadNoise ** 2)
            V_Sr = np.fft.ifft2(Kr2_hat * np.fft.fft2(EpsR))

        xrms = self.ArgsDict['SigmaX']
        yrms = self.ArgsDict['SigmaY']

        GradNy, GradNx = np.gradient(np.fft.ifft2(np.fft.fft2(Kn) * N_hat))
        GradRy, GradRx = np.gradient(np.fft.ifft2(np.fft.fft2(Kr) * R_hat))

        Vr_ast = xrms ** 2 * GradRx ** 2 + yrms ** 2 * GradRy ** 2
        Vn_ast = xrms ** 2 * GradNx ** 2 + yrms ** 2 * GradNy ** 2

        Snoise = V_Sn + V_Sr + Vr_ast + Vn_ast
        self.Snoise_ = np.real(Snoise)
        return np.real(Snoise)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform optimal image subtraction')
    parser.add_argument('-N', dest = 'input', help = 'New background subtracted image', required = True)
    parser.add_argument('-R', dest = 'template', help = 'Reference background subtracted image', required = True)
    parser.add_argument('--NewPSF', dest = 'input_PSF', default = '', help = 'New PSF')
    parser.add_argument('--RefPSF', dest = 'template_PSF', default = '', help = 'Reference FWHM')
    parser.add_argument('--NewFWHM', dest = 'input_FWHM', default = '', help = 'New FWHM')
    parser.add_argument('--RefFWHM', dest = 'template_FWHM', default = '', help = 'Reference PSF')
    parser.add_argument('-o', dest = 'output', help = 'Output image')
    parser.add_argument('--NewBackground', dest = 'NewBackground', help = 'New image with background')
    parser.add_argument('--RefBackground', dest = 'RefBackground', help = 'Reference image with background')
    parser.add_argument('--Mode', dest = 'mode', default = 'D', help = 'Subtraction mode')
    #parser.add_argument('--threshold', dest = 'threshold', help = 'Sigma threshold for detection')
    
    parser.add_argument('--Interactive', dest = 'interactive', action = 'store_true', default = False, help = 'Perform subtraction interactively')
    parser.add_argument('--Align', dest = 'align', action = 'store_true', default = False, help = 'Align images before subtraction')
    parser.add_argument('--Verbose', dest = 'verbose', action = 'store_true', default = False, help = 'Show IRAF output in PSF fitting')
    parser.add_argument('--Beta', dest = 'beta', default = 99999, help = 'Gain matching parameter beta')
    parser.add_argument('--Gamma', dest = 'gamma', default = 99999, help = 'Gain matching parameter gamma')
    parser.add_argument('--NewNoise', dest = 'NewNoise', default = 5, help = 'Standard deviation of new image background')
    parser.add_argument('--RefNoise', dest = 'RefNoise', default = 5, help = 'Standard deviation of ref image background')
    parser.add_argument('--Flatten', dest = 'flatten', action = 'store_true', default = False, help = 'Fix Background locally')
    parser.add_argument('--Normalize', default = '', dest = 'normalize', help = 'Normalize to which image')
    args = parser.parse_args()

    d = {
        'Align': args.align, 
        'Flatten': args.flatten, 
        'NewPSF': args.input_PSF, 
        'RefPSF': args.input_PSF, 
        'NewFWHM': args.input_FWHM, 
        'RefFWHM': args.template_FWHM, 
        'Interactive': args.interactive, 
        'Beta': float(args.beta), 
        'Gamma': float(args.gamma), 
        'Verbose': args.verbose, 
        'NewBackground': args.NewBackground, 
        'RefBackground': args.RefBackground, 
        'NewNoise': float(args.NewNoise), 
        'RefNoise': float(args.RefNoise)
        }

    if args.mode == 'D':
        OptimalSubtraction(args.input, args.template, d).SaveD(args.output, args.normalize)
    elif args.mode == 'Scorr':
        OptimalSubtraction(args.input, args.template, d).SaveScorr(args.output, args.normalize)
    elif args.mode == 'S':
        OptimalSubtraction(args.input, args.template, d).SaveS(args.output, args.normalize)
    elif args.mode == 'Thresh':
        OptimalSubtraction(args.input, args.template, d).SaveScorrThreshold(args.output, normalize = args.normalize, Thresh = args.threshold)
    elif args.mode == 'Find':
        OptimalSubtraction(args.input, args.template, d).FindTransient()
    else:
        pass

