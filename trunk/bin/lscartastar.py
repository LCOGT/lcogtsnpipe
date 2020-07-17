#!/usr/bin/env python
description = ">> add artificial stars "
usage = "%prog image [options] "

from optparse import OptionParser
import string
import re
import os
import sys
import glob
import time
import lsc
import numpy as np
import pylab as pl
from astropy.io import fits as pyfits


##################################################################3

def runsex(img, fwhm, thresh, pix_scale, mina=5.):  ## run_sextractor  fwhm arcsec
    temp_file0 = 'sex_'+next(tempfile._get_candidate_names())

    # seeing should be in pixels
    # fixed a bug fhwm (arcsec) / pixelscale (arcsec/pixel)
    seeing = fwhm / pix_scale

    cdef = open(lsc.__path__[0] + '/standard/sex/default2.param')
    riga = cdef.readlines()
    cparam = []
    for r in riga:
        if r[0] != '#' and len(r.strip()) > 0: \
                cparam.append(r.split()[0])

    pid = subprocess.Popen("sex " + img + " -catalog_name " + temp_file0 + \
                           " -c  " + lsc.__path__[0] + '/standard/sex/default2.sex' \
                                                       " -PARAMETERS_NAME " + lsc.__path__[
                               0] + "/standard/sex/default2.param" + \
                           " -STARNNW_NAME " + lsc.__path__[0] + "/standard/sex/default2.nnw" + \
                           " -PSF_NAME     " + lsc.__path__[0] + "/standard/sex/default2.psf" + \
                           " -PIXEL_SCALE " + str(pix_scale) + \
                           " -DETECT_MINAREA " + str(mina) + \
                           " -DETECT_THRESH  " + str(thresh) + \
                           " -ANALYSIS_THRESH  " + str(thresh) + \
                           " -PHOT_FLUXFRAC 0.5" + \
                           " -SEEING_FWHM " + str(seeing),
                           stdout=subprocess.PIPE, shell=True)

    output, error = pid.communicate()

    csex = open(temp_file0)
    tab = {}
    riga = csex.readlines()
    for k in cparam: tab[k] = []
    for r in riga:
        if r[0] != '#':
            for i in range(len(cparam)):
                tab[cparam[i]].append(float(r.split()[i]))
    for k in cparam: tab[k] = np.array(tab[k])

    lsc.util.delete(temp_file0)
    return tab


# #######################################################################

if __name__ == "__main__":
    start_time = time.time()
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-R", "--RA", dest="ra", default='', type="str", help='-R  ra    \t [%default]')
    parser.add_option("-D", "--DEC", dest="dec", default='', type="str", help='-D dec   \t [%default]')
    parser.add_option("-M", "--MAG", dest="mag", default='', type="str", help='-M mag   \t [%default]')
    parser.add_option("-n", "--num", dest="num", default='1', type="str", help='-n num   \t [%default]')
    parser.add_option("-v", "--verbose", action="store_true",
                      dest='verbose', default=False, help='verbose \t\t\t [%default]')

    option, args = parser.parse_args()
    if len(args) < 1: 
        sys.argv.append('--help')
    _ra = option.ra
    _dec = option.dec
    _mag = option.mag
    _num = option.num
    _verbose  = option.verbose
    if _mag: 
        _mag = float(_mag)
    else:
        _mag = -5
    if _ra:
        _ra = string.split(',',_ra)
    if _dec:
        _dec = string.split(',',_dec)


    option, args = parser.parse_args()
    imglist = lsc.util.readlist(args[0])
    for img in imglist:
      _input = []
      _output = []
      if _verbose:
          _mag = [13,14,15,16,17,18,19,13,14,15,16,17,18,19,18,18,18,18]
      else:
          _mag = [_mag]

      for mag in _mag:
        imgart = re.sub('.art.fits','.fits',img)
        _r,_d,_x,_y,_m = lsc.lscartstardef.artstar(img, _ra, _dec, int(_num), mag, _verbose=True)

        if _verbose:
            hdr1 = pyfits.getheader(img)
            _exptime = hdr1['exptime']
            _zeropoint = hdr1['ZP']
            artfile = re.sub('.fits','.art.fits',img)
            lsc.myloopdef.run_cosmic([artfile], 'photlco', 4.5, 0.2, 4, False)
            #   add zeropoint stage
#            dlt40.dlt40loopdef.run_zero([artfile], 'photlco', '', '', '', False)

            fwhm =5
            _threshold = 5
            ####I need to change this number
            pixelscale = 0.59 
            tab = runsex(artfile, fwhm, _threshold, pixelscale)
            fwhm = np.median(tab['FWHM_IMAGE'])
            print fwhm
            tab = runsex(artfile, fwhm, _threshold, pixelscale)

            xx = tab['X_IMAGE']
            yy = tab['Y_IMAGE']
            aaa = lsc.lscastrodef.crossmatchxy(np.array(xx),np.array(yy), _x, _y,10)

            print 'magnituded added'
            for jj,kk in enumerate(_x):
                print '%g %g %g %g %g' % (_r[jj],_d[jj],_x[jj],_y[jj],_m[jj])
            if len(aaa[0]):
                id = aaa[1][0]
                print 'magnitude measure on the image'
                print tab['X_WORLD'][id],tab['Y_WORLD'][id], tab['X_IMAGE'][id], tab['Y_IMAGE'][id], tab['MAG_AUTO'][id]
                _input.append(_m[0])
                _output.append(tab['MAG_AUTO'][id])

    if _verbose:
        _input = np.array(_input) #- np.log10(_exptime) * 2.5 + _zeropoint
        _output = np.array(_output) + _zeropoint
        pl.ion()
        pl.plot(_input,_output,'ro')
        pl.plot(_input,_input,'b-')

