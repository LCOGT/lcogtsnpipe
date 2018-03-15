#!/usr/bin/env python
description = ">> New automated psf version"
usage = "%prog image [options] "

import lsc
import time
import sys
import os
from pyraf import iraf
from optparse import OptionParser

if __name__ == "__main__":
    start_time = time.time()
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-f", "--fwhm", dest="fwhm", type=float, help='starting FWHM  \t\t\t %default')
    parser.add_option("-t", "--threshold", dest="threshold", default=10., type='float',
                      help='Source detection threshold \t\t\t %default')
    parser.add_option("-p", "--psfstars", dest="psfstars", default=6, type='int',
                      help='Maximum number of psf stars \t\t\t %default')
    parser.add_option("-d", "--distance", dest="distance", default=5,
                      type='int', help='Minimum star separation (in unit of FWHM) \t\t %default')
    parser.add_option("--function", dest="psffun", default='gauss', type='str',
                      help='psf function (gauss,auto,moffat15,moffat25,penny1,penny2,lorentz) \t\t %default')
    parser.add_option("-r", "--redo", action="store_true", dest='redo', default=False,
                      help='Re-do \t\t\t\t [%default]')
    parser.add_option("-i", "--interactive", action="store_true", dest='interactive', default=False,
                      help='Interactive \t\t\t [%default]')
    parser.add_option("-s", "--show", dest="show", action='store_true',
                      default=False, help='Show PSF output \t\t [%default]')
    parser.add_option("--use-sextractor", action="store_true", help="use souces from sextractor instead of catalog")
    parser.add_option("-c", "--catalog", dest="catalog", default='', type='str',
                      help='use input catalog  \t\t %default')
    parser.add_option("--fix", action="store_true", dest='fixaperture', default=False,
                      help='fixaperture \t\t\t [%default]')
    parser.add_option("--datamin", type=int, default=-100,
                      help="value below which pixels are considered bad")
    parser.add_option("--datamax", type=int,
                      help="value above which pixels are considered saturated (default = SATURATE in header)")

    option, args = parser.parse_args()
    if len(args) < 1: 
        sys.argv.append('--help')
    option, args = parser.parse_args()
    imglist = lsc.util.readlist(args[0])
    _catalog = option.catalog
    fixaperture = option.fixaperture
    psffun = option.psffun
    if psffun not in ['gauss', 'auto', 'lorentz', 'moffat15', 'moffat25', 'penny1', 'penny2']:
        sys.argv.append('--help')
    option, args = parser.parse_args()

    for img in imglist:
        if os.path.exists(img.replace('.fits', '.psf.fits')) and not option.redo:
            print img + ': psf already calculated'
        else:
            if 'optimal' in img: # PyZOGY difference images
                img_for_psf = img.replace('.fits', '.zogypsf')
                if option.fwhm:
                    fwhm0 = option.fwhm
                else:
                    fwhm0 = 5
                make_sn2 = False
                psfstars = 1
                catalog = lsc.__path__[0] + '/standard/cat/zero.cat'
            else:
                img_for_psf = img.replace('.fits', '')
                fwhm0 = option.fwhm
                make_sn2 = True
                psfstars = option.psfstars
                if option.use_sextractor:
                    catalog = ''
                elif option.catalog:
                    catalog = option.catalog
                else:
                    for system in ['sloan', 'apass', 'landolt']:
                        catalog = lsc.util.getcatalog(img, system)
                        if catalog:
                            break
            while True:
                result, fwhm = lsc.lscpsfdef.ecpsf(img_for_psf, fwhm0, option.threshold, psfstars,
                                                   option.distance, option.interactive, psffun, fixaperture,
                                                   catalog, option.datamin, option.datamax, option.show, make_sn2)
                print '\n### ' + str(result)
                if option.show:
#                    lsc.util.marksn2(img + '.fits', img + '.sn2.fits', 1, '')
                    iraf.delete('tmp.psf.fit?', verify=False)
                    iraf.seepsf(img.replace('.fits', '.psf'), '_psf.psf')
                    iraf.surface('_psf.psf')
                    aa = raw_input('>>>good psf [[y]/n] ? ')
                    if not aa or aa.lower()[0] == 'y':
                        break
                    if aa.lower()[0] == 'n':
                        result = 0
                        bb = raw_input('If you want to try again, type the new FWHM to try here. Otherwise press enter to continue. ')
                        if bb: fwhm0 = float(bb)
                        else: break
                else: 
                    break

            iraf.delete("tmp.*", verify="no")
            iraf.delete("_psf.*", verify="no")
            print  "********** Completed in ", int(time.time() - start_time), "sec"
            print result
            try:
                basename = os.path.basename(img)
                if result == 1:
                    lsc.mysqldef.updatevalue('photlco', 'psf', basename.replace('.fits', '.psf.fits'), basename)
                    lsc.mysqldef.updatevalue('photlco', 'fwhm', fwhm, basename)
                    lsc.mysqldef.updatevalue('photlco', 'mag', 9999, basename)
                    lsc.mysqldef.updatevalue('photlco', 'psfmag', 9999, basename)
                    lsc.mysqldef.updatevalue('photlco', 'apmag', 9999, basename)
                    if 'diff' not in img and os.path.isfile(img + '.diff.fits') and os.path.isfile(img + '.sn2.fits'):
                        # update diff info is file available
                        os.system('cp ' + img + '.sn2.fits ' + img + '.diff.sn2.fits')
                        lsc.mysqldef.updatevalue('photlco', 'psf', basename + '.psf.fits', basename + '.diff.fits')
                        lsc.mysqldef.updatevalue('photlco', 'fwhm', fwhm, basename + '.diff.fits')
                        lsc.mysqldef.updatevalue('photlco', 'mag', 9999, basename + '.diff.fits')
                        lsc.mysqldef.updatevalue('photlco', 'psfmag', 9999, basename + '.diff.fits')
                        lsc.mysqldef.updatevalue('photlco', 'apmag', 9999, basename + '.diff.fits')
                else:
                    print 'psf not good, database not updated '
                    lsc.mysqldef.updatevalue('photlco', 'psf', 'X', basename + '.fits')
                    if os.path.isfile(img + '.diff.fits'):
                        lsc.mysqldef.updatevalue('photlco', 'psf', 'X', basename + '.diff.fits')
            except:
                print 'module mysqldef not found'

