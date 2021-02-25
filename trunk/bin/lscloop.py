#!/usr/bin/env python
description = "> process lsc data  "

import re
import numpy as np
from argparse import ArgumentParser
from datetime import datetime, timedelta
import lsc
from multiprocessing import Pool, cpu_count
import os
import warnings

def multi_run_cosmic(args):
    return lsc.myloopdef.run_cosmic(*args)

# ###########################################################################

if __name__ == "__main__":   # main program
    parser = ArgumentParser(description=description)
    parser.add_argument("-e", "--epoch", help='args.epoch to reduce')
    parser.add_argument("-T", "--telescope", default='all')
    parser.add_argument("-I", "--instrument", default='', help='kb, fl, fs, sinistro, sbig, ep, muscat')
    parser.add_argument("-n", "--name", default='', help='object name')
    parser.add_argument("-d", "--id", default='')
    parser.add_argument("-f", "--filter", default='', nargs='+',
                        choices=['landolt', 'sloan', 'apass', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I', 'w'])
    parser.add_argument("-F", "--force", action="store_true")
    parser.add_argument("-b", "--bad", default='', help='bad stage',
                        choices=['wcs', 'psf', 'psfmag', 'zcat', 'abscat', 'mag', 'goodcat', 'quality', 'cosmic', 'diff'])
    parser.add_argument("-s", "--stage", default='',
                        choices=['wcs', 'psf', 'psfmag', 'zcat', 'abscat', 'mag', 'local', 'getmag', 'merge', 'mergeall', 'diff',
                                 'template', 'makestamp', 'apmag', 'cosmic', 'ingestsloan', 'ingestps1', 'fpack',
                                 'checkwcs', 'checkpsf', 'checkmag', 'checkquality', 'checkpos', 'checkcat',
                                 'checkmissing', 'checkfvd', 'checkcosmic', 'checkdiff'])
    parser.add_argument("-R", "--RA", default='')
    parser.add_argument("-D", "--DEC", default='')
    parser.add_argument("--RAS", default='')
    parser.add_argument("--DECS", default='')
    parser.add_argument("--RA0", default='')
    parser.add_argument("--DEC0", default='')
    parser.add_argument("-x", "--xord", default=3, type=int, help='x-order for bg fit')
    parser.add_argument("-y", "--yord", default=3, type=int, help='y-order for bg fit')
    parser.add_argument("--bkg", default=4., type=float, help='bkg radius for the fit')
    parser.add_argument("--size", default=7., type=float, help='size of the stamp for the fit')
    parser.add_argument("-t", "--threshold", default=5., type=float, help='Source detection threshold')
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("--show", action="store_true")
    parser.add_argument("-c", "--center", action="store_false", dest='recenter', help='do not recenter')
    parser.add_argument("--unfix", action="store_false", dest='fix', help='use a variable color term')
    parser.add_argument("--cutmag", default=99., type=float, help='magnitude instrumental cut for zeropoint')
    parser.add_argument("--field", default='', choices=['landolt', 'sloan', 'apass', 'gaia'],
            help='note: gaia only functional for psf stage')
    parser.add_argument("--ref", default='', help='get sn position from this file')
    parser.add_argument("--use-sextractor", action="store_true", help="use souces from sextractor for PSF instead of catalog")
    parser.add_argument("--catalogue", default='', help="filename of catalog (full path OR ./___apass.cat OR apass/___apass.cat)")
    parser.add_argument("--calib", default='', choices=['sloan', 'natural', 'sloanprime'])
    parser.add_argument("--sigma-clip", default=2., help='number of sigma at which to reject stars for zero point calibration')
    parser.add_argument("--type", choices=['fit', 'ph', 'mag'], default='', help='type of magnitude (PSF, aperture, apparent); default ph for filetype=3 and fit for everything else')
    parser.add_argument("--standard", default='', help='use the zeropoint from this standard')
    parser.add_argument("--xshift", default=0, type=int, help='x-shift in the guess astrometry')
    parser.add_argument("--yshift", default=0, type=int, help='y-shift in the guess astrometry')
    parser.add_argument("--fwhm", default='', help='fwhm (in pixel)')
    parser.add_argument("--mode", choices=['sv', 'astrometry'], default='sv', help='mode for wcs')
    parser.add_argument("--combine", default=1e-10, type=float, help='range to combine (in days)')
    parser.add_argument("--datamax", type=int, help='data max for saturation (counts)')
    parser.add_argument("--datamin", type=int, help='data min for saturation (counts)')
    parser.add_argument("--filetype", choices=range(5), default=1, type=int, help='filetype 1 [single], 2 [merge], 3 [difference], 0 [all]')
    parser.add_argument("-o", "--output", default='', help='output file')
    parser.add_argument("--tempdate", default='19990101-20080101', help='template date')
    parser.add_argument("--temptel", default='', help='--temptel  template instrument')
    parser.add_argument("--normalize", choices=['i', 't'], default='t', help='normalize to image [i], template [t] (hotpants parameter)')
    parser.add_argument("--convolve", choices=['i', 't'], default='', help='force convolution with image [i], template [t] (hotpants parameter)')
    parser.add_argument("--z1", type=float)
    parser.add_argument("--z2", type=float)
    parser.add_argument("--groupidcode", type=int)
    parser.add_argument("--targetid", default=None,   type=int)
    parser.add_argument("--ps1frames", default='', help='list of ps1 frames (download them manually)')
    parser.add_argument('--zcatold', action='store_true', help='use original zero point and color term routine')
    parser.add_argument("--bgo", default=3., type=float, help=' bgo parameter for hotpants')
    parser.add_argument("-p", "--psf", default='', help='psf image for template')
    parser.add_argument("--mag", type=float, default=0., help='mag to subtract from template image')
    parser.add_argument("--uncleaned", action='store_false', dest='clean', help='do not use cosmic ray cleaned image as template')
    parser.add_argument("--subtract-mag-from-header", action='store_true', help='automatically subtract mag from header of template image')
    parser.add_argument("--fixpix", action="store_true", help='Run fixpix on the images before doing image subtraction')
    parser.add_argument("--nstars", type=int, default=6, help="number of stars used to make the PSF")
    parser.add_argument("--minstars", type=int, default=0, help="minimum number of stars matching catalog (-s abscat/local)")
    parser.add_argument("--difftype", type=int, choices=[0, 1], help='Choose hotpants (0) or optimal (1) subtraction')
    parser.add_argument("--multicore", default=8, type=int, help='numbers of cores')
    parser.add_argument("--unmask", action='store_false', dest='use_mask', help='do not use mask for PyZOGY gain calculation')
    parser.add_argument("--banzai", action="store_true", help="Use souces from BANZAI catalogs stored in image headers for PSF")
    parser.add_argument("--b_sigma", type=float, default=3.0, help="value used to sigma-clip BANZAI sources")
    parser.add_argument("--b_crlim", type=float, default=3.0, help="lower limit used to reject CRs identified as BANZAI sources")
    parser.add_argument("--no_iraf", action="store_true", help="Don't use iraf (currently only option in checkpsf stage)")
    parser.add_argument("--max_apercorr", type=float, default=0.1, help="absolute value of the maximum aperture correction (mag) above which a PSF is marked as bad")

    args = parser.parse_args()
    if args.multicore >= cpu_count():
        args.multicore = cpu_count()-1 if cpu_count() > 1 else 1
        print('Attempting to run on too many cores, reducing to {n}'.format(n=args.multicore))

    if args.stage == 'checkdiff':
        filetype = 3
    elif args.bad == 'diff':
        filetype = 1
    elif args.stage == 'fpack':
        filetype = 0
    else:
        filetype = args.filetype
    
    #Set aperture photometry as default for difference imaging
    if filetype == 3 and not args.type:
        args.type = 'ph'
    elif not args.type:
        args.type = 'fit'

    filters = ','.join(args.filter)

    ll = lsc.myloopdef.get_list(args.epoch, args.telescope, filters, args.bad, args.name, args.id, args.RA, args.DEC,
                                'photlco', filetype, args.groupidcode, args.instrument, args.temptel, args.difftype, None, args.targetid)
    if ll:
        if args.stage != 'merge':
            print '##' * 50
            print '# IMAGE                                    OBJECT           FILTER           WCS            ' \
                  ' PSF           PSFMAG    APMAG       ZCAT          MAG      ABSCAT'
            for i in range(0, len(ll['filename'])):
                try:
                    print '%s\t%12s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s' % \
                          (ll['filename'][i].replace('.fits', ''), ll['objname'][i], ll['filter'][i],
                           str(ll['wcs'][i]), ll['psf'][i].replace('.fits', ''),
                           str(round(ll['psfmag'][i], 4)), str(round(ll['apmag'][i], 4)), ll['zcat'][i].replace('.cat', ''),
                           str(round(ll['mag'][i], 4)), ll['abscat'][i].replace('.cat', ''))
                except:
                    print '%s\t%12s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s' % \
                          (ll['filename'][i].replace('.fits', ''), ll['objname'][i], ll['filter'][i],
                           str(ll['wcs'][i]), ll['psf'][i],
                           str(ll['psfmag'][i]), str(ll['apmag'][i]), ll['zcat'][i].replace('.cat', ''),
                           str(ll['mag'][i]), ll['abscat'][i])
            print '\n###  total number = ' + str(len(ll['filename']))
            if args.standard:
                if args.standard == 'all':
                    mm = lsc.myloopdef.get_standards(args.epoch, args.name, filters)
                else:
                    mm = lsc.myloopdef.get_list(args.epoch, args.telescope, filters, _instrument=args.instrument, _name=args.standard)
                for i in range(len(mm['filename'])):
                    print '%s\t%12s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s' % \
                          (str(mm['filename'][i]), str(mm['objname'][i]), str(mm['filter'][i]),
                           str(mm['wcs'][i]), str(mm['psf'][i]),
                           str(mm['psfmag'][i]), str(mm['zcat'][i]), str(mm['mag'][i]),
                           str(mm['abscat'][i]))
                print '\n###  total number = ' + str(len(mm['filename']))
                if args.stage not in ['mag', 'abscat', 'local']:
                    ll = mm
            else:
                mm = {'filename': []}
            listfile = np.array([k + v for k, v in zip(ll['filepath'], ll['filename'])])
            # ####################################
            if args.stage not in ['', 'getmag', 'fpack']: # unpack any fpacked files before anything else
                packed_files = listfile[ll['lastunpacked'] == np.array(None)]
                def funpack(filename):
                    return os.system('funpack -v -D ' + filename)
                p = Pool(args.multicore)
                exitcodes = p.map(funpack, packed_files)
                p.close()
                p.join()
                unpacked = [os.path.basename(filepath) for filepath, exitcode in zip(packed_files, exitcodes) if exitcode == 0]
                query = 'update photlco set lastunpacked="{}" where filename="'.format(datetime.utcnow())
                query += '" or filename="'.join(unpacked) + '"'
                lsc.mysqldef.query([query], lsc.conn)

            if args.stage == 'fpack':
                is_unpacked = ll['lastunpacked'] != np.array(None)
                unpacked_files = listfile[is_unpacked]
                if not args.force: # only fpack month-old data unless run with -F
                    month_old = np.array(ll['lastunpacked'])[is_unpacked] < datetime.utcnow() - timedelta(31)
                    unpacked_files = unpacked_files[month_old]
                def fpack(filename):
                    return os.system('fpack -q 64 -v -D -Y ' + filename)
                p = Pool(args.multicore)
                exitcodes = p.map(fpack, unpacked_files)
                p.close()
                p.join()
                packed = [os.path.basename(filepath) for filepath, exitcode in zip(unpacked_files, exitcodes) if exitcode == 0]
                query = 'update photlco set lastunpacked=NULL where filename="' + '" or filename="'.join(packed) + '"'
                lsc.mysqldef.query([query], lsc.conn)
            elif args.stage == 'getmag':  # get final magnitude from mysql
                lsc.myloopdef.run_getmag(ll['filename'], args.output, args.interactive, args.show, args.combine, args.type)
            elif args.stage == 'psf':
                catalogue = lsc.util.getcatalog(args.name, args.field) if args.field else args.catalogue
                lsc.myloopdef.run_psf(ll['filename'], args.threshold, args.interactive, args.fwhm, args.show, args.force, args.fix, catalogue,
                                      'photlco', args.use_sextractor, args.datamin, args.datamax, args.nstars, args.banzai, args.b_sigma, args.b_crlim, 
                                      max_apercorr=args.max_apercorr)
            elif args.stage == 'psfmag':
                lsc.myloopdef.run_fit(ll['filename'], args.RAS, args.DECS, args.xord, args.yord, args.bkg, args.size, args.recenter, args.ref,
                                      args.interactive, args.show, args.force, args.datamax,args.datamin,'photlco',args.RA0,args.DEC0)
            elif args.stage == 'wcs':
                lsc.myloopdef.run_wcs(ll['filename'], args.interactive, args.force, args.xshift, args.xshift, args.catalogue,'photlco',args.mode)
            elif args.stage == 'makestamp':
                lsc.myloopdef.makestamp(ll['filename'], 'photlco', args.z1, args.z2, args.interactive, args.force, args.output)
            elif args.stage == 'apmag':
                lsc.myloopdef.run_apmag(ll['filename'], 'photlco')
            elif args.stage == 'cosmic':
                if args.multicore > 1:
                    p = Pool(args.multicore)
                    inp = [([i], 'photlco', 4.5, 0.2, 4,args.force) for i in listfile]
                    p.map(multi_run_cosmic, inp)
                    p.close()
                    p.join()
                else:
                    lsc.myloopdef.run_cosmic(listfile, 'photlco', 4.5, 0.2, 4, args.force)
            elif args.stage == 'ingestsloan':
                lsc.myloopdef.run_ingestsloan(listfile, 'sloan', show=args.show, force=args.force)
            elif args.stage == 'ingestps1':
                lsc.myloopdef.run_ingestsloan(listfile, 'ps1', args.ps1frames, show=args.show, force=args.force)
            elif args.stage == 'zcat':
                def run_absphot(img):
                    return lsc.lscabsphotdef.absphot(img, args.field, args.catalogue, args.fix, args.sigma_clip, args.interactive,
                                                     args.type, args.force, args.show, args.cutmag, args.calib, args.zcatold)
                args.multicore = 1 # parallel processing doesn't work yet; problem with too many mysql queries
                if args.multicore > 1:
                    p = Pool(args.multicore)
                    p.map(run_absphot, listfile)
                    p.close()
                    p.join()
                else:                    
                    for img in listfile:
                        if args.field == 'gaia':
                            print('Cannot use gaia catalog for zcat stage: ' 
                                  'use apass or sloan instead')
                            continue
                        run_absphot(img)
            elif args.stage in ['mag', 'abscat', 'local']:  # compute magnitudes for sequence stars or supernova
                # Don't allow PSF photometry to be performed on difference images becaues of complications with aperture correction
                if args.filetype == 3 and args.stage == 'mag' and args.type =='fit':
                    warnings.warn('Aperture photometry recommended for difference images (--type ph)', UserWarning)
                if args.catalogue:
                    catalogue = args.catalogue
                elif args.field:
                    catalogue = lsc.util.getcatalog(args.name, args.field)
                else:
                    catalogue = lsc.util.getcatalog(args.name, 'apass')
                if not args.field and args.filter and args.filter[0] in ['landolt', 'sloan', 'apass']:
                    field = args.filter[0]
                else:
                    field = args.field
                lsc.myloopdef.run_cat(ll['filename'], mm['filename'], args.interactive, args.stage, args.type, 'photlco', field, catalogue, args.force, args.minstars)
            elif args.stage == 'diff':  #    difference images using hotpants
                if not args.name and not args.targetid:
                    raise Exception('you need to select one object: use option -n/--name')
                if args.telescope=='all':
                    raise Exception('you need to select one type of instrument -T [fs, fl, kb, ep]')
                
                startdate = args.tempdate.split('-')[0]
                enddate   = args.tempdate.split('-')[-1]

                if args.difftype == 1:
                    suffix = '.optimal.{}.diff.fits'.format(args.temptel).replace('..', '.')
                else:
                    suffix = '.{}.diff.fits'.format(args.temptel).replace('..', '.')

                if args.temptel.upper() in ['SDSS', 'PS1']:
                    if args.telescope == 'kb':
                        fake_temptel = 'sbig'
                    elif args.telescope == 'fs':
                        fake_temptel = 'spectral'
                    elif args.telescope == 'fl':
                        fake_temptel = 'sinistro'
                    elif args.telescope == 'fa':
                        fake_temptel = 'sinistro'
                    elif args.telescope == 'ep':
                        fake_temptel = 'muscat'
                elif args.temptel:
                    fake_temptel = args.temptel
                else:
                    fake_temptel = args.telescope

                lista = lsc.mysqldef.getlistfromraw(lsc.myloopdef.conn, 'photlco', 'dayobs', startdate, enddate, '*', fake_temptel)
                if lista:
                    ll00 = {}
                    for jj in lista[0].keys():
                        ll00[jj] = []
                    for i in range(0, len(lista)):
                        for jj in lista[0].keys():
                            ll00[jj].append(lista[i][jj])
                    inds = np.argsort(ll00['mjd'])  #  sort by mjd
                    for i in ll00.keys():
                        ll00[i] = np.take(ll00[i], inds)
                    lltemp = lsc.myloopdef.filtralist(ll00, filters, '', args.name, args.RA, args.DEC, '', 4, args.groupidcode, '', '', '', None, args.targetid)

                if not lista or not lltemp:
                    raise Exception('template not found')

                listtemp = np.array([k + v for k, v in zip(lltemp['filepath'], lltemp['filename'])])

                lsc.myloopdef.run_diff(listfile, listtemp, args.show, args.force, args.normalize, args.convolve, args.bgo, args.fixpix, args.difftype, suffix, args.use_mask, args.no_iraf)
            elif args.stage == 'template':
                lsc.myloopdef.run_template(listfile, args.show, args.force, args.interactive, args.RA, args.DEC, args.psf, args.mag, args.clean, args.subtract_mag_from_header)
            elif args.stage == 'mergeall':  #    merge images using lacos and swarp
                lsc.myloopdef.run_merge(listfile, args.force)
            elif args.stage == 'checkpsf':
                lsc.myloopdef.checkpsf(ll['filename'],args.no_iraf)
            elif args.stage == 'checkmag':
                lsc.myloopdef.checkmag(ll['filename'], args.datamax)
            elif args.stage == 'checkwcs':
                lsc.myloopdef.checkwcs(ll['filename'], args.force, 'photlco', args.z1, args.z2)
            elif args.stage == 'checkfast':
                lsc.myloopdef.checkfast(ll['filename'], args.force)
            elif args.stage == 'checkquality':
                lsc.myloopdef.checkquality(ll['filename'])
            elif args.stage == 'checkpos':
                lsc.myloopdef.checkpos(ll['filename'], args.RA, args.DEC)
            elif args.stage == 'checkcat':
                lsc.myloopdef.checkcat(ll['filename'])
            elif args.stage == 'checkmissing':
                lsc.myloopdef.check_missing(ll['filename'])
            elif args.stage == 'checkfvd':
                lsc.myloopdef.checkfilevsdatabase(ll)
            elif args.stage == 'checkcosmic':
                lsc.myloopdef.checkcosmic(ll['filename'])
            elif args.stage == 'checkdiff':
                lsc.myloopdef.checkdiff(ll['filename'])
            elif args.stage:
                print args.stage + ' not defined'
    # ################################################
        else: # if args.stage == 'merge'
            for epo in np.unique(ll['dayobs']):
                print '\n#### ' + str(epo)
                ll0 = {key: val[ll['dayobs'] == epo] for key, val in ll.items()}
                for i in range(0, len(ll0['filename'])):
                    print '%s\t%12s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s' % \
                          (str(ll0['filename'][i]), str(ll0['objname'][i]), str(ll0['filter'][i]), str(ll0['wcs'][i]),
                           str(ll0['psf'][i]),
                           str(ll0['psfmag'][i]), str(ll0['zcat'][i]), str(ll0['mag'][i]), str(ll0['abscat'][i]))
                print '\n###  total number = ' + str(len(ll0['filename']))
                lsc.myloopdef.run_merge(listfile[ll['dayobs'] == epo], args.force) # merge images using lacos and swarp
    else:
        print '\n### no data selected'
