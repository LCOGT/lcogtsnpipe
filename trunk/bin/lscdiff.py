#!/usr/bin/env python
description = ">> make different image using hotpants"
usage = "%prog imagein  imagetem [options] "
import os
import lsc
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from argparse import ArgumentParser
from PyZOGY.subtract import run_subtraction
from pyraf import iraf
from reproject import reproject_interp
import shutil
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ZScaleInterval

def crossmatchtwofiles(img1, img2, radius=3):
    ''' This module is crossmatch two images:
        It run sextractor transform the pixels position of the the sources in coordinates and crossmatch them  
        The output is a dictionary with the objects in common
    '''

    hd1 = fits.getheader(img1)
    hd2 = fits.getheader(img2)
    wcs1 = WCS(hd1)
    wcs2 = WCS(hd2)

    xpix1, ypix1, fw1, cl1, cm1, ell1, bkg1, fl1 = lsc.lscastrodef.sextractor(img1)
    xpix2, ypix2, fw2, cl2, cm2, ell2, bkg2, fl2 = lsc.lscastrodef.sextractor(img2)
    xpix1, ypix1, xpix2, ypix2 = np.array(xpix1, float), np.array(ypix1, float), np.array(xpix2, float), np.array(ypix2, float)

    bb = wcs1.wcs_pix2world(zip(xpix1, ypix1), 1)  #   transform pixel in coordinate
    xra1, xdec1 = zip(*bb)
    bb = wcs2.wcs_pix2world(zip(xpix2, ypix2), 1)  #   transform pixel in coordinate
    xra2, xdec2 = zip(*bb)

    xra1, xdec1, xra2, xdec2 = np.array(xra1, float), np.array(xdec1, float), np.array(xra2, float), np.array(xdec2, float)
    distvec, pos1, pos2 = lsc.lscastrodef.crossmatch(xra1, xdec1, xra2, xdec2, radius)
    #dict={}
    dict = {'ra1': xra1[pos1], 'dec1': xdec1[pos1], 'ra2': xra2[pos2], 'dec2': xdec2[pos2], \
            'xpix1': xpix1[pos1], 'ypix1': ypix1[pos1], 'xpix2': xpix2[pos2], 'ypix2': ypix2[pos2]}
    np.savetxt('substamplist', zip(xpix1[pos1], ypix1[pos1]), fmt='%10.10s\t%10.10s')
    return 'substamplist', dict

###################################################
if __name__ == "__main__":
    parser = ArgumentParser(usage=usage, description=description)
    parser.add_argument("targlist", help='list of images with the target in them')
    parser.add_argument("templist", help='list of template images')
    parser.add_argument("-f", "--force", dest="force", action="store_true",
                      default=False, help=' force archiving \t\t\t [%default]')
    parser.add_argument("--show", dest="show", action="store_true",
                      default=False, help=' show result  \t\t\t [%default]')
    parser.add_argument("--fixpix", dest="fixpix", action="store_true", default=False,
                      help='Run fixpix on the images before doing the subtraction')
    parser.add_argument('--suffix', default='.diff.fits', help='suffix for difference images')

    hotpants = parser.add_argument_group("hotpants parameters")
    hotpants.add_argument("--nrxy", dest="nrxy", default='1,1',
                        help='Number of image region in x y directions \t [%default]')
    hotpants.add_argument("--nsxy", dest="nsxy", default='8,8',
                        help="Number of region's stamps in x y directions\t [%default]")
    hotpants.add_argument("--ko", dest="ko", default='2',
                        help='spatial order of kernel variation within region\t [%default]')
    hotpants.add_argument("--bgo", dest="bgo", default='2',
                        help='spatial order of background variation within region \t [%default]')
    hotpants.add_argument("--afssc", dest="afssc", default=False,
                        action="store_true", help='use selected stamps \t\t\t [%default]')
    hotpants.add_argument("--normalize", dest="normalize", default='i', choices=['i', 't'],
                        help='normalize zero point to image [i] or template [t] \t [%default]')
    hotpants.add_argument("--convolve", dest="convolve", default='', choices=['i', 't', ''],
                        help='convolve direction to image [i] or template [t] \t [%default]')
    hotpants.add_argument("--interpolation", dest="interpolation", default='drizzle',
                          choices=['drizzle', 'nearest', 'linear', 'poly3', 'poly5', 'spline3'],
                        help='interpolation algorithm  [drizzle,nearest,linear,poly3,poly5,spline3]\t [%default]')
    parser.add_argument("--difftype", type=str, default='0', help='Choose hotpants (0) or optimal (1) subtraction \t [%(default)s]')
    parser.add_argument("--unmask", action='store_false', dest='use_mask', help='do not use mask for PyZOGY gain calculation')
    parser.add_argument("--no-iraf", action='store_true', help='transform images in Python instead of IRAF'
                                                               '(IRAF is still used for fixpix if run with --fixpix)')

    args = parser.parse_args()
    imglisttar = lsc.util.readlist(args.targlist)
    imglisttemp = lsc.util.readlist(args.templist)

    listatar = {}
    try:
        hostname, username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
        conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
    except:
        print not 'connected to the database'
        conn=0
    # divide targets by filter and targetid
    for img in imglisttar:
        hdr = lsc.util.readhdr(img)
        try:
            _targetid = lsc.mysqldef.targimg(img)
        except:
            _targetid = 1

        _filt = lsc.util.readkey3(hdr, 'filter')
        _filter = lsc.sites.filterst1[_filt]

        _obj = lsc.util.readkey3(hdr, 'object')
        if _filter not in listatar:
            listatar[_filter] = {}
        if _targetid not in listatar[_filter]:
            listatar[_filter][_targetid] = []
        listatar[_filter][_targetid].append(img)

    # divide template by filter and targetid
    listatemp = {}
    for img in imglisttemp:
        hdr = lsc.util.readhdr(img)
        try:
            _targetid = lsc.mysqldef.targimg(img)
        except:
            _targetid = 1
        _filt = lsc.util.readkey3(hdr, 'filter')
        _filter = lsc.sites.filterst1[_filt]

        _obj = lsc.util.readkey3(hdr, 'object')
        if _filter not in listatemp:
            listatemp[_filter] = {}
        if _targetid not in listatemp[_filter]:
            listatemp[_filter][_targetid] = []
        listatemp[_filter][_targetid].append(img)

    for f in listatar:
        for o in listatar[f]:
            if f in listatemp:
                if o in listatemp[f]:
                    imglist1 = listatar[f][o]
                    imglist2 = listatemp[f][o]
                    for imgtarg_path in imglist1:
                        _dir, imgtarg0 = os.path.split(imgtarg_path)
                        if _dir: _dir += '/'
                        imgtemp_path = imglist2[0]
                        _dirtemp, imgtemp0 = os.path.split(imgtemp_path)
                        if _dirtemp: _dirtemp += '/'

                        imgout0 = imgtarg0.replace('.fits', args.suffix)

                        if os.path.isfile(_dir + imgout0) and not args.force:
                            print 'file', imgout0, 'already there'
                            continue
                        targmask0 = imgtarg0.replace('.fits', '.mask.fits')
                        if not os.path.isfile(_dir + targmask0):
                            print "no cosmic ray mask for target image, run 'lscloop.py -s cosmic' first"
                            continue
                        tempmask0 = imgtemp0.replace('.fits', '.mask.fits')
                        if not os.path.isfile(_dirtemp + tempmask0):
                            print "no cosmic ray mask for template image, run 'lscloop.py -s cosmic' first"
                            continue
                        tempnoise0 = imgtemp0.replace('.fits', '.var.fits')
                        outmask0 = imgout0.replace('.fits', '.mask.fits')
                        artar, hdtar = fits.getdata(_dir+imgtarg0, header=True)

                        if os.path.isfile(_dirtemp+imgtemp0.replace('.fits', '.sn2.fits')):
                            hdtempsn = fits.getheader(_dirtemp+imgtemp0.replace('.fits', '.sn2.fits'))
                        else:
                            hdtempsn = {}

                        imgtemp = '_temp.fits'
                        imgtarg = '_targ.fits'
                        imgout = '_out.fits'
                        tempmask = '_tempmask.fits'
                        targmask = '_targmask.fits'
                        tempnoise = 'tempnoise.fits'

                        if args.no_iraf:
                            fits.writeto(imgtarg, artar, hdtar, overwrite=True)
                            shutil.copy(_dir + targmask0, targmask)

                            temp_reproj, temp_foot = reproject_interp(_dirtemp + imgtemp0, hdtar)
                            fits.writeto(imgtarg, temp_reproj, hdtar)

                            tempmask_reproj, tempmask_foot = reproject_interp(_dirtemp + tempmask0, hdtar)
                            fits.writeto(tempmask, tempmask_reproj, hdtar)

                            tempnoise_reproj, tempnoise_foot = reproject_interp(_dirtemp + tempnoise0, hdtar)
                            fits.writeto(tempnoise, tempnoise_reproj, hdtar)

                            if args.show:
                                fig, (ax1, ax2) = plt.subplots(1, 2)
                                norm1 = ImageNormalize(artar, interval=ZScaleInterval())
                                ax1.imshow(artar, norm=norm1, origin='lower')
                                norm2 = ImageNormalize(temp_reproj, interval=ZScaleInterval())
                                ax2.imshow(temp_reproj, norm=norm2, origin='lower')

                        else:
                            substamplist, dict = crossmatchtwofiles(_dir + imgtarg0, _dirtemp + imgtemp0, 4)
                            xra1, xdec1, xra2, xdec2, xpix1, ypix1, xpix2, ypix2 = dict['ra1'], dict['dec1'], dict['ra2'], \
                                                                                   dict['dec2'], dict['xpix1'], \
                                                                                   dict['ypix1'], dict['xpix2'], \
                                                                                   dict['ypix2']

                            vector4 = [str(k) + ' ' + str(v) + ' ' + str(j) + ' ' + str(l) for k, v, j, l in
                                       zip(xpix1, ypix1, xpix2, ypix2)]
                            if len(vector4) >= 12:
                                num = 3
                            else:
                                num = 2
                            np.savetxt('tmpcoo', vector4, fmt='%1s')
                            iraf.immatch.geomap('tmpcoo', "tmp$db", 1, hdtar['NAXIS1'], 1, hdtar['NAXIS2'],
                                                fitgeom="general", functio="legendre", xxor=num, xyor=num, xxterms="half",
                                                yxor=num, yyor=num, yxterms="half", calctype="real", inter='No')

                            lsc.util.delete(imgtemp)
                            lsc.util.delete(imgtarg)
                            lsc.util.delete(imgout)
                            lsc.util.delete(tempmask)
                            lsc.util.delete(targmask)
                            lsc.util.delete('tempmask.fits')
                            iraf.imcopy(_dir + imgtarg0 + '[0]', imgtarg, verbose='yes')
                            iraf.imcopy(_dir + targmask0, targmask, verbose='yes')
    #                        try:
                            iraf.immatch.gregister(_dirtemp + imgtemp0, imgtemp, "tmp$db", "tmpcoo", geometr="geometric",
                                                   interpo=args.interpolation, boundar='constant', constan=0, flux='yes', verbose='yes')
                            print 'here'
                            try:
                                iraf.immatch.gregister(_dirtemp + tempmask0, tempmask, "tmp$db", "tmpcoo", geometr="geometric",
                                                   interpo=args.interpolation, boundar='constant', constan=0, flux='yes', verbose='yes')
                            except:
                                print 'try again'
                                # this is strange, sometime the registration of the msk fail the first time, but not the second time
                                iraf.immatch.gregister(_dirtemp + tempmask0, tempmask, "tmp$db", "tmpcoo", geometr="geometric",
                                                   interpo=args.interpolation, boundar='constant', constan=0, flux='yes', verbose='yes')

                            if os.path.isfile(_dirtemp + imgtemp0.replace('.fits','.var.fits')):
                                print 'variance image already there, do not create noise image'
                                iraf.immatch.gregister(_dirtemp + tempnoise0, tempnoise, "tmp$db", "tmpcoo", geometr="geometric",
                                                   interpo=args.interpolation, boundar='constant', constan=0, flux='yes', verbose='yes')

                            if args.show:
                                iraf.set(stdimage='imt2048')
                                iraf.display(_dir + imgtarg0 + '[0]', frame=1, fill='yes')
                                iraf.display(imgtarg, frame=2, fill='yes')
                                iraf.display(_dirtemp + imgtemp0 + '[0]', frame=3, fill='yes')
                                iraf.display(imgtemp, frame=4, fill='yes')

                        data_targ, head_targ = fits.getdata(imgtarg, header=True)
                        exp_targ  = lsc.util.readkey3(head_targ, 'exptime')
                        sat_targ = lsc.util.readkey3(head_targ, 'datamax')
                        gain_targ = lsc.util.readkey3(head_targ, 'gain')
                        rn_targ = lsc.readkey3(head_targ,'ron')

                        data_temp, head_temp = fits.getdata(imgtemp, header=True)
                        exp_temp = lsc.util.readkey3(head_temp, 'exptime')
                        sat_temp = lsc.util.readkey3(head_temp, 'datamax')
                        gain_temp = lsc.util.readkey3(head_temp, 'gain')

                        if 'rdnoise' in head_temp:
                            rn_temp   = head_temp['rdnoise']
                        else:
                            rn_temp   = 1

                        if 'L1FWHM' in head_targ: 
                            max_fwhm = head_targ['L1FWHM']  # to be check
                        else:                     
                            max_fwhm = 3.0

                        targmask_data = fits.getdata(targmask)

                        # round all values in template mask up to 1 (unless they are 0)
                        tempmask_data, tempmask_header = fits.getdata(tempmask, header=True)
                        tempmask_int = (tempmask_data > 0).astype('uint8')
                        tempmask_fits = fits.PrimaryHDU(header=tempmask_header, data=tempmask_int)
                        tempmask_fits.writeto(tempmask, overwrite=True, output_verify='fix')

                        # create the noise images
                        median = np.median(data_targ)
                        noise = 1.4826*np.median(np.abs(data_targ - median))
                        pssl_targ = gain_targ*noise**2 - rn_targ**2/gain_targ - median
                        #noiseimg = (data_targ - median)**2
                        noiseimg = data_targ + pssl_targ + rn_targ**2
                        noiseimg[targmask_data > 0] = sat_targ
                        fits.writeto('targnoise.fits', noiseimg, output_verify='fix', overwrite=True)

#                        print 'variance image already there, do not create noise image'
                        if not os.path.isfile(_dirtemp + tempnoise0):
                            median = np.median(data_temp)
                            noise = 1.4826*np.median(np.abs(data_temp - median))
                            pssl_temp = gain_temp*noise**2 - rn_temp**2/gain_temp - median
                            #noiseimg = (data_temp - median)**2
                            noiseimg = data_temp + pssl_temp + rn_temp**2
                            noiseimg[tempmask_data > 0] = sat_temp
                            fits.writeto(tempnoise, noiseimg, output_verify='fix', overwrite=True)
                        else:
                            pssl_temp = 0
                            print 'variance image already there, do not create noise image'

                        #  if skylevel is in the header, swarp with bkg subtraction has been applyed
                        if 'SKYLEVEL' in head_temp:
                            pssl_temp = head_temp['SKYLEVEL']

                        # create mask image for template
                        mask = np.abs(data_temp) < 1e-6
                        fits.writeto('tempmask.fits',mask.astype('i'))

                        if args.fixpix:
                            try:
                                iraf.flpr(); iraf.flpr()
                                iraf.unlearn(iraf.fixpix)
                                cwd = os.getcwd()
                                iraf.fixpix(os.path.join(cwd, imgtarg),
                                            os.path.join(cwd, targmask), verbose='yes')
                                iraf.flpr(); iraf.flpr()
                                iraf.unlearn(iraf.fixpix)
                                iraf.fixpix(os.path.join(cwd, imgtemp),
                                            os.path.join(cwd, tempmask), verbose='yes')
                                iraf.flpr(); iraf.flpr()
                                iraf.unlearn(iraf.fixpix)
                            except:
                                print 'FIXPIX ERROR, continuing without fixpix'
                        # hotpants parameters
                        iuthresh = str(sat_targ)                        # upper valid data count, image
                        iucthresh = str(0.95*sat_targ)                   # upper valid data count for kernel, image
                        tuthresh = str(sat_temp)                        # upper valid data count, template
                        tucthresh = str(0.95*sat_temp)                   # upper valid data count for kernel, template
                        rkernel = str(np.median([10, 2.*max_fwhm, 20])) # convolution kernel half width
                        radius = str(np.median([15, 3.0*max_fwhm, 25])) # HW substamp to extract around each centroid
                        sconv = '-sconv'                                # all regions convolved in same direction (0)

                        if args.convolve in ['t','i']:
                            _convolve=' -c ' + args.convolve + ' '
                        else:
                            _convolve=''

                        if args.afssc:
                            substamplist, xpix1, ypix1, xpix2, ypix2 = crossmatchtwofiles(imgtarg, imgtemp)
                            _afssc = ' -cmp ' + str(substamplist) + ' -afssc 1 '
                        else:
                            _afssc = ''

                        if args.difftype == '1':
                            psftarg = '_targpsf.fits'
                            psftemp = '_temppsf.fits'
                            if args.no_iraf:
                                lsc.myloopdef.seepsf(imgtarg_path.replace('.fits', '.psf.fits'), psftarg)
                                lsc.myloopdef.seepsf(imgtemp_path.replace('.fits', '.psf.fits'), psftemp)
                            else:
                                iraf.noao()
                                iraf.digiphot()
                                iraf.daophot(_doprint=0)
                                os.system('rm {0} {1}'.format(psftarg, psftemp))
                                iraf.seepsf(imgtarg_path.replace('.fits','.psf.fits'), psftarg)
                                iraf.seepsf(imgtemp_path.replace('.fits','.psf.fits'), psftemp)
                            try:
                                print 'Passing images to PyZOGY'
                                run_subtraction(imgtarg, imgtemp, psftarg, psftemp,
                                                           science_mask=targmask,
                                                           reference_mask=tempmask,
                                                           science_saturation=sat_targ,
                                                           reference_saturation=sat_temp,
                                                           n_stamps=1,
                                                           output=imgout,
                                                           normalization=args.normalize,
                                                           show=args.show,
                                                           use_mask_for_gain=args.use_mask)
                            except Exception as e:
                                print e
                                print 'PyZOGY failed on', imgtarg0
                                continue

                            # create fields in header that hotpants does
                            hotpants_fields = {'TARGET': (imgtarg_path, 'target image'),
                                               'TEMPLATE': (imgtemp_path, 'template image'),
                                               'DIFFIM': (imgout, 'Difference Image'),
                                               'NREGION': (1, 'Number of independent regions'),
                                               'MASKVAL': (1e-30, 'Value of Masked Pixels')}
                            lsc.util.updateheader(imgout, 0,  hotpants_fields)


                        else:
                            line = ('hotpants -inim ' + imgtarg + ' -tmplim ' + imgtemp + ' -outim ' + imgout +
                                    ' -tu ' + tuthresh + ' -tuk ' + tucthresh +
                                    ' -tl ' + str(min(-pssl_temp,0)) + ' -tg ' + str(gain_temp) +
                                    ' -tr ' + str(rn_temp) + ' -tp ' + str(min(-pssl_temp,0)) +
                                    ' -tni tempnoise.fits ' +
                                    ' -iu ' + str(iuthresh) + ' -iuk ' + str(iucthresh) + 
                                    ' -il ' + str(min(-pssl_targ,0)) + ' -ig ' + str(gain_targ) +
                                    ' -ir ' + str(rn_targ) + ' -ip ' + str(min(-pssl_targ,0)) +
                                    ' -ini targnoise.fits ' +
                                    ' -r ' + str(rkernel) +
                                    ' -nrx ' + args.nrxy.split(',')[0] + ' -nry ' + args.nrxy.split(',')[1] +
                                    ' -nsx ' + args.nsxy.split(',')[0] + ' -nsy ' + args.nsxy.split(',')[1] +
                                    _afssc + ' -rss ' + str(radius) +
                                    _convolve + ' -n ' + args.normalize + ' ' + sconv +
                                    ' -ko ' + args.ko + ' -bgo ' + args.bgo+' -tmi tempmask.fits  -okn -ng 4 9 0.70 6 1.50 4 3.00 2 5 ')


                            print line
                            os.system(line)


                        # delete temporary files
                        os.system('rm -f ' + imgout0.replace('.fits', '.conv.fits') +
                                  ' ' + imgout0.replace('.fits', 'xy') + 
                                  ' ' + imgout0.replace('.fits', 'xy.skipped') + 
                                  ' ' + imgout0.replace('.fits', 'xy.all'))

                        hd = fits.getheader(imgout)
                        lsc.util.updateheader(imgout, 0,
                                              {'template': (imgtemp0, 'template image'),
                                               'exptemp': (exp_temp,'exposure time template')})
                        lsc.util.updateheader(imgout, 0,
                                              {'target': (imgtarg0, 'target image'),
                                               'exptarg': (exp_targ,'exposure time terget')})

                        if args.normalize == 't':
                            lsc.util.updateheader(imgout, 0, {'EXPTIME': (exp_temp, '[s] Exposure length'),
                                                              'SATURATE': (sat_temp, '[ADU] Saturation level used')})

                        if 'CONVOL00' not in hd:
                            print '\n ### PSF computed by PyZOGY'
                        elif hd['CONVOL00'] == 'TEMPLATE':
                            print '\n ### image to compute  psf: '+imgtarg0
                            lsc.util.updateheader(imgout, 0, {'PSF': (imgtarg0, 'image to compute  psf')})
                        else:
                            print '\n ### image to compute  psf: '+imgtemp0
                            lsc.util.updateheader(imgout, 0, {'PSF': (imgtemp0, 'image to compute  psf')})

                            #                    copy all information from target
                        hd = fits.getheader(imgout)
                        dictionary = {}

                        try:
                            ggg0 = lsc.mysqldef.getfromdataraw(conn, 'photlco', 'filename', imgtarg0, '*')
                            for voce in ggg0[0].keys():
                                if voce not in ['id']:
                                    dictionary[voce] = ggg0[0][voce]
                        except:
                            dictionary = {'dateobs': lsc.util.readkey3(hd, 'date-obs'),
                                          'exptime': lsc.util.readkey3(hd, 'exptime'),
                                          'dayobs': lsc.util.readkey3(hd, 'day-obs'),
                                          'filter': lsc.util.readkey3(hd, 'filter'),
                                          'telescope': lsc.util.readkey3(hd, 'telescop'),
                                          'airmass': lsc.util.readkey3(hd, 'airmass'),
                                          'objname': lsc.util.readkey3(hd, 'object'),
                                          'wcs': lsc.util.readkey3(hd, 'wcserr'), 'ut': lsc.util.readkey3(hd, 'ut'),
                                          'mjd': lsc.util.readkey3(hd, 'mjd'),
                                          'instrument': lsc.util.readkey3(hd, 'instrume'),
                                          'ra0': lsc.util.readkey3(hd, 'RA'), 'dec0': lsc.util.readkey3(hd, 'DEC')}

                        dictionary['exptime'] =  lsc.util.readkey3(hd, 'exptime')
                        dictionary['psf'] = 'X'
                        dictionary['mag'] = 9999
                        dictionary['psfmag'] = 9999
                        dictionary['apmag'] = 9999
                        dictionary['z1'] = 9999
                        dictionary['z2'] = 9999
                        dictionary['c1'] = 9999
                        dictionary['c2'] = 9999
                        dictionary['filename'] = imgout0
                        dictionary['filepath'] = _dir
                        dictionary['filetype'] = 3
                        dictionary['difftype'] = args.difftype
                        if dictionary['filepath']:
                            if not os.path.isdir(dictionary['filepath']):
                                print dictionary['filepath']
                                os.mkdir(dictionary['filepath'])
                            if not os.path.isfile(dictionary['filepath'] + imgout0) or args.force in ['yes', True]:
                                os.system('mv -v ' + imgout + ' ' + dictionary['filepath'] + imgout0)
                                os.system('mv -v ' + imgtemp + ' ' + dictionary['filepath'] + imgout0.replace('.diff.', '.ref.'))
                                if args.difftype == '1':
                                    hdulist = fits.open(imgout.replace('.fits', '.psf.fits'))
                                    imgdata = hdulist[0].data
                                    yctr, xctr = np.array(imgdata.shape) / 2
                                    cutsize = 100
                                    hdulist[0].data = imgdata[yctr - cutsize : yctr + cutsize, xctr - cutsize : xctr + cutsize]
                                    psffile_fields = {'PIXSCALE': head_targ['PIXSCALE'],
                                                      'CRPIX1': cutsize, 'CRPIX2': cutsize, # make a fake WCS solution
                                                      'CRVAL1': 0, 'CRVAL2': 0,             # where we know the PSF star
                                                      'CD1_1': 1, 'CD2_2': 1,               # is at (0, 0), which is the
                                                      'CD1_2': 0, 'CD2_1': 0}               # center of the image
                                    if args.normalize == 't':
                                        psffile_fields['EXPTIME'] = head_temp['EXPTIME']
                                        psffile_fields['SATURATE'] = head_temp['SATURATE']
                                    elif args.normalize == 'i':
                                        psffile_fields['EXPTIME'] = head_targ['EXPTIME']
                                        psffile_fields['SATURATE'] = head_targ['SATURATE']
                                    hdulist[0].header.update(psffile_fields)
                                    hdulist.writeto(dictionary['filepath'] + imgout0.replace('.fits', '.zogypsf.fits'), overwrite=True)
                                    hdulist.close()
                        ###########################################################################################################
                        #                           choose sn2 file depending on
                        #                           normalization parameter
                        #
                        ##########################################################################################################
                        if args.normalize == 'i':
                            print '\n ### scale to target'
                            imgscale = imgtarg0
                            pathscale = _dir
                        elif args.normalize == 't':
                            print '\n ### scale to reference'
                            imgscale = imgtemp0
                            pathscale = _dirtemp

                        if os.path.isfile(pathscale + imgscale.replace('.fits', '.sn2.fits')):
                            line = ('cp ' + pathscale + imgscale.replace('.fits', '.sn2.fits') + ' ' +
                                    dictionary['filepath'] + imgout0.replace('.fits', '.sn2.fits'))
                            os.system(line)
                            print line
                            lsc.util.updateheader(dictionary['filepath'] + imgout0.replace('.fits', '.sn2.fits'), 0,
                                                  {'mag': (9999., 'apparent'), 'psfmag': (9999., 'inst mag'), 'apmag': (9999., 'aperture mag')})
                            #
                            # this is to keep track on the zeropoint and the flux measurement on the template image 
                            # in the difference image. the final flux is the sum of the flux in the difference image 
                            # + the flux on the reference image 
                            # -2.5 * log10( Ft/ exp_t ) = -2.5 log10 (Delta /exp_t + Fr/exp_r  * 10**(Zr - Zt)/(-2.5)) 
                            #
                            #  Ft  flux target
                            #  Fr flux reference
                            #  Zr zeropoint rerefence
                            #  Zt zeropoint target        
                            #  Delta   flux on difference image
                            #
                            if 'ZN' in hdtempsn:
                                lsc.util.updateheader(dictionary['filepath'] + imgout0.replace('.fits', '.sn2.fits'), 0,
                                                      {'ZNref': (hdtempsn['ZN'], 'ZN reference image')})
                            else:
                                print 'not ZN'
                            if 'apflux' in hdtempsn:
                                lsc.util.updateheader(dictionary['filepath'] + imgout0.replace('.fits', '.sn2.fits'), 0,
                                                      {'apfl1re': (hdtempsn['apflux'], 'flux reference image'),
                                                       'dapfl1re': (hdtempsn['dapflux'], 'error flux reference image')})
                            else:
                                print 'not apflux'

                        else:
                            print 'fits table not found ' + dictionary['filepath'] + imgscale.replace('.fits', '.sn2.fits')
                        if conn:
                            ###################    insert in photlco
                            ggg = lsc.mysqldef.getfromdataraw(conn, 'photlco', 'filename', imgout0, '*')
                            if ggg and args.force:
                                lsc.mysqldef.deleteredufromarchive(imgout0, 'photlco', 'filename')
                            if not ggg or args.force:
                                print 'insert'
                                print dictionary
                                lsc.mysqldef.insert_values(conn, 'photlco', dictionary)
                            else:
                                for voce in ggg[0].keys():
                                    if voce not in ['id']:
                                        lsc.mysqldef.updatevalue('photlco', voce, dictionary[voce], imgout0)
                            ggg = lsc.mysqldef.getfromdataraw(conn, 'photpairing', 'nameout', imgout0, '*')
                            if ggg:
                                lsc.mysqldef.deleteredufromarchive(imgout0, 'photpairing', 'nameout')
                            dictionary = {'namein': imgtarg0, 'nameout': imgout0,
                                          'nametemplate': imgtemp0, 'tablein': 'photlco',
                                          'tableout': 'photlco', 'tabletemplate': 'photlco'}
                            print 'insert in out'
                            print dictionary
                            lsc.mysqldef.insert_values(conn, 'photpairing', dictionary)
                        else:
                            print 'warning: mysql database not updated, you are not connected to the database'
