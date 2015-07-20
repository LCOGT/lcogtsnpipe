#!/usr/bin/env python
description = ">> make template "
usage = "%prog image [options] "

from numpy import mean, array, compress, std, average, median, abs
import lsc
import os
import sys
import string
import re
from pyfits import getheader
from optparse import OptionParser
import math

hostname, username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)

if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-R", "--RA", dest="RA", default='',
                      type='str', help='RA coordinate \t\t\t %default')
    parser.add_option("-D", "--DEC", dest="DEC", default='',
                      type='str', help='DEC coordinate \t\t\t %default')
    parser.add_option("-p", "--psf", dest="psf", default='',
                      type='str', help='psf image \t\t\t %default')
    parser.add_option("-s", "--show", dest="show", action='store_true',
                      default=False, help='Show output \t\t [%default]')
    parser.add_option("-f", "--force", dest="force", action='store_true',
                      default=False, help='force archiving \t\t [%default]')
    parser.add_option("--mag", dest="mag", action='store_true',
                      default=False, help='chose mag interactively \t\t [%default]')
    parser.add_option("--clean", dest="clean", action='store_true',
                      default=False, help='clean from cosmic \t\t [%default]')
    parser.add_option("--pos", dest="pos", action='store_true',
                      default=False, help='chose mag interactively \t\t [%default]')
    option, args = parser.parse_args()
    if len(args) < 1: sys.argv.append('--help')
    option, args = parser.parse_args()
    imglist = lsc.util.readlist(args[0])
    imgdic = {}
    for img in imglist:
        hdr0 = lsc.util.readhdr(img)
        _filter0 = lsc.util.readkey3(hdr0, 'filter')
        if _filter0 not in imgdic:
            imgdic[_filter0] = {}
            imgdic[_filter0]['img'] = []
            imgdic[_filter0]['psf'] = []
        imgdic[_filter0]['img'].append(img)

    if option.psf:
        psflist = lsc.util.readlist(option.psf)
        for img in psflist:
            hdr0 = lsc.util.readhdr(img)
            _filter0 = lsc.util.readkey3(hdr0, 'filter')
            if _filter0 in imgdic:
                imgdic[_filter0]['psf'].append(img)
    else:
        psflist = ''
    _ra = option.RA
    _dec = option.DEC
    _show = option.show
    _force = option.force
    chosemag = option.mag
    chosepos = option.pos
    _clean = option.clean
    # #################################
    from pyraf import iraf
    from iraf import digiphot
    from iraf import daophot
    import pywcs
    ##################################
    goon = False
    for fil in imgdic:
        imglist1 = imgdic[fil]['img']
        for img0 in imglist1:
            if not imgdic[fil]['psf']:
                if os.path.exists(re.sub('.fits', '.psf.fits', img0)):
                    imgdic[fil]['psf'].append(re.sub('.fits', '.psf.fits', img0))

            if not imgdic[fil]['psf']:
                psfimg = ''
                print 'psf not found'
            else:
                print '\### found psffile'
                goon = True

            if goon:
                psfimg = imgdic[fil]['psf'][0]
                print img0, psfimg, _ra, _dec
                if os.path.exists(re.sub('.fits', '.sn2.fits', img0)):
                    hdr1 = lsc.util.readhdr(re.sub('.fits', '.sn2.fits', img0))
                    _xpos = lsc.util.readkey3(hdr1, 'PSFX1')
                    _ypos = lsc.util.readkey3(hdr1, 'PSFY1')
                    _mag = lsc.util.readkey3(hdr1, 'PSFMAG1')
                    if 'ZPN' in hdr1.keys():
                        _ZPN = lsc.util.readkey3(hdr1, 'ZPN')
                    else:
                        _ZPN = ''
                    _exptime = lsc.util.readkey3(hdr1, 'exptime')
                    print 'found fits table'
                    print _ZPN
                #######################  chose mag  ##############
                if not chosemag:
                    if _mag:
                        print 'use mag from fits table'
                        try:
                            _mag = float(_mag) + 2.5 * math.log10(float(_exptime))
                        except:
                            _mag = 0.0
                    else:
                        print 'mag not found, mag not input, no modification to the template '
                        _mag = 0
                else:
                    if _ZPN:
                        print '2 magnitude difference'
                        _magobj = raw_input('which is the mag of the object  ? ')
                        if float(_magobj) == 0.0:
                            _mag = 0.0
                        else:
                            _mag = float(_magobj) - (float(_ZPN))  #+2.5*math.log10(float(_exptime)))
                    else:
                        _mag0 = raw_input('which mag do you want to subtract  ' + str(_mag) + ' ? ')
                        if not _mag0:
                            _mag = float(_mag)
                        elif _mag0 == 0.0:
                            _mag = 0.0
                        else:
                            _mag = float(_mag0)
                    print _mag
                #######################  Chose Ra   and    Dec  ##############
                if not chosepos:
                    if not _ra and not _dec:
                        _ra, _dec, _SN0, _type = lsc.util.checksndb(img0, 'targets')
                    if _ra and _dec:
                        print 'use ra and dec from input database !!! '
                        hdr0 = lsc.util.readhdr(img0)
                        wcs = pywcs.WCS(hdr0)
                        pix1 = wcs.wcs_sky2pix(array(zip([_ra], [_dec]), float), 1)
                        _xpos, _ypos = pix1[0][0], pix1[0][1]
                else:
                    _xpos, _ypos = '', ''
                if _xpos and _ypos:
                    print 'use xpos and ypos from fits table'
                    print _xpos, _ypos
                else:
                    print 'no xpos and ypos define, do that interactively'
                    jj = 0
                    while not _xpos or not _ypos:
                        print "  MARK SN REGION WITH - x -, EXIT  - q -"
                        try:
                            lsc.util.delete('tmp.log')
                            zz1, zz2, goon = lsc.util.display_image(img, 1, '', '', True)
                            iraf.imexamine(img, 1, wcs='logical', logfile='tmp.log', keeplog=True)
                            xytargets = iraf.fields('tmp.log', '1,2', Stdout=1)
                            _xpos, _ypos = string.split(xytargets[0])[0], string.split(xytargets[0])[1]
                        except:
                            _xpos, _ypos = '', ''
                        jj = jj + 1
                        if jj > 10:
                            goon = False
                            break

            if goon:
                print _xpos, _ypos, _mag
                print img0, psfimg, _xpos, _ypos, _mag
                imgout = re.sub('.fits', '.temp.fits', string.split(img0, '/')[-1])
                lsc.util.delete('_tmp.fits,_tmp2.fits,_tmp2.fits.art,' + imgout)
                _targetid = lsc.mysqldef.targimg(img0)
                if _clean:
                    if os.path.isfile(re.sub('.fits', '.clean.fits', img0)):
                        img0 = re.sub('.fits', '.clean.fits', img0)
                if _show:
                    _z11, _z22, good = lsc.util.display_image(img0, 1, '', '', False)
                    z11 = float(_z11)
                    z22 = float(_z22)
                    answ = 'y'
                    answ = raw_input(">>>>> Cuts OK [y/n] [y]?")
                    if not answ:
                        answ = 'y'
                    elif answ == 'no':
                        answ = 'n'
                    while answ == 'n':
                        z11 = raw_input('>>> z1 = ? [' + str(_z11) + '] ? ')
                        z22 = raw_input('>>> z2 = ? [' + str(_z22) + '] ? ')
                        if not z11:
                            z11 = _z11
                        else:
                            z11 = float(z11)
                        if not z22:
                            z22 = _z22
                        else:
                            z22 = float(z22)
                        _z11, _z22, goon = lsc.util.display_image(img0, 1, z11, z22, False)
                        answ = raw_input(">>>>> Cuts OK [y/n] [y]?")
                        if not answ:
                            answ = 'y'
                        elif answ == 'no':
                            answ = 'n'
                        z11 = float(_z11)
                        z22 = float(_z22)
                else:
                    z11, z22 = '', ''
                answ = 'n'
                while answ == 'n':
                    coordlist = str(_xpos) + '   ' + str(_ypos) + '    ' + str(_mag)
                    os.system('echo ' + coordlist + ' > ddd')
                    iraf.imarith(img0, '-', img0, '_tmp.fits', verbose='no')
                    if float(_mag) != 0.0:
                        iraf.daophot.addstar("_tmp.fits", 'ddd', psfimg, "_tmp2.fits", nstar=1, veri='no', simple='yes',
                                             verb='no')
                        iraf.imarith(img0, '-', '_tmp2.fits', imgout, verbose='yes')
                    else:
                        print '\####  copy file '
                        iraf.imcopy(img0, imgout, verbose='yes')
                    if _show:
                        _z11, _z22, goon = lsc.util.display_image(imgout, 2, z11, z22, False)
                        answ = raw_input('ok  ? [[y]/n]')
                        if not answ:
                            answ = 'y'
                    else:
                        answ = 'y'
                    lsc.util.delete('_tmp.fits,_tmp2.fits,_tmp2.fits.art,ddd')
                    if answ == 'n':
                        lsc.util.delete(imgout)
                        _mag0 = raw_input('which magnitude  ' + str(_mag) + ' ?')
                        if _mag0:
                            _mag = _mag0
                print 'insert in the archive'
                hd = lsc.util.readhdr(imgout)
                dictionary = {'dateobs': lsc.util.readkey3(hd, 'date-obs'), 'exptime': lsc.util.readkey3(hd, 'exptime'),
                              'dayobs': lsc.util.readkey3(hd, 'day-obs'), 'filter': lsc.util.readkey3(hd, 'filter'),
                              'targetid': _targetid, 'mjd': lsc.util.readkey3(hd, 'mjd'),
                              'telescope': lsc.util.readkey3(hd, 'telescop'),
                              'airmass': lsc.util.readkey3(hd, 'airmass'),
                              'objname': lsc.util.readkey3(hd, 'object'), 'ut': lsc.util.readkey3(hd, 'ut'),
                              'wcs': lsc.util.readkey3(hd, 'wcserr'), 'instrument': lsc.util.readkey3(hd, 'instrume'),
                              'ra0': lsc.util.readkey3(hd, 'RA'), 'dec0': lsc.util.readkey3(hd, 'DEC')}
                dictionary['filename'] = string.split(imgout, '/')[-1]
                dictionary['filepath'] = '/science/supernova/data/lsc/' + lsc.util.readkey3(hd, 'date-night') + '/'
                dictionary['filetype'] = 4

                ###################    insert in photlco
                ggg = lsc.mysqldef.getfromdataraw(conn, 'photlco', 'filename', string.split(imgout, '/')[-1], '*')
                if ggg and _force:
                    lsc.mysqldef.deleteredufromarchive(string.split(imgout, '/')[-1], 'photlco', 'filename')
                if not ggg or _force:
                    print 'insert'
                    print dictionary
                    lsc.mysqldef.insert_values(conn, 'photlco', dictionary)
                else:
                    for voce in ggg[0].keys():
                        #                for voce in ['filetype','ra0','dec0']:
                        if voce in dictionary.keys():
                            lsc.mysqldef.updatevalue('photlco', voce, dictionary[voce], string.split(imgout, '/')[-1])
                if not os.path.isdir(dictionary['filepath']):
                    print dictionary['filepath']
                    os.mkdir(dictionary['filepath'])
                if not os.path.isfile(dictionary['filepath'] + imgout) or _force:
                    print 'mv ' + imgout + ' ' + dictionary['filepath'] + imgout
                    os.system('mv ' + imgout + ' ' + dictionary['filepath'] + imgout)
                print 'Warning: if this template still contain the target you want to measure, you need to run'
                print 'stages psf and apmag on the reference image'
