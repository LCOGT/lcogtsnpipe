#!/usr/bin/env python

import datetime
import os
import sys
import re
import string
import time
import lsc
from numpy import compress

start = time.time()
from optparse import OptionParser

description = "> automatic process od lsc data "
usage = "%prog  fileter \t [listfile]"

coomandsn = {'LSQ12fxd': '-c -x 3 -y 3 --bkg 3 --size 6',
             'PSNJ081753': '-c -x 3 -y 3 --bkg 4 --size 8',
             'SN2013ak:': '-c -x 4 -y 4 --bkg 3 --size 4',
             'SN2013L': '-c -x 4 -y 4 --bkg 3 --size 4 --datamax 47000',
             'PSNJ0621': '-c -x 4 -y 4 --bkg 5 --size 8',
             'PSNJ1053': '-x 3 -y 3 --bkg 4 --size 8  --datamax 51000 --ref elp1m008-kb74-20130115-0322-e90.sn2.fits',
             'PSNJ0333': '-x 4 -y 4 --size 8 --bkg 4 --datamax 65000 --ref lsc1m005-kb78-20121102-0107-e90.sn2.fits',
             'PSNJ083745': '-c -x 3 -y 3 --bkg 3 --size 6',
             'PSNJ182501': '-c -x 5 -y 5 --bkg 3 --size 5',
             'PSNPGC027573': '-c -x 3 -y 3 --bkg 4 --size 6',
             'PSNJ081753': '-c -x 5 -y 5 --bkg 3 --size 5',
             'OGLE16': '-c -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13zm': '-c -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13aco': '-x 4 -y 4 --bkg 2 --size 4 --ref elp1m008-kb74-20130427-0131-e90.sn2.fits',
             'LSQ13acf': '-c -x 5 -y 5 --bkg 3 --size 4',
             'LSQ13aav': '-c -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13xh': '-c -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13xy': '-c -x 3 -y 3 --bkg 3 --size 4',
             'LSQ13xy': '-c -x 4 -y 4 --bkg 3 --size 4',
             'LSQ13abf': '-c -x 3 -y 3 --bkg 2 --size 4',
             'LSQ13aje': '-c -x 5 -y 5 --bkg 2 --size 4',
             'LSQ13aiz': '-c -x 3 -y 3 --bkg 4 --size 6',
             'LSQ13ajb': '-c -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13ajp': '-c -x 3 -y 3 --bkg 3 --size 6',
             'LSQ13ajo': '-c -x 3 -y 3 --bkg 3 --size 4',
             'LSQ13ajg': '-x 5 -y 5 --bkg 2 --size 4 -ref lsc1m005-kb78-20130515-0027-e90.sn2.fits',
             'PSN223702': '-x 3 -y 3 --bkg 2 --size 4 --ref elp1m008-kb74-20130502-0292-e90.sn2.fits',
             'PSN205408': '-c -x 3 -y 3 --bkg 3 --size 6',
             'PSN171722': '-c -x 3 -y 3 --bkg 3 --size 5',
             'PSNJ1739': '-x 5 -y 5 --bkg 2 --size 3 --ref lsc1m004-kb77-20130215-0225-e90.sn2.fits',
             'PSN132651': '-c -x 4 -y 4 --bkg 3 --size 5',
             'PSNJ111856': '-c -x 4 -y 4 --bkg 3 --size 4',
             'PSN205753': '-c -x 4 -y 4 --bkg 3 --size 4',
             'LSQ13bhl': '-c -x 3 -y 3 --bkg 3 --size 4',
             'LSQ13bgk': '-c -x 3 -y 3 --bkg 3 --size 4',
             'LSQ13bgg': '-c -x 3 -y 3 --bkg 3 --size 4',
             'SN2013en': '-c -x 3 -y 3 --bkg 3 --size 5',
             'PSN142115': '-c -x 4 -y 4 --bkg 3 --size 4',
             'PTF13ayw': '-c -x 4 -y 4 --bkg 3 --size 4',
             'PSN024508': '-c -x 4 -y 4 --bkg 3 --size 4',
             'PSN092656': ' -x 5 -y 5 --bkg 2 --size 3 --ref lsc1m009-kb73-20130515-0095-e90.sn2.fits',
             'PSN220221': '-c -x 5 -y 5 --bkg 3 --size 4',
             'PSN165902': '-c -x 3 -y 3 --bkg 4 --size 6',
             'PSN233746': '-c -x 3 -y 3 --bkg 4 --size 6',
             'PTF13dqy': '-c -x 3 -y 3 --bkg 4 --size 6',
             'PSNJ17492705+36083599': '-c -x 3 -y 3 --bkg 4 --size 6',
             'PTF14atg': '-c -x 3 -y 3 --bkg 4 --size 6',
             'PTF14ayo': '-c -x 3 -y 3 --bkg 4 --size 6',
             'lsq14eer': '-c -x 3 -y 3 --bkg 3 --size 5',
             'PSN09554214': '-c -x 3 -y 3 --bkg 3 --size 5',
             'asassn-15aj': '-c -x 4 -y 4 --bkg 3 --size 5',
             'LSQ13cnl': '-c -x 4 -y 4 --bkg 3 --size 5 --ref cpt1m012-kb75-20131008-0126-e90.sn2.fits',
             'SN 2016bkv': '-c -x 3 -y 3 --RAS 154.58051,154.58255 --DECS 41.427596,41.42898 --RA0 154.58144 --DEC0 41.428297',
}

###########################################################
if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog lsc")
    parser.add_option("-f", "--filter", dest="filter", default='', type="str",
                      help='-f filter [sloan,landolt,u,g,r,i,z,U,B,V,R,I] \t [%default]')
    parser.add_option("-T", "--telescope", dest="telescope", default='all', type="str")
    parser.add_option("-e", "--epoch", dest="epoch", default='', type="str",
                      help='epoch to reduce  \t [%default]')
    parser.add_option("-i", "--ingest", action="store_true",
                      dest='ingest', default=False, help='ingest new data \t\t\t [%default]')
    parser.add_option("--field", dest="field", default='', type="str",
                      help='--field  [landolt,sloan]  \t [%default]')
    parser.add_option("--filetype", dest="filetype", default=1, type="int",
                      help='filetype  1 [single], 2 [merge], 3 differences \t [%default]')

    option, args = parser.parse_args()
    #     _instrument=option.instrument
    _telescope = option.telescope
    _ingest = option.ingest
    _field = option.field
    _filetype = option.filetype
    _filter = option.filter

    if _filter:
        if _filter not in ['landolt', 'sloan', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I', 'ug',
                           'gr', 'gri', 'griz', 'riz', 'BVR', 'BV', 'BVRI', 'VRI']:
            sys.argv.append('--help')
    option, args = parser.parse_args()
    epoch = option.epoch
    if epoch:
        if '-' not in str(epoch):
            epoch0 = datetime.date(int(epoch[0:4]), int(epoch[4:6]), int(epoch[6:8]))
        else:
            epoch1, epoch2 = string.split(epoch, '-')
            epoch = epoch1 + '-' + epoch2
    else:
        d = datetime.date.today() + datetime.timedelta(1)
        g = d - datetime.timedelta(4)
        epoch = g.strftime("%Y%m%d") + '-' + d.strftime("%Y%m%d")

    if _field:
        fil = [_field]
    else:
        fil = ['landolt', 'sloan']

    if _filter:
        ff = ' -f ' + _filter
    else:
        ff = ''

    tt = '  -T ' + _telescope + ' '
    print '\n####  compute  astrometry, when missing '
    #  compute astrometry when tim astrometry failed
    os.system('lscloop.py -e ' + epoch + ' -b wcs -s wcs --mode astrometry ' + ff + tt)
    #  try again or set to bad image
    os.system('lscloop.py -e ' + epoch + ' -b wcs -s wcs --xshift 1 --yshift 1 ' + ff + tt)

    print '\n####  compute  psf, when missing, with catalog'
    os.system('lscloop.py -e ' + epoch + ' -b psf -s psf ' + ff + tt)  #  compute psf

    print '\n####  compute  psf, when missing, with sextractor'
    os.system('lscloop.py -e ' + epoch + ' -b psf -s psf --use-sextractor ' + ff + tt)  #  compute psf
#    print '\n####  compute  psf, when missing with fix fwhm'
#    os.system('lscloop.py -e ' + epoch + ' -b psf -s psf --fwhm 6 ' + ff + tt)  #  compute psf

    #############################################################################
    import glob

    ll = lsc.myloopdef.get_list(epoch, _telescope, '', '', '', '', '', '', 'photlco', _filetype)
    lista = list(set(ll['objname']))
    sloancal = []
    landoltcal = []
    apasscal = []
    standard = []
    notinthelist = []
    for obj in lista:
        img = compress(ll['objname'] == obj, ll['filename'])[0]
        _dir = compress(ll['objname'] == obj, ll['filepath'])[0]
        _SN0 = ''
        _Std = ''
        _ra0, _dec0, _Std = lsc.util.checksnlist(_dir + img, 'standardlist.txt')
        if not _Std:
            _ra0, _dec0, _SN0 = lsc.util.checksnlist(_dir + img, 'supernovaelist.txt')
        if not _Std and not _SN0:
            _ra0, _dec0, _SS, _type = lsc.util.checksndb(_dir + img, 'targets')
            if _SS:
                if _type in ['STD', 'std']:
                    _Std = _SS
                else:
                    _SN0 = _SS
        if _Std:
            print '\n### object in the standard list ' + _Std
            standard.append(obj)
            _temp = _Std
        else:
            _temp = _SN0
        if _temp:
            print '\n### object in the supernova list ' + _temp
            _catlandolt = lsc.util.getcatalog(_temp,'landolt')
            _catsloan = lsc.util.getcatalog(_temp,'sloan')
            _catapass = lsc.util.getcatalog(_temp,'apass')

            if len(_catlandolt) >= 1:
                _catlandolt = _catlandolt
                print '\n### landolt catalogue available ' + _catlandolt
            else:
                _catlandolt = ''
            if len(_catsloan) >= 1:
                _catsloan = _catsloan
                print '\n### sloan catalogue available ' + _catsloan
            else:
                _catsloan = ''
            if len(_catapass) >= 1:
                _catapass = _catapass
                print '\n### apass catalogue available ' + _catapass
            else:
                _catapass = ''
        if not _Std and not _SN0:
            print '\n###warning:  object NOT in the supernova and standard lists ' + obj
            notinthelist.append(obj)
            _catlandolt = ''
            _catsloan = ''
        if _catlandolt:
            landoltcal.append(obj)
        if _catsloan:
            sloancal.append(obj)
        if _catapass:
            apasscal.append(obj)

    print '\n### standard fields: ' + str(standard)
    print '\n### landolt fields:\n' + str(landoltcal)
    print '\n### sloan field:\n' + str(sloancal)
    print '\n### apass field:\n' + str(apasscal)
    print '\n### not in the lists:\n ' + str(notinthelist)
    #                 compute zero point for different fields
    for field in fil:
        if field == 'landolt':
            listacampi = standard + landoltcal
        if field == 'sloan':
            listacampi = standard + sloancal
        for _standard in listacampi:
            cmd = 'lscloop.py --type ph -b zcat -F -e ' + epoch + ' -n ' + _standard + ' -f ' + field + \
                  ' -s zcat --cutmag -6 -T ' + _telescope
            print cmd
            os.system(cmd)

    #    compute zeropoint for apass field not in landolt or sloan
    for field in fil:
        if field == 'landolt':
            for _name in apasscal:
                if _name not in landoltcal:
                    cmd = 'lscloop.py --type ph -b zcat -F -e ' + epoch + ' -n ' + _name + ' -f apass' + \
                          ' -s zcat --cutmag -6 -T ' + _telescope
                    print cmd
                    os.system(cmd)
        if field == 'sloan':
            for _name in apasscal:
                if _name not in sloancal:
                    cmd = 'lscloop.py --type ph -b zcat -F -e ' + epoch + ' -n ' + _name + ' -f apass' + \
                          ' -s zcat --cutmag -6 -T ' + _telescope
                    print cmd
                    os.system(cmd)

    #               compute catalogues for SN fields
    for field in fil:
        for obj in lista:
            if field == 'landolt':
                #                for obj in lista:
                #                if obj not in landoltcal:
                print obj + ': object not calibrated in landolt'
                for _std in standard:
                    os.system('lscloop.py --type ph  -F -e ' + epoch + ' -n ' + obj + ' -f ' + field + \
                              ' -b abscat -s abscat --standard ' + _std + ' -T ' + _telescope)
            if field == 'sloan':
                #               if obj not in sloancal:
                print obj + ': object not calibrated in sloan'
                for _std in standard:
                    os.system('lscloop.py --type ph  -F -e ' + epoch + ' -n ' + obj + ' -f ' + field + \
                              ' -b abscat -s abscat --standard ' + _std + ' -T ' + _telescope)
    # run psffit on all objects
    for field in fil:
        for obj in lista:
            if obj not in standard:
                print '###', obj
                if obj in coomandsn.keys():
                    os.system('lscloop.py -b psfmag -n ' + obj + ' -e ' + epoch + ' -f ' + field + \
                              ' -s psfmag -c  ' + coomandsn[obj])
                else:
                    os.system('lscloop.py -b psfmag -n ' + obj + ' -e ' + epoch + ' -f ' + field + \
                              ' -s psfmag -c  -x 3 -y 3 --bkg 4 --size 7 ')

    print '\n### landolt fields:\n' + str(landoltcal)
    print '\n### sloan field:\n' + str(sloancal)
    print '\n### apass field:\n' + str(apasscal)

    for field in fil:
        for obj in lista:
            if obj not in standard:
                if field == 'landolt':
                    if obj in landoltcal:
                        print '###', obj
                        try:
                            os.system('lscloop.py -b mag -n ' + obj + ' -e ' + epoch + ' -f ' + field + \
                                      ' -s mag --type fit')
                        except:
                            pass
                    elif obj in apasscal:
                        print '###', obj
                        try:
                            os.system('lscloop.py -b mag -n ' + obj + ' -e ' + epoch + ' -f apass' + \
                                      ' -s mag --type fit')
                        except:
                            pass
                elif field == 'sloan':
                    if obj in sloancal:
                        print '###', obj
                        try:
                            os.system('lscloop.py -b mag -n ' + obj + ' -e ' + epoch + ' -f ' + field + \
                                      ' -s mag --type fit')
                        except:
                            pass
                    elif obj in apasscal:
                        print '###', obj
                        try:
                            os.system('lscloop.py -b mag -n ' + obj + ' -e ' + epoch + ' -f apass' + \
                                      ' -s mag --type fit')
                        except:
                            pass

# cannot run this stage as a cronjob, because matplotlib needs $DISPLAY to be set
#    # make stamps for all new images
#    try:
#        if _telescope == 'all':
#            os.system('lscloop.py -e ' + str(epoch) + ' -s makestamp')
#        else:
#            os.system('lscloop.py -e ' + str(epoch) + ' -s makestamp' + ' -T '+_telescope)
#    except:
#        print 'warning makestap did not work'

    stop = time.time()
    print 'time to process all  data ', str(stop - start)
            
