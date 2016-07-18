#!/usr/bin/env python
description = "> Ingest raw data "
usage = "%prog  instrument -e epoch"

import re
import sys
import string
import os
from optparse import OptionParser
import datetime
import lsc
from numpy import take, argsort

if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog 1.0")
    parser.add_option("-e", "--epoch", dest="epoch", default='20121212', type="str",
                      help='epoch to reduce  \t [%default]')
    parser.add_option("-T", "--telescope", dest="telescope", default='all', type="str",
                      help='-T telescope ' + ', '.join(lsc.telescope0['all']) + ', '.join(
                          lsc.site0) + ', fts, ftn, 1m0, '
                                       'kb, fl \t [%default]')
    parser.add_option("-R", "--RA", dest="ra", default='', type="str",
                      help='-R  ra    \t [%default]')
    parser.add_option("-D", "--DEC", dest="dec", default='', type="str",
                      help='-D dec   \t [%default]')
    parser.add_option("-n", "--name", dest="name", default='', type="str",
                      help='-n image name   \t [%default]')
    parser.add_option("-d", "--id", dest="id", default='', type="str",
                      help='-d identification id   \t [%default]')
    parser.add_option("-f", "--filter", dest="filter", default='', type="str",
                      help='-f filter [sloan,landolt,u,g,r,i,z,U,B,V,R,I] \t [%default]')
    parser.add_option("-F", "--force", dest="force", action="store_true")
    parser.add_option("-b", "--bad", dest="bad", default='', type="str",
                      help='-b bad stage [wcs,psf,psfmag,zcat,abscat,mag,goodcat,getmag,quality,cosmic] \t [%default]')
    parser.add_option("-s", "--stage", dest="stage", default='', type="str",
                      help='-s stage [checkwcs,checkpsf,checkmag,checkquality,checkpos,checkcat,'
                           'checkmissing,checkfvd,checkcosmic] \t [%default]')
    parser.add_option("--filetype", dest="filetype", default=1, type="int",
                      help='filetype  1 [single], 2 [merge], 3 differences \t [%default]')
    parser.add_option("--z1", dest="z1", default=None, type="int",
                      help='z1 \t [%default]')
    parser.add_option("--z2", dest="z2", default=None, type="int",
                      help='z2 \t [%default]')
    parser.add_option("--temptel", dest="temptel", default='', type="str",
                      help='--temptel  template instrument \t [%default]')

    option, args = parser.parse_args()
    _telescope = option.telescope
    if _telescope not in lsc.telescope0['all'] + lsc.site0 + ['all', 'ftn', 'fts', '1m0', 'kb', 'fl']:
        sys.argv.append('--help')
    option, args = parser.parse_args()
    _id = option.id
    _filter = option.filter
    _ra = option.ra
    _dec = option.dec
    _name = option.name
    _bad = option.bad
    _stage = option.stage
    _z1 = option.z1
    _z2 = option.z2
    _filetype = option.filetype
    _temptel = option.temptel
    if option.force == None:
        _redo = False
    else:
        _redo = True
    if _stage:
        if _stage not in ['checkfast', 'checkwcs', 'checkpsf', 'checkmag', 'checkquality', 'checkpos',
                          'checkcat', 'checkmissing', 'checkfvd', 'checkcosmic', 'checkdiff']:
            sys.argv.append('--help')
        if _stage=='checkdiff' and _filetype!=3: sys.exit('must run checkdiff on difference images (--filetype 3)')
    if _bad:
        if _bad not in ['wcs', 'psf', 'psfmag', 'zcat', 'abscat', 'mag', 'goodcat', 'quality', 'cosmic', 'diff']:
            sys.argv.append('--help')
    if _filter:
        if _filter not in ['landolt', 'sloan', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
            sys.argv.append('--help')
        else:
            try:
                _filter = lsc.sites.filterst(_telescope)[_filter]
            except:
                pass
    option, args = parser.parse_args()
    epoch = option.epoch
    if '-' not in str(epoch):
        # epoch0 = datetime.date(int(epoch[0:4]), int(epoch[4:6]), int(epoch[6:8]))
        lista = lsc.mysqldef.getlistfromraw(lsc.myloopdef.conn, 'photlco', 'dayobs', epoch, '', '*', _telescope)
        # lista=getfromdataraw(lsc.src.myloopdef.conn, 'photlco', 'dateobs',str(epoch0), 'all')
    else:
        epoch1, epoch2 = string.split(epoch, '-')
        start = datetime.date(int(epoch1[0:4]), int(epoch1[4:6]), int(epoch1[6:8]))
        stop = datetime.date(int(epoch2[0:4]), int(epoch2[4:6]), int(epoch2[6:8]))
        listepoch = [re.sub('-', '', str(i)) for i in
                     [start + datetime.timedelta(days=x) for x in range(0, 1 + (stop - start).days)]]
        lista = lsc.mysqldef.getlistfromraw(lsc.myloopdef.conn, 'photlco', 'dayobs', str(listepoch[0]),
                                            str(listepoch[-1]), '*', _telescope)

    if lista:
        ll0 = {}
        for jj in lista[0].keys():
            ll0[jj] = []
        for i in range(0, len(lista)):
            for jj in lista[0].keys():
                ll0[jj].append(lista[i][jj])
        inds = argsort(ll0['mjd'])  # sort by mjd
        for i in ll0.keys():
            ll0[i] = take(ll0[i], inds)
        ll = lsc.myloopdef.filtralist(ll0, _filter, _id, _name, _ra, _dec, _bad, _filetype, _temptel=_temptel)
        print '##' * 50
        print '# IMAGE                                    OBJECT           FILTER           WCS             PSF   ' \
              '        PSFMAG          ZCAT          MAG      ABSCAT'
        for i in range(0, len(ll['filename'])):
            print '%s\t%12s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s' % \
                  (str(ll['filename'][i]), str(ll['objname'][i]), str(ll['filter'][i]), str(ll['wcs'][i]),
                   str(ll['psf'][i]), str(ll['psfmag'][i]), str(ll['zcat'][i]), str(ll['mag'][i]), str(ll['abscat'][i]))

        print '\n###  total number = ' + str(len(ll['filename']))
        if _stage:
            print '##' * 50
            print _stage
            if _stage == 'checkpsf':
                lsc.myloopdef.checkpsf(ll['filename'])
            elif _stage == 'checkmag':
                lsc.myloopdef.checkmag(ll['filename'])
            elif _stage == 'checkwcs':
                lsc.myloopdef.checkwcs(ll['filename'], _redo, 'photlco', _z1, _z2)
            elif _stage == 'checkfast':
                lsc.myloopdef.checkfast(ll['filename'], _redo)
            elif _stage == 'checkquality':
                lsc.myloopdef.checkquality(ll['filename'])
            elif _stage == 'checkpos':
                lsc.myloopdef.checkpos(ll['filename'], _ra, _dec)
            elif _stage == 'checkcat':
                lsc.myloopdef.checkcat(ll['filename'])
            elif _stage == 'checkmissing':
                lsc.myloopdef.check_missing(ll['filename'])
            elif _stage == 'checkfvd':
                lsc.myloopdef.checkfilevsdatabase(ll)
            elif _stage == 'checkcosmic':
                lsc.myloopdef.checkcosmic(ll['filename'])
            elif _stage == 'checkdiff':
                lsc.myloopdef.checkdiff(ll['filename'])
            else:
                print _stage + ' not defined'
    else:
        print '\n### no data selected'
