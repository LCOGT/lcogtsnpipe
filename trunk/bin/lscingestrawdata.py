#!/usr/bin/env python

import re
import sys
import glob
import string
import os
from optparse import OptionParser
import MySQLdb
import datetime
import lsc

description = "> Ingest raw data "
usage = "%prog  instrument -e epoch"

if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog 1.0")
    parser.add_option("-T", "--telescope", dest="telescope", default='all', type="str",
                      help='name -T telescope elp,fts, ftn, lsc, cpt, all, coj, ogg, tar   \t [%default]')
    parser.add_option("--list", dest="lista", default='', type="str",
                      help='list of files to ingest    \t [%default]')
    parser.add_option("-i", "--instrument", dest="instrument", default='', type="str",
                      help='name -i instrument ' + ', '.join(lsc.util.instrument0['all']) + ' \t [%default]')
    parser.add_option("-f", "--force", dest="force", action="store_true")
    parser.add_option("-t", "--type", dest="type", default='oracproc', type="str",
                      help='type -t type \t [oracproc/quicklook/pylcogt] \n')
    parser.add_option("-e", "--epoch", dest="epoch", default='20121212', type="str",
                      help='-e epoch \t [%default]')
    parser.add_option("--object", dest="object", default='', type="str",
                      help='type --object object \t [name] \n')

    option, args = parser.parse_args()
    _type = option.type
    _list = option.lista 
    _instrument = option.instrument
    _telescope = option.telescope
    _object = option.object

#    if option.type not in ['raw', 'redu']:
#        sys.argv.append('--help')

    if len(sys.argv) < 2:
        sys.argv.append('--help')

    if _type not in ['oracproc','quicklook','pylcogt']:
        sys.argv.append('--help')
    if _list:
        lista = string.split(_list,',')
    else:
        lista=''

    option, args = parser.parse_args()

    epoch = option.epoch
    if '-' not in str(epoch):
        epoch0 = datetime.date(int(epoch[0:4]), int(epoch[4:6]), int(epoch[6:8]))
        listepoch = [re.sub('-', '', str(epoch0))]
    else:
        epoch1, epoch2 = string.split(epoch, '-')
        start = datetime.date(int(epoch1[0:4]), int(epoch1[4:6]), int(epoch1[6:8]))
        stop = datetime.date(int(epoch2[0:4]), int(epoch2[4:6]), int(epoch2[6:8]))
        listepoch = [re.sub('-', '', str(i)) for i in
                     [start + datetime.timedelta(days=x) for x in range(0, (stop - start).days)]]

    _force = option.force
    if not _force:
        _force = False

    print _instrument, listepoch, _force, _type, _telescope
    lsc.mysqldef.ingestdata(_telescope, _instrument, listepoch, _force, _type,_object,lista)
