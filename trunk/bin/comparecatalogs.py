#!/usr/bin/env python

import lsc
import os
import argparse

default_catdir = os.path.join(os.getenv('LCOSNDIR', lsc.util.workdirectory), 'standard', 'cat')

parser = argparse.ArgumentParser()
parser.add_argument('-F', '--force', action='store_true', help="try to download catalog even if we've tried before")
parser.add_argument('-R', '--radius', type=float, default=20., help="query for stars within this radius of the coordinates")
parser.add_argument('-f', '--field', nargs='+', choices=['landolt', 'apass', 'sloan', 'gaia'],
        default=['landolt', 'apass', 'sloan', 'gaia'], help="catalogs to check (all by default)")
parser.add_argument('-p', '--panstarrs', action='store_true', help='use a Pan-STARRS1 3pi catalog for the SDSS filters')
args = parser.parse_args()

if args.force:
    cond = '=""' # empty string means we've checked before and it wasn't there
else:
    cond = ' is null' # null means we haven't checked yet

for system in args.field:
    targets = lsc.mysqldef.query(['select id, ra0, dec0 from targets where ' + system + '_cat' + cond], lsc.conn)
    print len(targets), 'targets with no', system, 'catalog'
    print
    source = 'panstarrs' if system == 'sloan' and args.panstarrs else system
    for target in targets:
        tid = str(target['id'])
        print 'Target ID', tid
        names = lsc.mysqldef.query(['select name from targetnames where targetid=' + tid], lsc.conn)
        print ' aka '.join([n['name'] for n in names])
        filepath = os.path.join(default_catdir, system)
        for name in names: # see en.wikipedia.org/wiki/Filename#Comparison_of_filename_limitations
            filename = name['name'].translate(None, ' "*:<>?/|\\') + '_' + source + '.cat'
            fullpath = os.path.join(filepath,  filename)
            fileexists = os.path.isfile(fullpath)
            if fileexists:
                print 'Adding', filename, 'to database'
                lsc.mysqldef.query(['update targets set ' + system + '_cat="' + filename + '", '+ system + '_cat_path  where id=' + tid], lsc.conn)
                break
        if not fileexists and system != 'landolt':
            print 'Querying for catalog...'
            if system == 'apass':
                os.system('queryapasscat.py -r {ra0} -d {dec0} -R {radius} -o '.format(radius=args.radius, **target) + fullpath)
            elif system == 'sloan' and args.panstarrs:
                lsc.lscabsphotdef.panstarrs2file(target['ra0'], target['dec0'], args.radius, output=fullpath)
            elif system == 'sloan':
                lsc.lscabsphotdef.sloan2file(target['ra0'], target['dec0'], args.radius, output=fullpath)
            elif system == 'gaia':
                lsc.lscabsphotdef.gaia2file(target['ra0'], target['dec0'], output=fullpath)
            fileexists = os.path.isfile(fullpath)
            if fileexists:
                print 'Adding', filename, 'to database'
                lsc.mysqldef.query(['update targets set ' + system + '_cat="' + filename + '" where id=' + tid], lsc.conn)
        if not fileexists:
            print 'No catalog exists'
        print
        lsc.mysqldef.query(['update targets set ' + system + '_cat="" where ' + system + '_cat is NULL'], lsc.conn) # change NULL to '' if not in field
