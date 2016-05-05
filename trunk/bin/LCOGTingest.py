#!/usr/bin/env python
import requests
import lsc
import re
import os
import sys
from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u
from matplotlib.image import imsave
import numpy as np
from datetime import datetime
from glob import glob

def authenticate(username, password):
    '''Get the authentication token'''
    response = requests.post('https://archive-api.lcogt.net/api-token-auth/',
                             data = {'username': username, 'password': password}).json()
    token = response.get('token')
    authtoken = {'Authorization': 'Token ' + token}
    return authtoken

def get_metadata(authtoken={}, limit=None, **kwargs):
    '''Get the list of files meeting criteria in kwargs'''
    url = 'https://archive-api.lcogt.net/frames/?' + '&'.join(
            [key + '=' + str(val) for key, val in kwargs.items() if val is not None])
    url = url.replace('False', 'false')
    url = url.replace('True', 'true')
    print url

    response = requests.get(url, headers=authtoken).json()
    frames = response['results']
    while response['next'] and (limit is None or len(frames) < limit):
        print response['next']
        response = requests.get(response['next'], headers=authtoken).json()
        frames += response['results']
    return frames[:limit]

def download_frame(frame, force=False):
    '''Download a single image from the LCOGT archive and put it in the right directory'''
    filename = frame['filename']
    dayobs = re.search('(20\d\d)(0\d|1[0-2])([0-2]\d|3[01])', filename).group()
    if 'fs' in frame['INSTRUME']:
        daydir = 'data/fts/' + dayobs + '/'
    elif frame['INSTRUME'] == 'en06':
        daydir = 'data/floyds/' + dayobs + '_ftn/'
    elif frame['INSTRUME'] == 'en05':
        daydir = 'data/floyds/' + dayobs + '_fts/'
    else:
        daydir = 'data/lsc/' + dayobs + '/'
    filepath = lsc.util.workdirectory + daydir

    if not os.path.isdir(filepath):
        os.mkdir(filepath)

    basename = filename.replace('.fits.fz', '')[:-1] + '?.fits'
    matches = glob(filepath + basename) + glob(filepath + 'bad/' + basename)
    if not matches or force:
        print 'downloading', filename, 'to', filepath
        with open(filepath + filename, 'wb') as f:
            f.write(requests.get(frame['url']).content)
    else:
        matches_filenames = [fullpath.replace(filepath, '') for fullpath in matches]
        print matches_filenames, 'already in', filepath

    if os.path.isfile(filepath + filename) and os.stat(filepath + filename).st_size == 0:
        print filename, 'has size 0. Redownloading.'
        with open('filesize0.log', 'a') as l:
            l.write(str(datetime.utcnow()) + '\t' + filename)
        with open(filepath + filename, 'wb') as f:
            f.write(requests.get(frame['url']).content)

    if (filename[-3:] == '.fz' and os.path.isfile(filepath + filename)
            and (not os.path.isfile(filepath + filename[:-3]) or force)):
        print 'unpacking', filename
        os.system('funpack ' + filepath + filename)
        filename = filename[:-3]
    elif filename[-3:] == '.fz':
        print filename, 'already unpacked'
        filename = filename[:-3]

    return filepath, filename

hostname, username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)

programs = lsc.mysqldef.query(['select idcode, groupidcode from programs'], conn)
groupidcodes = {prog['idcode']: prog['groupidcode'] for prog in programs}

telescopes = lsc.mysqldef.query(['select id, name from telescopes'], conn)
telescopeids = {tel['name']: tel['id'] for tel in telescopes}

instruments = lsc.mysqldef.query(['select id, name from instruments'], conn)
instrumentids = {inst['name']: inst['id'] for inst in instruments}

photlcoraw_to_hdrkey = {'objname': 'OBJECT',
                        'dayobs': 'DAY-OBS',
                        'dateobs': 'DATE-OBS', # ignores everything after T
                        'ut': 'UTSTART', # ignores everything after .
                        'mjd': 'MJD-OBS',
                        'exptime': 'EXPTIME',
                        'filter': 'FILTER',
                        'telescope': 'TELESCOP',
                        'instrument': 'INSTRUME',
                        'airmass': 'AIRMASS',
                        'ra0': 'RA',
                        'dec0': 'DEC',
                        'cat_ra': 'CAT-RA',
                        'cat_dec': 'CAT-DEC',
                        'temperature': 'CCDATEMP',
                        'propid': 'PROPID',
                        'obid': 'BLKUID',
                        'userid': 'USERID',
                        'fwhm': 'L1FWHM',
                        'tracknumber': 'TRACKNUM',
                        'moonfrac': 'MOONFRAC',
                        'moondist': 'MOONDIST'}

speclcoraw_to_hdrkey = {'objname': 'OBJECT',
                        'dayobs': 'DAY-OBS',
                        'dateobs': 'DATE-OBS', # ignores everything after T
                        'ut': 'UTSTART', # ignores everything after .
                        'mjd': 'MJD-OBS',
                        'exptime': 'EXPTIME',
                        'filter': 'FILTER',
                        'telescope': 'TELESCOP',
                        'instrument': 'INSTRUME',
                        'type': 'OBSTYPE',
                        'airmass': 'AIRMASS',
                        'slit': 'APERWID',
                        'ra0': 'RA',
                        'dec0': 'DEC',
                        'cat_ra': 'CAT-RA',
                        'cat_dec': 'CAT-DEC',
                        'temperature': 'CCDATEMP',
                        'propid': 'PROPID',
                        'obid': 'BLKUID',
                        'rotskypa': 'ROTSKYPA',
                        'userid': 'USERID',
                        'fwhm': 'AGFWHM',
                        'tracknumber': 'TRACKNUM'}

def db_ingest(filepath, filename, force=False):
    '''Read an image header and add a row to the database'''
    global telescopeids, instrumentids
    if '-en0' in filename:
        table = 'speclcoraw'
        db_to_hdrkey = speclcoraw_to_hdrkey
    else:
        table = 'photlcoraw'
        db_to_hdrkey = photlcoraw_to_hdrkey
    fileindb = lsc.mysqldef.getfromdataraw(conn, table, 'filename', filename, column2='filepath')
    if fileindb:
        filepath = fileindb[0]['filepath'] # could be marked as bad
    if not fileindb or force:
        if filename[-3:] == '.fz':
            hdr = fits.getheader(filepath + filename, 1)
        else:
            hdr = fits.getheader(filepath + filename)
        dbdict = {'targetid': lsc.mysqldef.targimg(hdrt=hdr),
                    'filename': filename,
                    'filepath': filepath,
                    'groupidcode': groupidcodes[hdr['PROPID']]}
        for dbcol, hdrkey in db_to_hdrkey.items():
            if hdrkey in hdr and hdr[hdrkey] not in ['NaN', 'UNKNOWN', None, '']:
                if hdrkey in ['RA', 'CAT-RA']:
                    dbdict[dbcol] = Angle(hdr[hdrkey], u.hourangle).to_string(u.deg, decimal=True, precision=7)
                elif hdrkey in ['DEC', 'CAT-DEC']:
                    dbdict[dbcol] = Angle(hdr[hdrkey], u.deg).to_string(decimal=True, precision=7)
                else:
                    dbdict[dbcol] = hdr[hdrkey]
        if hdr['TELESCOP'] not in telescopeids:
            print hdr['TELESCOP'], 'not recognized. Adding to telescopes table.'
            lsc.mysqldef.insert_values(conn, 'telescopes', {'name': hdr['TELESCOP']})
            telescopes = lsc.mysqldef.query(['select id, name from telescopes'], conn)
            telescopeids = {tel['name']: tel['id'] for tel in telescopes}
        dbdict['telescopeid'] = telescopeids[hdr['TELESCOP']]
        if hdr['INSTRUME'] not in instrumentids:
            print hdr['INSTRUME'], 'not recognized. Adding to instruments table.'
            lsc.mysqldef.insert_values(conn, 'instruments', {'telescopeid': telescopeids[hdr['TELESCOP']], 'name': hdr['INSTRUME']})
            instruments = lsc.mysqldef.query(['select id, name from instruments'], conn)
            instrumentids = {inst['name']: inst['id'] for inst in instruments}
        dbdict['instrumentid'] = instrumentids[hdr['INSTRUME']]
        if fileindb:
            lsc.mysqldef.query(["delete from " + table + " where filename='" + filename + "'"], conn)
        print 'ingesting', filename
        lsc.mysqldef.insert_values(conn, table, dbdict)
    else:
        print filename, 'already ingested'

def fits2png(filename, zclip=5):
    data = fits.getdata(filename)
    z1 = np.percentile(data, zclip)
    z2 = np.percentile(data, 100-zclip)
    imsave(filename.replace('.fits', '.png'), data, cmap='gray', vmin=z1, vmax=z2, origin='lower')

#################################################################
if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Downloads LCOGT data from archive.lcogt.net')
    parser.add_argument("-u", "--username")
    parser.add_argument("-p", "--password")
    parser.add_argument("-l", "--limit", type=int, help="maximum number of frames to return")
    parser.add_argument("-F", "--force-dl", action="store_true", help="download files even if they already exist")
    parser.add_argument("-G", "--force-db", action="store_true", help="reingest files even if they already exist")
    parser.add_argument("-H", "--force-tn", action="store_true", help="regenerate thumbnails even if they already exist")

    parser.add_argument("-S", "--site", choices=['bpl', 'coj', 'cpt', 'elp', 'lsc', 'ogg', 'sqa', 'tfn'])
    parser.add_argument("-T", "--telescope", choices=['0m4a', '0m4b', '0m4c', '0m8a', '1m0a', '2m0a'])
    parser.add_argument("-I", "--instrument")
    parser.add_argument("-f", "--filter", choices=['up', 'gp', 'rp', 'ip', 'zs', 'U', 'B', 'V', 'R', 'I'])
    parser.add_argument("-P", "--proposal", help="proposal ID (PROPID in the header)")
    parser.add_argument("-n", "--name", help="target name")
    parser.add_argument("-s", "--start", help="start date")
    parser.add_argument("-e", "--end", help="end date")

    parser.add_argument("-t", "--obstype", choices=['ARC', 'BIAS', 'CATALOG', 'DARK', 'EXPERIMENTAL',
                                        'EXPOSE', 'LAMPFLAT', 'SKYFLAT', 'SPECTRUM', 'STANDARD'])
    parser.add_argument("-r", "--reduction", choices=['raw', 'quicklook', 'reduced'])
    parser.add_argument("-x", "--experimental", action='store_true', help="get products from Curtis's pipeline")
    parser.add_argument("--public", action='store_true', help="include public data")

    args = parser.parse_args()

    if args.username and args.password:
        authtoken = authenticate(args.username, args.password)
    else:
        authtoken = {}

    if args.reduction == 'raw':
        rlevel = 0
    elif args.reduction == 'quicklook' and args.experimental:
        rlevel = 11
    elif args.reduction == 'quicklook':
        rlevel = 10
    elif args.reduction == 'reduced' and args.experimental:
        rlevel = 91
    elif args.reduction == 'reduced':
        rlevel = 90
    else:
        rlevel = None

    frames = get_metadata(authtoken, limit=args.limit, SITEID=args.site, TELID=args.telescope,
                        INSTRUME=args.instrument, FILTER=args.filter, PROPID=args.proposal, OBJECT=args.name,
                        start=args.start, end=args.end, OBSTYPE=args.obstype, RLEVEL=rlevel, public=args.public)

    print 'Total number of frames:', len(frames)

    for frame in frames:
        filepath, filename = download_frame(frame, args.force_dl)
        db_ingest(filepath, filename, args.force_db)
        if '-en0' in filename and '-e00.fits' in filename and (not os.path.isfile(filepath + filename.replace('.fits', '.png')) or args.force_tn):
            fits2png(filepath + filename)
