#!/usr/bin/env python

import lsc
import base64
from LCOGTingest import *
from datetime import datetime, timedelta
import os
import sys
import traceback
import logging

logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')

handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

if len(sys.argv) > 1:
    daterange = sys.argv[1]
    start, end = daterange.split('-')
    start = start[0:4] + '-' + start[4:6] + '-' + start[6:8]
    end = end[0:4] + '-' + end[4:6] + '-' + end[6:8]
else:
    now = datetime.utcnow()
    start = datetime.strftime(now - timedelta(days=7), '%Y-%m-%d') # check for data in the last week
    end = datetime.strftime(now, '%Y-%m-%d %H:%M:%S')
    daterange = datetime.strftime(now - timedelta(days=7), '%Y%m%d') + '-' + datetime.strftime(now, '%Y%m%d')

login = lsc.mysqldef.query(["select username, userpw from programs where idcode='KEY2014A-003'"], conn)
username = login[0]['username']
password = base64.decodestring(login[0]['userpw'])
try:
    authtoken = authenticate(username, password)
except ValueError as e:
    logger.info('Getting throttled: {error}'.format(error = str(e)))
    sys.exit()

frames = get_metadata(authtoken, start=start, end=end, OBSTYPE='EXPOSE', RLEVEL=91, public=False)       # all images where SNEx is a co-I
for telid in ['2m0a', '1m0a', '0m4a', '0m4b', '0m4c']:
    frames += get_metadata(authtoken, start=start, end=end, PROPID='Photometric standards', OBSTYPE='STANDARD', TELID=telid, RLEVEL=91)  # all photometric standards (except SQA)
    frames += get_metadata(authtoken, start=start, end=end, PROPID='Photometric standards', OBSTYPE='EXPOSE', TELID=telid, RLEVEL=91)  # all photometric standards (except SQA)
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en06', basename = 'e91-1d', public=False)        # all FTN spectra SNEx is a co-I, checks basename for Banzai-floyds 1ds
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en05',  basename = 'e91-1d', public=False)        # all FTS spectra SNEx is a co-I, checks basename for Banzai-floyds 1ds
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en12',  basename = 'e91-1d', public=False)        # all FTS spectra SNEx is a co-I, checks basename for Banzai-floyds 1ds

frames = get_metadata(authtoken, start=start, end=end, INSTRUME='en06', RLEVEL=0, public=False)        # all FTN spectra SNEx is a co-I, keeping e00 for floyds inbox (IRAF) reductions
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en05',  RLEVEL=0, public=False)        # all FTS spectra SNEx is a co-I, keeping e00 for floyds inbox (IRAF) reductions
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en12',  RLEVEL=0, public=False)        # all FTS spectra SNEx is a co-I, keeping e00 for floyds inbox (IRAF) reductions

frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en06', RLEVEL=0, PROPID='FLOYDS standards')  # FTN standard star spectra
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en05', RLEVEL=0, PROPID='FLOYDS standards')  # FTS standard star spectra
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en12', RLEVEL=0, PROPID='FLOYDS standards')  # FTS standard star spectra
frames += get_metadata(authtoken, start=start, end=end, PROPID='LCOEPO2016B-001', RLEVEL=91)            # supernova tracker images

logger.info('Total number of frames: {:d}'.format(len(frames)))

fullpaths = []
for frame in frames:
    try:
        filepath, filename = download_frame(frame)
        if '-en' not in filename:
            fullpaths.append(filepath + filename)
    except:
        logger.error('!!! FAILED TO DOWNLOAD {}'.format(frame['filename']))
        traceback.print_exc()
        continue
    try:
        if '-en' in filename: #ingesting the 1d spectra into both speclcoraw and spec
            table_raw = 'speclcoraw' # this table is what the scheduler checks, shouldn't be updated with versions
            dbdict = db_ingest(filepath, filename,table_raw)
        if 'e91-1d' in filename: # only ingest the banzai floyds directly into spec
            table_reduced = 'spec' # this table will track version control of spectra that get re-reduced
            dbdict = db_ingest(filepath, filename,table_reduced)
        if '-en' not in filename and 'e91-1d' not in filename: # all others are photlcoraw
            table = 'photlcoraw'
            dbdict = db_ingest(filepath, filename,table)
    except:
        logger.error('!!! FAILED TO INGEST {}'.format(filename))
        traceback.print_exc()
        continue
    if '-en' in filename and '-e00' in filename: #ignoring banzai files
        try:
            fits2png(filepath + filename)
        except:
            logger.error('!!! FAILED TO MAKE PNG FOR {}'.format(filename))
            traceback.print_exc()

lsc.mysqldef.ingestredu(fullpaths) # ingest new data into photlco

# add links to FLOYDS guider frames in database
frames = get_metadata(authtoken, start=start, end=end, INSTRUME='en06', RLEVEL=90, public=False)  # all FTN spectra SNEx is a co-I
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en05', RLEVEL=90, public=False) # all FTS spectra SNEx is a co-I
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en12', RLEVEL=90, public=False) # all FTS spectra SNEx is a co-I
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en06', RLEVEL=90, PROPID='FLOYDS standards')  # FTN standard star spectra
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en05', RLEVEL=90, PROPID='FLOYDS standards')  # FTS standard star spectra
frames += get_metadata(authtoken, start=start, end=end, INSTRUME='en12', RLEVEL=90, PROPID='FLOYDS standards')  # FTS standard star spectra

for frame in frames:
    try:
        record_floyds_tar_link(authtoken, frame)
    except:
        logger.error('!!! FAILED TO RECORD GUIDER LINK FOR {}'.format(frame['filename']))
        traceback.print_exc()
