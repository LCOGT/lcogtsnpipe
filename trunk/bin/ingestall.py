#!/usr/bin/env python

import lsc
import base64
from LCOGTingest import *
from datetime import datetime, timedelta
import os

now = datetime.utcnow()
start = datetime.strftime(now - timedelta(days=7), '%Y-%m-%d') # check for data in the last week
end = datetime.strftime(now, '%Y%m%d')

login = lsc.mysqldef.query(["select username, userpw from programs where idcode='KEY2014A-003'"], conn)
username = login[0]['username']
password = base64.decodestring(login[0]['userpw'])
authtoken = authenticate(username, password)

frames = get_metadata(authtoken, start=start, OBSTYPE='EXPOSE', RLEVEL=91, public=False)    # all images where SNEx is a co-I
frames += get_metadata(authtoken, start=start, TELID='1m0a', OBSTYPE='STANDARD', RLEVEL=91) # all 1m standard fields
frames += get_metadata(authtoken, start=start, TELID='2m0a', OBSTYPE='STANDARD', RLEVEL=91) # all 2m standard fields
frames += get_metadata(authtoken, start=start, INSTRUME='en06', RLEVEL=0, public=False)     # all FTN spectra SNEx is a co-I
frames += get_metadata(authtoken, start=start, INSTRUME='en05', RLEVEL=0, public=False)     # all FTS spectra SNEx is a co-I
frames += get_metadata(authtoken, start=start, PROPID='OGG_calib', RLEVEL=0)                # FTN standard star spectra
frames += get_metadata(authtoken, start=start, PROPID='COJ_calib', RLEVEL=0)                # FTS standard star spectra

print 'Total number of frames:', len(frames)

for frame in frames:
    try:
        filepath, filename = download_frame(frame)
    except Exception as e:
        print e
        print '!!! FAILED TO DOWNLOAD', frame['filename']
        continue
    try:
        db_ingest(filepath, filename)
    except Exception as e:
        print e
        print '!!! FAILED TO INGEST', filename
        continue
    if '-en0' in filename and '-e00.fits' in filename and not os.path.isfile(filepath + filename.replace('.fits', '.png')):
        try:
            fits2png(filepath + filename)
        except Exception as e:
            print e
            print '!!! FAILED TO MAKE PNG FOR', filename

os.system('lscingestredudata.py -e ' + start.replace('-', '') + '-' + end) # ingest new data into photlco
