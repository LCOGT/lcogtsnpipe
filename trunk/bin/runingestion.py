#!/opt/epd/bin/python

import datetime
import os,sys,re,string
import time
from optparse import OptionParser

description="> ingest phtometric data  " 
usage= "%prog  -e epoch \n"

def runin(telescope,instrument,epoch,_type='raw'):
    import os
    if instrument: ii=' -i '+instrument
    else:          ii=' '
    try:
        if _type=='raw':
            print 'lscingestrawdata.py -T '+telescope+ii+' -e '+epoch
            os.system('lscingestrawdata.py -T '+telescope+ii+' -e '+epoch)
        elif _type=='redu':
            print 'lscingestredudata.py -e '+epoch+' -m -T '+telescope+ii
            os.system('lscingestredudata.py -e '+epoch+' -m -T '+telescope+ii)
        else: print 'type not defined'
    except:
        print 'problem with ingestion'

if __name__ == "__main__":
    d=datetime.date.today()+datetime.timedelta(1)
    g=d-datetime.timedelta(4)
    epoch=g.strftime("%Y%m%d")+'-'+d.strftime("%Y%m%d")
    tel=['1m0-03', '1m0-04','1m0-05', '1m0-08', '1m0-09', '1m0-10', '1m0-11', '1m0-12', '1m0-13','2m0-01','2m0-02']
    ###########################   ingest   1m data
    for i in range(0,len(tel)):    runin(tel[i],'',epoch,'raw')  # ingest raw date
    print '\n### ingest redu data'
    for i in range(0,len(tel)):    runin(tel[i],'',epoch,'redu') # ingest redu data
    ###########################   ingest   FT data
