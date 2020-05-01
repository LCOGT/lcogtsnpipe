#!/usr/bin/env python

import lsc
import sys
import os
import re
import string
import glob
import subprocess
import zipfile
from optparse import OptionParser
_dir = os.environ['LCOSNDIR']

def unpack_zip(_zipfile):
    with zipfile.ZipFile(_zipfile) as myzip:
        imglist = myzip.namelist()
        myzip.extractall()
        return imglist

def ingest_files(imglist,force=False):
    if force:       
        lsc.mysqldef.ingestredu(imglist, 'yes') # ingest new data into photlco
    else:
        lsc.mysqldef.ingestredu(imglist, 'no') # ingest new data into photlco

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Downloads data from SNEX')
    parser.add_argument("-f", "--file")
    parser.add_argument("-G", "--force-db", action="store_true", help="reingest files even if they already exist")
    args = parser.parse_args()
    if args.file:
        if 'zip' in args.file:
            imglist = unpack_zip(args.file)
        else:
            imglist_fz = glob.glob(os.path.join(args.file, '*.fz'))
            if len(imglist_fz)>0:
                for ifile in imglist_fz:
                    subprocess.Popen("echo 'funpack {}'".format(ifile), shell=True)
                    process = subprocess.Popen("funpack {}".format(ifile), shell=True)
                    process.wait()
            imglist = glob.glob(os.path.join(args.file, '*.fits'))
        ingest_files(imglist, force=args.force_db)
    else:
        print("you must specify a tar or zip file or directory with the -f or --file")
    
        
