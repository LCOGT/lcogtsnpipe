#!/usr/bin/env python

from runlsc import run_cmd
from argparse import ArgumentParser


if __name__ == "__main__":
    basecmd = 'lscloop.py'
    logfile = open('ingesttemplates.log', 'w+')
   
    parser = ArgumentParser()
    parser.add_argument("-f", "--filter", default='')
    parser.add_argument("--targetid", default=None, type=int)
    parser.add_argument("-t", "--timeout", type=float, default=86400.)
    args = parser.parse_args()

    if args.filter:
        basecmd += ' -f ' + args.filter
    if args.targetid is not None:
        basecmd += ' --targetid ' + str(args.targetid)

    run_cmd(basecmd + ' -s ingestsloan', logfile, args.timeout)
    #run_cmd(basecmd + ' -s ingestps1', logfile, args.timeout)
    #run_cmd(basecmd + ' -s ingestdecam', logfile, args.timeout)

    ### Generate PSFs and run cosmic ray removal

    run_cmd(basecmd + ' -e 19990101-20191231 -s psf --filetype 4 --use-sextractor --fwhm 5 --max_apercorr 2', logfile, args.timeout)
    run_cmd(basecmd + ' -e 19990101-20191231 -s cosmic --filetype 4', logfile, args.timeout)
    
