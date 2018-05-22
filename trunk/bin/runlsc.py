#!/usr/bin/env python

from datetime import datetime, timedelta
import os
import time
import lsc
import subprocess32
from argparse import ArgumentParser

commandsn = {'LSQ12fxd': ' -x 3 -y 3 --bkg 3 --size 6',
             'PSNJ081753': ' -x 3 -y 3 --bkg 4 --size 8',
             'SN2013ak:': ' -x 4 -y 4 --bkg 3 --size 4',
             'SN2013L': ' -x 4 -y 4 --bkg 3 --size 4 --datamax 47000',
             'PSNJ0621': ' -x 4 -y 4 --bkg 5 --size 8',
             'PSNJ1053': ' -x 3 -y 3 --bkg 4 --size 8  --datamax 51000 --ref elp1m008-kb74-20130115-0322-e90.sn2.fits',
             'PSNJ0333': ' -x 4 -y 4 --size 8 --bkg 4 --datamax 65000 --ref lsc1m005-kb78-20121102-0107-e90.sn2.fits',
             'PSNJ083745': ' -x 3 -y 3 --bkg 3 --size 6',
             'PSNJ182501': ' -x 5 -y 5 --bkg 3 --size 5',
             'PSNPGC027573': ' -x 3 -y 3 --bkg 4 --size 6',
             'PSNJ081753': ' -x 5 -y 5 --bkg 3 --size 5',
             'OGLE16': ' -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13zm': ' -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13aco': ' -x 4 -y 4 --bkg 2 --size 4 --ref elp1m008-kb74-20130427-0131-e90.sn2.fits',
             'LSQ13acf': ' -x 5 -y 5 --bkg 3 --size 4',
             'LSQ13aav': ' -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13xh': ' -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13xy': ' -x 3 -y 3 --bkg 3 --size 4',
             'LSQ13xy': ' -x 4 -y 4 --bkg 3 --size 4',
             'LSQ13abf': ' -x 3 -y 3 --bkg 2 --size 4',
             'LSQ13aje': ' -x 5 -y 5 --bkg 2 --size 4',
             'LSQ13aiz': ' -x 3 -y 3 --bkg 4 --size 6',
             'LSQ13ajb': ' -x 3 -y 3 --bkg 3 --size 5',
             'LSQ13ajp': ' -x 3 -y 3 --bkg 3 --size 6',
             'LSQ13ajo': ' -x 3 -y 3 --bkg 3 --size 4',
             'LSQ13ajg': ' -x 5 -y 5 --bkg 2 --size 4 -ref lsc1m005-kb78-20130515-0027-e90.sn2.fits',
             'PSN223702': ' -x 3 -y 3 --bkg 2 --size 4 --ref elp1m008-kb74-20130502-0292-e90.sn2.fits',
             'PSN205408': ' -x 3 -y 3 --bkg 3 --size 6',
             'PSN171722': ' -x 3 -y 3 --bkg 3 --size 5',
             'PSNJ1739': ' -x 5 -y 5 --bkg 2 --size 3 --ref lsc1m004-kb77-20130215-0225-e90.sn2.fits',
             'PSN132651': ' -x 4 -y 4 --bkg 3 --size 5',
             'PSNJ111856': ' -x 4 -y 4 --bkg 3 --size 4',
             'PSN205753': ' -x 4 -y 4 --bkg 3 --size 4',
             'LSQ13bhl': ' -x 3 -y 3 --bkg 3 --size 4',
             'LSQ13bgk': ' -x 3 -y 3 --bkg 3 --size 4',
             'LSQ13bgg': ' -x 3 -y 3 --bkg 3 --size 4',
             'SN2013en': ' -x 3 -y 3 --bkg 3 --size 5',
             'PSN142115': ' -x 4 -y 4 --bkg 3 --size 4',
             'PTF13ayw': ' -x 4 -y 4 --bkg 3 --size 4',
             'PSN024508': ' -x 4 -y 4 --bkg 3 --size 4',
             'PSN092656': ' -x 5 -y 5 --bkg 2 --size 3 --ref lsc1m009-kb73-20130515-0095-e90.sn2.fits',
             'PSN220221': ' -x 5 -y 5 --bkg 3 --size 4',
             'PSN165902': ' -x 3 -y 3 --bkg 4 --size 6',
             'PSN233746': ' -x 3 -y 3 --bkg 4 --size 6',
             'PTF13dqy': ' -x 3 -y 3 --bkg 4 --size 6',
             'PSNJ17492705+36083599': ' -x 3 -y 3 --bkg 4 --size 6',
             'PTF14atg': ' -x 3 -y 3 --bkg 4 --size 6',
             'PTF14ayo': ' -x 3 -y 3 --bkg 4 --size 6',
             'lsq14eer': ' -x 3 -y 3 --bkg 3 --size 5',
             'PSN09554214': ' -x 3 -y 3 --bkg 3 --size 5',
             'asassn-15aj': ' -x 4 -y 4 --bkg 3 --size 5',
             'LSQ13cnl': ' -x 4 -y 4 --bkg 3 --size 5 --ref cpt1m012-kb75-20131008-0126-e90.sn2.fits',
             'SN 2016bkv': ' -x 3 -y 3 --RAS 154.58051,154.58255 --DECS 41.427596,41.42898 --RA0 154.58144 --DEC0 41.428297',
}

def run_cmd(cmd, logfile=None, timeout=3600):
    logfile.write(cmd + '\n')
    logfile.flush()
    proc = subprocess32.Popen(cmd, shell=True, stdout=logfile, stderr=logfile)
    try:
        proc.wait(timeout)
    except subprocess32.TimeoutExpired as e:
        proc.kill()
        logfile.close()
        with open(logfile.name, 'r') as f:
            print f.read()
        raise e

###########################################################
if __name__ == "__main__":
    start = time.time()
    
    parser = ArgumentParser()
    parser.add_argument("-f", "--filter")
    parser.add_argument("-T", "--telescope", default='all')
    parser.add_argument("-e", "--epoch")
    parser.add_argument("--field", nargs='+', default=['landolt', 'sloan'], choices=['landolt', 'sloan'])
    parser.add_argument("-t", "--timeout", type=float, default=86400.,
                        help="maximum time for each lscloop.py command in seconds")
    args = parser.parse_args()
    
    basecmd = 'lscloop.py'
    if args.epoch is not None:
        basecmd += ' -e ' + args.epoch
    if args.filter is not None:
        basecmd += ' -f ' + args.filter
    if args.telescope == 'all':
        basecmd_tel = basecmd
        telescopes = ['elp', 'coj', 'cpt', 'lsc', 'ogg', 'tfn']
    else:
        basecmd_tel = basecmd + ' -T ' + args.telescope
        telescopes = [args.telescope]
    
    logfile = open('lsc.log', 'w')
    #  compute astrometry, when missing, with astrometry.net
    run_cmd(basecmd_tel + ' -b wcs -s wcs --mode astrometry', logfile, args.timeout)
    #  try again with SV astrometry or set to bad image
    run_cmd(basecmd_tel + ' -b wcs -s wcs --xshift 1 --yshift 1', logfile, args.timeout)
    # compute  psf, when missing, with catalog
    run_cmd(basecmd_tel + ' -b psf -s psf', logfile, args.timeout)
    # compute  psf, when missing, with sextractor
    run_cmd(basecmd_tel + ' -b psf -s psf --use-sextractor', logfile, args.timeout)

    #############################################################################

    for tel in telescopes:
        logfile.write('####  starting reduction for ' + tel + '\n')
        logfile.flush()
        ll = lsc.myloopdef.get_list(args.epoch, tel)
        if not ll:
            continue
        lista = set(ll['objname'])
        standard = []
        for obj in lista:
            img = ll['filepath'][ll['objname'] == obj][0] + ll['filename'][ll['objname'] == obj][0]
            _, _, objtype = lsc.util.checksndb(img)
            if objtype == 1:
                standard.append(obj)
        logfile.write('####  standard fields: ' + str(standard) + '\n')
        logfile.flush()

        for obj in lista:
            for field in args.field:
                # compute zeropoint
                if lsc.util.getcatalog(obj, field):
                    run_cmd(basecmd + ' -b zcat -s zcat -F --cutmag 6 --field ' + field + ' -n ' + obj + ' -T ' + tel, logfile, args.timeout)
                # produce catalogs
                run_cmd(basecmd + ' -b abscat -s abscat -F -n ' + obj + ' --field ' + field + ' -T ' + tel + ' --standard all', logfile, args.timeout)
            # compute zeropoint with apass if not in landolt/sloan
            run_cmd(basecmd + ' -b zcat -s zcat -F --cutmag 6 --field apass -n ' + obj + ' -T ' + tel, logfile, args.timeout)
            if obj not in standard:
                # compute psf & apparent magnitudes
                run_cmd(basecmd_tel + ' -b psfmag -s psfmag -n ' + obj + commandsn.get(obj, ''), logfile, args.timeout)
                run_cmd(basecmd + ' -b mag -s mag -n ' + obj + ' -T ' + tel, logfile, args.timeout)
    
    epoch_to_fpack = datetime.strftime(datetime.utcnow() - timedelta(31), '%Y%m%d')
    run_cmd(basecmd + ' -e {} -s fpack'.format(epoch_to_fpack), logfile, args.timeout)
    run_cmd(basecmd + ' -e {} -s fpack -b quality'.format(epoch_to_fpack), logfile, args.timeout)

    stop = time.time()
    logfile.write('####  time to process all data: ' + str(stop - start) + ' seconds\n')
    logfile.close()
