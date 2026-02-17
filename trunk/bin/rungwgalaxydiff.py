#!/usr/bin/env python

from runlsc import run_cmd
from lsc.mysqldef import query
from lsc import conn

if __name__ == "__main__":
    logfile = open('gwgalaxydiff.log', 'w+')
    galaxies_to_subtract = query(['select targetid, SDSS_g, SDSS_r, SDSS_i, PS1_g, PS1_r, PS1_i, DECam_g, DECam_r, DECam_i, LCO_g, LCO_r, LCO_i from o4_galaxies'], conn)
    if galaxies_to_subtract:
        for galaxy in galaxies_to_subtract:
            targetid = int(galaxy['targetid'])
            images_to_subtract = query(['select distinct filter from photlco where targetid={} and filetype=1'.format(targetid)], conn)
            for filt_dict in images_to_subtract:
                filt = filt_dict['filter'][0]
                run_cmd('lscloop.py --targetid {} -f {} -e 20200101-20300101 --filetype 1 -s cosmic'.format(targetid, filt), logfile, 86400.0)
                if galaxy['SDSS_'+filt]:
                    run_cmd('lscloop.py --targetid {} -f {} -e 20200101-20300101 -s diff --temptel SDSS --tempdate 19990101-20191231 --normalize t -T 1m0'.format(targetid, filt), logfile, 86400.0)
                elif galaxy['PS1_'+filt]:
                    run_cmd('lscloop.py --targetid {} -f {} -e 20200101-20300101 -s diff --temptel PS1 --tempdate 19990101-20191231 --normalize t -T 1m0'.format(targetid, filt), logfile, 86400.0)
                elif galaxy['DECam_'+filt]:
                    run_cmd('lscloop.py --targetid {} -f {} -e 20200101-20300101 -s diff --temptel decam --tempdate 19990101-20191231 --normalize t -T 1m0'.format(targetid, filt), logfile, 86400.0)
                elif galaxy['LCO_'+filt]:
                    run_cmd('lscloop.py --targetid {} -f {} -e 20200101-20300101 -s diff --tempdate 19990101-20191231 --normalize t -T 1m0'.format(targetid, filt), logfile, 86400.0)
    


