#!/usr/bin/env python

from lsc.mysqldef import query
from lsc import conn
from runlsc import run_cmd


if __name__ == "__main__":
    ### Get list of galaxies with no templates ingested
    logfile = open('runo4templateingest.log', 'w+')
    undownloaded_galaxies = query(['select id, targetid from o4_galaxies where ingested = 0'], conn)
    if undownloaded_galaxies:
        for galaxy in undownloaded_galaxies:
            ### Run the template ingestion
            for filt in ['g', 'r', 'i']:
                run_cmd('ingesttemplates.py -f {} --targetid {}'.format(filt, galaxy['targetid']), logfile)
    
                ### Look for templates for this filter in each survey
                template_query = query(["select filepath, filename from photlco where targetid={} and filetype=4 and filter='{}p'".format(galaxy['targetid'], filt)], conn)
                if template_query:
                    ### Update the filepath and ingested columns in the o4_galaxies table
                    for temp in template_query:
                        if 'SDSS' in temp['filename']:
                            update_row = query(["update o4_galaxies set SDSS_{} = '{}', ingested = 1 where id = {}".format(filt, temp['filepath']+temp['filename'], int(galaxy['id']))], conn)
                        elif 'PS' in temp['filename']:
                            update_row = query(["update o4_galaxies set PS1_{} = '{}', ingested = 1 where id = {}".format(filt, temp['filepath']+temp['filename'], int(galaxy['id']))], conn)
                        elif 'dec' in temp['filename'].lower():
                            update_row = query(["update o4_galaxies set DECam_{} = '{}', ingested = 1 where id = {}".format(filt, temp['filepath']+temp['filename'], int(galaxy['id']))], conn)
                        else:
                            update_row = query(["update o4_galaxies set LCO_{} = '{}', ingested = 1 where id = {}".format(filt, temp['filepath']+temp['filename'], int(galaxy['id']))], conn)
                            
                                    
