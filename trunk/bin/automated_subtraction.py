#!/usr/bin/env python


"""
====================================

* **Filename**:         run_subtraction.py 
* **Author**:              Joseph Farah 
* **Description**:       Automates subtraction on a schedule.

====================================

**Notes**
*  Observation created on target with ID at MJD. Cronjob runs. 
*  Two things to check: 
*    - does template exist for object in filter?
*    - has observation in filter been conducted after MJD?
* If yes: get filter, mjd,  filetag. run subtraction with that info


Steps:
- 


"""

#------------- IMPORTS -------------#
import lsc
import datetime
from subprocess import call
from mysql.connector import connect


#------------- SETTINGS -------------#
template_order_labels = ['PS1', 'LCO', 'DECam', 'Skymapper', 'SDSS']
template_order = []
for t in template_order_labels:
    for f in ['g', 'r', 'i']:
        template_order.append(t+'_'+f)

# INGESTION_MODE = True will download a new set of templates
# INGESTION_MODE = False will try to use Giacomo's templates
INGESTION_MODE = True


#------------- CLASSES -------------#
class PotentialTarget(object):
    '''
    Object to organize data from table. Will be used to check if template
    exists for unsubtracted observation. Bundles information and passes to 
    subtraction.
    '''

    def __init__(self, unpacked_attr):
        self.targetid = unpacked_attr[0]
        self.objname = unpacked_attr[1]
        self.filepath = unpacked_attr[2]
        self.filename = unpacked_attr[3]
        self.filetype = unpacked_attr[4]
        self.dayobs = unpacked_attr[5]
        self.filter = unpacked_attr[6]
        self.PS1_g = unpacked_attr[7]
        self.PS1_r = unpacked_attr[8]
        self.PS1_i = unpacked_attr[9]
        self.DECam_g = unpacked_attr[10]
        self.DECam_r = unpacked_attr[11]
        self.DECam_i = unpacked_attr[12]
        self.LCO_g = unpacked_attr[13]
        self.LCO_r = unpacked_attr[14]
        self.LCO_i = unpacked_attr[15]
        self.Skymapper_g = unpacked_attr[16]
        self.Skymapper_r = unpacked_attr[17]
        self.Skymapper_i = unpacked_attr[18]
        self.instrument = unpacked_attr[19]
        self.TEMP_SRC = None
        self.TEMP_FILT = None




#------------- FUNCTIONS -------------#
def send_query(query, conn):
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    return result


def get_target_list():
    """
        Checks all recent observations for those relevant to GW events
        and organizes their information in preparation for subtraction. 
        Designed to run on schedule.
    
        **Args**:
    
        * none (none): none
    
        **Returns**:
    
        * target_list (list): list of targets to perform subtraction on
    
    """
    

    ## minimum to date to examine is 3 days ago ##
    # d = datetime.datetime(2023, 9, 01) - datetime.timedelta(3) ## DEBUGGING
    d = datetime.date.today() - datetime.timedelta(3)

    ## convert to pipeline-friendly fmt (YYYYMMDD) ##
    d = d.strftime("%Y%m%d")

    '''
     grab ID, name, obs fpath/fname, filetype, 
     obsdate, and filter from obs in both 
     photlco AND o5 as well as any template info 
     for associated tampltes
    '''
    query  = "SELECT p.targetid, p.objname, p.filepath, p.filename, p.filetype, p.dayobs, p.filter, "
    query += "o.PS1_g, o.PS1_r, o.PS1_i, "
    query += "o.DECam_g, o.DECam_r, o.DECam_i, "
    query += "o.LCO_g, o.LCO_r, o.LCO_i, "
    query += "o.Skymapper_g, o.Skymapper_r, o.Skymapper_i, p.instrument "
    query += "FROM photlco p INNER JOIN o4_galaxies o ON p.targetid = o.targetid " ## ensures that only targets in both are considered
    query += "WHERE p.dayobs >= {0};".format(d)

    ## connect to database ##
    conn = connect(
        user = 'supernova',
        password = 'supernova',
        host = 'supernovadb',
        database = 'supernova'
    )
 
    result = send_query(query, conn)

    ##
    potential_targets = [PotentialTarget(r) for r in result]

    target_list = []
    for target in potential_targets:
        for temp in template_order:
            if getattr(target, temp) is None or getattr(target, temp) == '--':
                continue
            else:
                target.TEMP_SRC = temp.split('_')[0]
                target.TEMP_FILT = temp.split('_')[1]
                if target.TEMP_FILT == target.filter or target.TEMP_FILT+'p' == target.filter:
                    print("Template found with properties: {0},{1}".format(target.TEMP_SRC, target.TEMP_FILT))
                    target_list.append(target)
                    break 
    

    print("List of targets: ", target_list)
    return target_list





def run_subtraction(targets, ingestion_mode=False):
    """
        Runs subtraction on narrowed list of targets provided by
        get_target_list().

        Info needed for subtraction: 
         - instrument
         - target id
         - filter
         - date
         - template date
    
        **Args**:
    
        * targets (list): list of targets to perform subtraction on
        * ingestion_mode (bool): whether to perform on-the-fly template download
    
        **Returns**:
    
        * none (none): none
    
    """


    def subtract_LCO(target):
        print("Perfoming LCO subtraction on {0}.".format(target.objname))
        print("No commands available.")


    def subtract_DECam(target):
        print("Perfoming DEcam subtraction on {0}.".format(target.objname))
        print("No commands available.")

    def subtract_PS1(target, ingestion_mode=False):
        print("Perfoming PS1 subtraction on {0}.".format(target.objname))
        objname = target.objname
        date = target.dayobs
        target_id = target.targetid
        conn = connect(
            user = 'supernova',
            password = 'supernova',
            host = 'supernovadb',
            database = 'supernova'
        )
        filename = getattr(target, target.TEMP_SRC+'_'+target.TEMP_FILT).split('/')[-1]
        if ingestion_mode==True:
            tempdate = '20100101-20160101'
            filt = target.filter.replace('p', '')
            instrument = filter(str.isalpha, target.instrument.encode('utf-8'))
            command = "lscloop.py -n {0} -e {1} && lscloop.py -n {0} -e {1} -s ingestps1 -f {4} && lscloop.py -n {0} -e {3} --filetype 4 -s psf --fwhm 5 --use-sextractor && lscloop.py -n {0} -e {3} --filetype 4 -s cosmic && lscloop.py -n {0} -e {1} -s cosmic && lscloop.py -n {0} -e {1} --normalize i --convolve t -T {5} --tempdate {3} --temptel PS1 --fixpix --difftype 0 -s diff -F".format(objname, date, target_id, tempdate, filt, instrument)
            print(command)
            # input()

        else:
            query = "SELECT dayobs FROM photlco WHERE targetid={1} AND filename='{0}'".format(filename, target_id)
            print(query)
            tempdate = send_query(query, conn)[0][0].encode('utf-8')
            filt = target.filter.replace('p', '')
            instrument = filter(str.isalpha, target.instrument.encode('utf-8'))
            command = "lscloop.py -n {0} -e {1} && lscloop.py -n {0} -e {3} --filetype 4 -s psf --fwhm 5 --use-sextractor && lscloop.py -n {0} -e {3} --filetype 4 -s cosmic && lscloop.py -n {0} -e {1} -s cosmic && lscloop.py -n {0} -e {1} --normalize i --convolve t -T {5} --tempdate {3} --temptel PS1 --fixpix --difftype 0 -s diff -F".format(objname, date, target_id, tempdate, filt, instrument)

        # command = "lscloop.py -n {0} -e {1} && lscloop.py -n {0} -e {1} -d {2} -s ingestps1 -f {4} && lscloop.py -n {0} -e {3} --filetype 4 -s psf --fwhm 5 --use-sextractor && lscloop.py -n {0} -e {3} --filetype 4 -s cosmic && lscloop.py -n {0} -e {1} -s cosmic && lscloop.py -n {0} -e {1} --normalize i --convolve t -T {5} --tempdate {3} --temptel PS1 --fixpix --difftype 0 -s diff --no_iraf".format(objname, date, tel_id, tempdate, filt, instrument)


        call(command, shell=True)



    def subtract_Skymapper(target):
        print("Perfoming Skymapper subtraction on {0}.".format(target.objname))
        print("No commands available.")


    def subtract_SDSS(target):
        print("Perfoming SDSS subtraction on {0}.".format(target.objname))
        print("No commands available.")

        objname = target.objname
        date = target.dayobs
        target_id = target.targetid
        conn = connect(
            user = 'supernova',
            password = 'supernova',
            host = 'supernovadb',
            database = 'supernova'
        )
        filename = getattr(target, target.TEMP_SRC+'_'+target.TEMP_FILT).split('/')[-1]
        if ingestion_mode==True:
            tempdate = '20060101-20070101'
            filt = target.filter.replace('p', '')
            instrument = filter(str.isalpha, target.instrument.encode('utf-8'))
            command = "lscloop.py -n {0} -e {1} && lscloop.py -n {0} -e {1} -s ingestsdss -f {4} && lscloop.py -n {0} -e {3} --filetype 4 -s psf --fwhm 5 --use-sextractor && lscloop.py -n {0} -e {3} --filetype 4 -s cosmic && lscloop.py -n {0} -e {1} -s cosmic && lscloop.py -n {0} -e {1} --normalize i --convolve t -T {5} --tempdate {3} --temptel PS1 --fixpix --difftype 0 -s diff -F".format(objname, date, target_id, tempdate, filt, instrument)
            print(command)
            # input()


        # command = "lscloop.py -n {0} -e {1} && lscloop.py -n {0} -e {1} -d {2} -s ingestps1 -f {4} && lscloop.py -n {0} -e {3} --filetype 4 -s psf --fwhm 5 --use-sextractor && lscloop.py -n {0} -e {3} --filetype 4 -s cosmic && lscloop.py -n {0} -e {1} -s cosmic && lscloop.py -n {0} -e {1} --normalize i --convolve t -T {5} --tempdate {3} --temptel PS1 --fixpix --difftype 0 -s diff --no_iraf".format(objname, date, tel_id, tempdate, filt, instrument)


        call(command, shell=True)



    subfunc_switchboard = {
        'LCO':subtract_LCO,
        'DECam':subtract_DECam,
        'PS1':subtract_PS1,
        'Skymapper':subtract_Skymapper,
        'SDSS':subtract_SDSS,
    }
    
    print("Running subtraction on: ", targets)
    
    ## iterate over targets ##
    for target in targets:

        print("Loading target {0} with target ID {1} and object name {2}".format(target, target.targetid, target.objname))

        ## ingest template with filetype=4
        print(target.TEMP_SRC)
        print(target.TEMP_FILT)
        lsc.mysqldef.ingestredu([getattr(target, target.TEMP_SRC+'_'+target.TEMP_FILT)], 'no', filetype=4)

        ## associate with correct subtraction function ##
        subtraction_func = subfunc_switchboard[target.TEMP_SRC]

        ## perform subtraction ##
        subtraction_func(target, ingestion_mode=ingestion_mode)





    



def main():
    global INGESTION_MODE

    target_list = get_target_list()
    run_subtraction(target_list, ingestion_mode=INGESTION_MODE)


#------------- SWITCHBOARD -------------#
if __name__ == '__main__':
    main()

