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
template_order_labels = ['PS1', 'LCO', 'DECam', 'Skymapper']
template_order = []
for t in template_order_labels:
    for f in ['g', 'r', 'i']:
        template_order.append(t+'_'+f)





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
    d = datetime.datetime(2023, 9, 01) - datetime.timedelta(3) ## DEBUGGING
    # d = datetime.date.today() - datetime.timedelta(3)

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
