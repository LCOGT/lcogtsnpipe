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

