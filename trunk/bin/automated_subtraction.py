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


