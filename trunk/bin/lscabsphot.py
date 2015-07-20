#!/usr/bin/env python
description = ">> zeropoint from table"
usage = "%prog image [options] "

import os
import string
import re
import sys
from optparse import OptionParser
import time
import lsc

# #######################################################################

if __name__ == "__main__":
    start_time = time.time()
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-i", "--interactive", action="store_true", dest='interactive', default=False,
                      help='Interactive \t\t\t [%default]')
    parser.add_option("-f", "--fix", action="store_true", dest='fix', default=False,
                      help='fix color \t\t\t [%default]')
    parser.add_option("-R", "--rejection", dest="rejection", default=2,
                      type='float', help='rejection \t\t\t %default')
    parser.add_option("-C", "--color", dest="color", default='',
                      type='str', help='colore slope [BV, VR, RI, gr, ri, iz] \t\t %default')
    parser.add_option("-s", "--system", dest="field", default='',
                      type='str', help='photometric system [sloan, landolt] \t\t %default')
    parser.add_option("-c", "--catalogue", dest="catalogue", default='',
                      type='str', help='catalogue to use (specific catalogue for calibration) \t\t %default')
    parser.add_option("-t", "--type", dest="type", default='fit',
                      type='str', help='type of magnitude (fit, ph) \t\t %default')
    parser.add_option("-r", "--redo", action="store_true", dest='redo', default=False,
                      help='redo zeropoint  \t\t\t [%default]')
    parser.add_option("--show", action="store_true", dest='show', default=False,
                      help='show plot in not interactive case  \t\t\t [%default]')
    parser.add_option("--cutmag", dest="cutmag", default=99., type="float",
                      help='--cutmag  [magnitude instrumental cut for zeropoint ]  \t [%default]')
    parser.add_option("--calib", dest="calibration", default='sloan',
                      type='str', help='calibration to  (sloan,sloanprime,natural) \t\t %default')
    parser.add_option('--zcatnew', action='store_true', help='use fitcol3 for the zero point and color term')

    option, args = parser.parse_args()
    if len(args) < 1: sys.argv.append('--help')
    _type = option.type
    if _type not in ['fit', 'ph']: sys.argv.append('--help')
    option, args = parser.parse_args()
    imglist = lsc.util.readlist(args[0])
    _field = option.field
    _catalogue = option.catalogue
    _interactive = option.interactive
    _fix = option.fix
    _rejection = option.rejection
    _color = option.color
    _redo = option.redo
    _show = option.show
    _cutmag = option.cutmag
    _calib = option.calibration
    zcatnew = option.zcatnew
    for img in imglist:
        lsc.absphot(img, _field, _catalogue, _fix, _color, _rejection, _interactive, _type, _redo, _show, _cutmag,
                    'photlco', _calib, zcatnew)
