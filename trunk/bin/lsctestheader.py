#!/usr/bin/env python

import lsc
import glob
import sys

if len(sys.argv) <= 1:
    listfits = glob.glob('*fits')
    for img in listfits:
        print img
    img = raw_input('Which image do you want to test [' + str(listfits[0]) + '] ? ')
    if not img: img = listfits[0]
else:
    img = sys.argv[1]

hdr = lsc.util.readhdr(img)

_imagetype = lsc.util.readkey3(hdr, 'type')
_object = lsc.util.readkey3(hdr, 'object')
_mjd = lsc.util.readkey3(hdr, 'mjd')
_airmass = lsc.util.readkey3(hdr, 'airmass')
_filter = lsc.util.readkey3(hdr, 'filter')
_grism = lsc.util.readkey3(hdr, 'grism')
_exptime = lsc.util.readkey3(hdr, 'exptime')
_date = lsc.util.readkey3(hdr, 'date-obs')
_gain = lsc.util.readkey3(hdr, 'gain')
_ron = lsc.util.readkey3(hdr, 'ron')
_lampid = lsc.util.readkey3(hdr, 'lampid')
_RA = lsc.util.readkey3(hdr, 'RA')
_DEC = lsc.util.readkey3(hdr, 'DEC')
_ccdmax = lsc.util.readkey3(hdr, 'datamax')
_ccdmin = lsc.util.readkey3(hdr, 'datamin')
_cenwav = lsc.util.readkey3(hdr, 'cenw')
_slitw = lsc.util.readkey3(hdr, 'slit')
_UT = lsc.util.readkey3(hdr, 'ut')
_xdimen = lsc.util.readkey3(hdr, 'NAXIS1')
_ydimen = lsc.util.readkey3(hdr, 'NAXIS2')
_instrument = lsc.util.readkey3(hdr, 'instrume')
_obsmode = lsc.util.readkey3(hdr, 'obsmode')

if not _gain:
    _gain = '########'
if not _ron:
    _ron = '########'
if not _instrument:
    _instrument = '########'
if not _ydimen:
    _ydimen = '########'
if not _xdimen:
    _xdimen = '########'
if not _filter:
    _filter = '########'
if not _RA:
    _RA = '########'
if not _grism:
    _grism = '########'
if not _slitw:
    _slitw = '########'
if not _lampid:
    _lampid = '#######'
if not _date:
    _date = '#######'
if not _cenwav:
    _cenwav = '#######'
if not _UT:
    _UT = '#######'
if not _ccdmin:
    _ccdmin = '#######'
if not _ccdmax:
    _ccdmax = '#######'
if not _obsmode:
    _obsmode = '#######'
if not _object:
    _object = '#######'
_system = '#######'

print '####################################################################'
print 'IMG                OBJECT  IMAGETYPE    EXPTIME    FILTER        GRISM      '
print str(img) + '\t' + str(_object) + '\t' + str(_imagetype) + '\t' + str(_exptime) + '\t' + str(_filter) + '\t' + str(
    _grism)
print '####################################################################'
print 'AIRMASS             MJD             DATE          XDIM   YDIM    GAIN   RON '
print str(_airmass) + '\t' + str(_mjd) + '\t' + str(_date) + '\t' + str(_xdimen) + '\t' + str(_ydimen) + '\t' + str(
    _gain) + '\t' + str(_ron)
print '####################################################################'
print 'LAMP_ID     slitw         RA          DEC     CCDMIN    CCDMAX   CENWAV '
print str(_lampid) + '\t' + str(_slitw) + '\t' + str(_RA) + '\t' + str(_DEC) + '\t' + str(_ccdmin) + '\t' + str(
    _ccdmax) + '\t' + str(_cenwav)
print '####################################################################'
print ' UT     xdimension    ydimension        instrument      SYSTEM   OBSMODE '
print str(_UT) + '\t' + str(_xdimen) + '\t' + str(_ydimen) + '\t' + str(_instrument) + '\t' + str(_system) + '\t' + str(
    _obsmode)
print '####################################################################'

