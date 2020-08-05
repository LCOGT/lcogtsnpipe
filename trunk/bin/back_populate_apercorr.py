#!/usr/bin/env python
"""
The aperture correction used to be stored in the header of the SN2 files and not in the database
This script reads the aperture correction from the SN2 files and uses it to populate the new
apercorr column in the photlco table. This should be run once on existing databases to populate 
already processed data. Future reductions will populate this keyword automatically
"""
import os
from astropy.io import fits
import lsc

# Query the database to get the filename and filepath
epoch = '20100101-20210101' #Make this big enough that it really gets all of the data
## Standard reduction
standard_file_dict = lsc.myloopdef.get_list(epoch=epoch, _telescope='all', _filter='', _bad='', _name='', _id='', _ra='', _dec='', database='photlco',
             filetype=1, _groupid=None, _instrument='', _temptel='', _difftype='', classid=None)
merged_file_dict = lsc.myloopdef.get_list(epoch=epoch, _telescope='all', _filter='', _bad='', _name='', _id='', _ra='', _dec='', database='photlco',
             filetype=2, _groupid=None, _instrument='', _temptel='', _difftype='', classid=None)
diff_file_dict = lsc.myloopdef.get_list(epoch=epoch, _telescope='all', _filter='', _bad='', _name='', _id='', _ra='', _dec='', database='photlco',
             filetype=3, _groupid=None, _instrument='', _temptel='', _difftype='', classid=None)
template_file_dict = lsc.myloopdef.get_list(epoch=epoch, _telescope='all', _filter='', _bad='', _name='', _id='', _ra='', _dec='', database='photlco',
             filetype=4, _groupid=None, _instrument='', _temptel='', _difftype='', classid=None)

# Modify the filename to be the sn2 filename
for file_dict in [standard_file_dict, merged_file_dict, diff_file_dict, template_file_dict]:
    for filename, filepath in zip(file_dict['filename'], file_dict['filepath']):
        img_filename = os.path.join(filepath, filename)
        sn2_filename = os.path.join(filepath, filename.replace('.fits', '.sn2.fits'))
        # Read the aperture correction from the header
        if os.path.exists(sn2_filename):
            apercorr = fits.getval(sn2_filename, 'APCO', 0)
            # Update the database with the aperture correction
            lsc.mysqldef.updatevalue('photlco', 'apercorr', apercorr, filename)
        #Notify user if PSF stage was run but no sn2 file is found
        elif os.path.exists(os.path.join(filepath, filename.replace('.fits', '.psf.fits'))):
            print(sn2_filename)

