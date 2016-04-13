###########################################################################
# 
#                                   lsc Pipeline
#
#                              INSTALLATION
###########################################################################


NTT is written in python and requires the following package:

- Python 2.5 or Python 2.6 or Python 2.7
   these modules have to be installed:
            - numpy
	    - pyraf
	    - matplotlib
	    - pyfits
	    - mysqldb
	    - astropy

- Iraf
	    
##############################################################################
extract the files from the tarball
> tar -xvf lsc-version.tar

> cd lsc-version
> python setup.py install  (--record files.txt) 

##########################################################################
To uninstall a previus version 

- delete the lsc directory in your site-package path
- delete the lsc****.egg-info from the same directory
- delete the lsc executable: lsc 

or if during installation  you used the option: --record files.txt
you can run the following command in theterminal:
> cat files.txt | xargs sudo rm -rf
