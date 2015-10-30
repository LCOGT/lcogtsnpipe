from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
from os import sys, path
import os,shutil,re
from glob import glob
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']


from imp import find_module
try: find_module('numpy')
except: sys.exit('### Error: python module numpy not found')
    
try: find_module('pyfits')
except: sys.exit('### Error: python module pyfits not found')

#try: find_module('pyraf')
#except: sys.exit('### Error: python module pyraf not found')

try: find_module('matplotlib')
except: sys.exit('### Error: python module matplotlib not found')

#try: find_module('MySQLdb')
#except: sys.exit('### Error: python module MySQLdb not found')


verstr = "unknown"
try:
    parentdir=os.getcwd()+'/'
    verstrline = open(parentdir+'/src/lsc/_version.py', "rt").read()
except EnvironmentError:
    pass # Okay, there is no version file.
else:
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        raise RuntimeError("unable to find version in "+parentdir+"+src/floyds/_version.py")


setup(
    name='lsc',
    version=verstr,#'0.1.3',
    author='S. Valenti',
    author_email='svalenti@lcogt.net',
    scripts=['bin/lscingestrawdata.py','bin/lscingestredudata.py','bin/lscabsphot.py','bin/runlsc.py','bin/lscdiff.py',\
                 'bin/lscastro.py','bin/lsccatalogue.py','bin/lsccheck.py','bin/lscloop.py','bin/lscmerge.py',\
                 'bin/lscmaketempl.py','bin/runingestion.py', 'bin/lscmaglocal.py', 'bin/lscingestsloan.py',\
                 'bin/lscmag.py','bin/lscpsf.py','bin/lscsn.py','bin/lsctestheader.py','bin/lscnewcalib.py'],
    url='lcogt.net',
    license='LICENSE.txt', 
    description='lsc is a package to reduce PHTOMETRIC SN data',
    long_description=open('README.txt').read(),
    requires=['numpy','pyfits','matplotlib','MySQLdb'],
    packages=['lsc'],
    package_dir={'':'src'},
    package_data = {'lsc' : ["standard/astrometry/*cat","standard/*txt","standard/stdlist/*txt",\
                                 "standard/cat/*dat","standard/cat/*cat",\
                                 "standard/cat/sloan/*cat","standard/cat/landolt/*cat",\
                                 "standard/cat/sloan/*cat","standard/cat/apass/*cat",\
                                 "standard/cat/sloanprime/*cat","standard/cat/sloannatural/*cat","standard/cat/landoltnatural/*cat",\
                                 "standard/sex/*"]}
)
