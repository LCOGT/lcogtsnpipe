import lsc.util
#from util import *
import lsc.mysqldef
#from mysqldef import *
import lsc.lscastrodef
#from lscastrodef import *
import lsc.lscabsphotdef
#from lscabsphotdef import *
import lsc.lscsnoopy
#from lscsnoopy import *
import lsc.sqlcl
#from sqlcl import *
import lsc.sites
#from sites import *
import lsc.externaldata
#from externaldata import *
import lsc.myloopdef
#from myloopdef import *
import lsc.cosmics
#from cosmics import *
import lsc.lscpsfdef
#from lscpsfdef import *
import lsc.banzaicat
#from banzaicat import *

__version__ = "unknown"
try:
    from _version import __version__
except ImportError:
    # We're running in a tree that doesn't have a _version.py, so we don't know what our version is.
    pass
