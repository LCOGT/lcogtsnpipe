from mysqldef import *
from util import *
from lscastrodef import *
from lscabsphotdef import *
from lscsnoopy import *
from sqlcl import *
from sites import *
from externaldata import *
from myloopdef import *
from cosmics import *

__version__ = "unknown"
try:
    from _version import __version__
except ImportError:
    # We're running in a tree that doesn't have a _version.py, so we don't know what our version is.
    pass
