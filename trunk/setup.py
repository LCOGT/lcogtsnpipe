# from distutils.core import setup
from setuptools import setup
from os import path
from glob import glob

setup(
    name='lsc',
    version='1.0.0',
    author='S. Valenti',
    author_email='svalenti@lcogt.net',
    scripts=glob(path.join('bin', '*.py')),
    url='lcogt.net',
    license='LICENSE.txt', 
    description='lsc is a package to reduce PHTOMETRIC SN data',
    long_description=open('README.txt').read(),
    requires=['numpy','astropy','matplotlib','MySQLdb', 'pyraf', 'reproject', 'mysql-connector-python-rf'],
    packages=['lsc'],
    package_dir={'': 'src'},
    package_data={'lsc': ["standard/astrometry/*cat","standard/*txt","standard/stdlist/*txt",
                         "standard/cat/*dat","standard/cat/*cat",
                         "standard/cat/sloan/*cat","standard/cat/landolt/*cat","standard/cat/apass/*cat",
                         "standard/cat/sloanprime/*cat","standard/cat/sloannatural/*cat","standard/cat/landoltnatural/*cat",
                         "standard/sex/*", 'configure']}
)
