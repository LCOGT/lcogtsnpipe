"""
Shared pytest fixtures for the lcogtsnpipe test suite.

Critical note: pyraf (and the iraf binary) is NOT required at test time.
We stub the entire pyraf namespace in sys.modules before any lsc import so
that module-level statements like ``from pyraf import iraf`` and
``iraf.noao(_doprint=0)`` are silently absorbed by MagicMock objects.
"""
import sys
import types
from unittest.mock import MagicMock

# ---------------------------------------------------------------------------
# Stub pyraf / iraf BEFORE any lsc module is imported
# ---------------------------------------------------------------------------
_iraf_mock = MagicMock()
_pyraf_mock = types.ModuleType("pyraf")
_pyraf_mock.iraf = _iraf_mock
sys.modules.setdefault("pyraf", _pyraf_mock)
sys.modules.setdefault("pyraf.iraf", _iraf_mock)
# iraf itself is sometimes imported as a bare name inside functions
sys.modules.setdefault("iraf", _iraf_mock)

# ---------------------------------------------------------------------------
# Also stub MySQLdb / pymysql so mysqldef imports cleanly without a DB server
# ---------------------------------------------------------------------------
_msql_mock = MagicMock()
_cursor_mock = MagicMock()
_cursor_mock.fetchall.return_value = ()
_cursor_mock.rowcount = 0
_msql_mock.connect.return_value.cursor.return_value = _cursor_mock
_msql_mock.cursors = MagicMock()
_msql_mock.cursors.DictCursor = MagicMock()
sys.modules.setdefault("MySQLdb", _msql_mock)
sys.modules.setdefault("pymysql", _msql_mock)

import os
import datetime
import numpy as np
import pytest
from astropy.io import fits


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def simple_fits(tmp_path_factory):
    """A minimal 100×100 FITS file with standard header keywords."""
    tmp = tmp_path_factory.mktemp("fits")
    path = tmp / "test.fits"
    data = np.random.default_rng(42).normal(1000, 50, (100, 100)).astype(np.float32)
    hdr = fits.Header()
    hdr["SIMPLE"] = True
    hdr["BITPIX"] = -32
    hdr["NAXIS"] = 2
    hdr["NAXIS1"] = 100
    hdr["NAXIS2"] = 100
    hdr["EXPTIME"] = 120.0
    hdr["GAIN"] = 2.0
    hdr["RDNOISE"] = 10.0
    hdr["PIXSCALE"] = 0.389
    hdr["RA"] = 150.0
    hdr["DEC"] = 2.2
    hdr["AIRMASS"] = 1.2
    hdr["FILTER"] = "r"
    hdr["OBJECT"] = "SN2024abc"
    hdr["TELESCOP"] = "1m0-01"
    hdr["SITEID"] = "lsc"
    hdr["MJD"] = 58970.5
    hdr["MJD-OBS"] = 58970.5
    hdr["DATE-OBS"] = "2020-05-01T12:00:00"
    hdr["INSTRUME"] = "fa15"
    hdr["NAXIS1"] = 100
    hdr["NAXIS2"] = 100
    hdr["CCDSUM"] = "1 1"
    hdr["ROLLERDR"] = 0.0
    hdr["WCSERR"] = 0
    fits.writeto(str(path), data, hdr, overwrite=True)
    return str(path)


@pytest.fixture(scope="session")
def banzai_fits(tmp_path_factory):
    """A FITS file with a synthetic BANZAI-style CAT binary table extension."""
    tmp = tmp_path_factory.mktemp("banzai")
    path = tmp / "banzai.fits"

    # Primary extension (science image)
    data = np.ones((50, 50), dtype=np.float32) * 1000.0
    primary = fits.PrimaryHDU(data)

    # CAT extension - matching BANZAI column names
    n = 20
    rng = np.random.default_rng(0)
    col_ra   = fits.Column(name="RA",          format="D", array=rng.uniform(149, 151, n))
    col_dec  = fits.Column(name="DEC",         format="D", array=rng.uniform(1.5, 2.5, n))
    col_ell  = fits.Column(name="ELLIPTICITY", format="E", array=rng.uniform(0.0, 0.3, n))
    col_bg   = fits.Column(name="BACKGROUND",  format="E", array=rng.uniform(800, 1200, n))
    col_fwhm = fits.Column(name="FWHM",        format="E", array=rng.uniform(3.0, 8.0, n))
    col_peak = fits.Column(name="PEAK",        format="E", array=rng.uniform(2000, 60000, n))
    col_flag = fits.Column(name="FLAG",        format="I", array=np.zeros(n, dtype=np.int16))
    cat_hdu = fits.BinTableHDU.from_columns(
        [col_ra, col_dec, col_ell, col_bg, col_fwhm, col_peak, col_flag],
        name="CAT"
    )
    hdul = fits.HDUList([primary, cat_hdu])
    hdul.writeto(str(path), overwrite=True)
    return str(path)


@pytest.fixture(scope="session")
def image_array():
    """A synthetic 100×100 image as a numpy array (sky + Gaussian stars + cosmics)."""
    rng = np.random.default_rng(7)
    arr = rng.normal(1000.0, 30.0, (100, 100)).astype(np.float32)
    # Add a few fake 'cosmic ray' hot pixels
    arr[30, 40] = 50000.0
    arr[70, 20] = 45000.0
    arr[10, 90] = 48000.0
    return arr


@pytest.fixture
def mock_db_conn():
    """
    A MagicMock representing a MySQLdb/pymysql connection.
    Callers can customise cursor.fetchall.return_value for specific queries.
    """
    conn = MagicMock()
    cursor = MagicMock()
    cursor.rowcount = 1
    cursor.fetchall.return_value = ({"filename": "test.fits"},)
    conn.cursor.return_value = cursor
    # DictCursor class attribute
    conn.cursor.return_value.__class__ = MagicMock
    return conn
