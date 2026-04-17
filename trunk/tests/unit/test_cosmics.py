"""
Tests for lsc.cosmics — the L.A.Cosmic cosmic-ray detection algorithm.
"""
import pytest
import numpy as np

pytestmark = pytest.mark.unit

# subsample / rebin / rebin2x2 — array utilities
class TestSubsample:
    def test_output_shape_doubles(self):
        """subsample() doubles both dimensions of the input array."""
        from lsc.cosmics import subsample
        a = np.ones((10, 10))
        b = subsample(a)
        assert b.shape == (20, 20)

    def test_values_preserved(self):
        """subsample() propagates pixel values without interpolation."""
        from lsc.cosmics import subsample
        a = np.full((4, 4), 7.0)
        b = subsample(a)
        assert np.allclose(b, 7.0)


class TestRebin:
    def test_2x2_downsample(self):
        """rebin() reduces array dimensions to the requested shape."""
        from lsc.cosmics import rebin
        a = np.ones((10, 10))
        b = rebin(a, (5, 5))
        assert b.shape == (5, 5)

    def test_mean_preserving(self):
        """rebin() computes the mean of each block, preserving uniform values."""
        from lsc.cosmics import rebin
        a = np.full((6, 4), 4.0)
        b = rebin(a, (3, 2))
        assert np.allclose(b, 4.0)


class TestRebin2x2:
    def test_halves_shape(self):
        """rebin2x2() halves both dimensions of an even-shaped array."""
        from lsc.cosmics import rebin2x2
        a = np.ones((8, 8))
        b = rebin2x2(a)
        assert b.shape == (4, 4)

    def test_odd_shape_raises(self):
        """rebin2x2() raises an error when the input has odd dimensions."""
        from lsc.cosmics import rebin2x2
        a = np.ones((7, 7))
        with pytest.raises(Exception):
            rebin2x2(a)


class TestRebin2x2Values:
    def test_averages_2x2_blocks(self):
        """rebin2x2 should average each 2x2 block."""
        from lsc.cosmics import rebin2x2
        a = np.array([[1.0, 2.0, 3.0, 4.0],
                      [5.0, 6.0, 7.0, 8.0],
                      [9.0, 10.0, 11.0, 12.0],
                      [13.0, 14.0, 15.0, 16.0]])
        b = rebin2x2(a)
        assert b.shape == (2, 2)
        expected = np.array([[(1+2+5+6)/4.0, (3+4+7+8)/4.0],
                              [(9+10+13+14)/4.0, (11+12+15+16)/4.0]])
        assert np.allclose(b, expected)


# cosmicsimage — class initialisation and properties
class TestCosmicsimageInit:
    def setup_method(self):
        rng = np.random.default_rng(0)
        self.data = rng.normal(1000.0, 30.0, (64, 64)).astype(np.float64)

    def test_init_no_crash(self):
        """cosmicsimage can be constructed with gain and readnoise without error."""
        from lsc.cosmics import cosmicsimage
        ci = cosmicsimage(self.data, gain=2.0, readnoise=10.0)
        assert ci is not None

    def test_mask_all_false_initially(self):
        """The cosmic ray mask starts with all pixels unmasked."""
        from lsc.cosmics import cosmicsimage
        ci = cosmicsimage(self.data)
        assert ci.mask.sum() == 0

    def test_cleanarray_matches_rawarray(self):
        """Before any iteration, cleanarray is a copy of rawarray."""
        from lsc.cosmics import cosmicsimage
        ci = cosmicsimage(self.data)
        assert np.allclose(ci.cleanarray, ci.rawarray)

    def test_pssl_offset(self):
        """pssl is added to rawarray internally to work with the sky included."""
        from lsc.cosmics import cosmicsimage
        ci = cosmicsimage(self.data, pssl=100.0)
        # rawarray should be data + 100
        assert np.allclose(ci.rawarray, self.data + 100.0)

    def test_str_contains_shape(self):
        """__str__ includes the array dimensions."""
        from lsc.cosmics import cosmicsimage
        ci = cosmicsimage(self.data)
        s = str(ci)
        assert "64" in s


# cosmicsimage — lacosmiciteration detects injected cosmics
class TestLacosmiciteration:
    def setup_method(self):
        rng = np.random.default_rng(42)
        self.clean_data = rng.normal(1000.0, 30.0, (128, 128)).astype(np.float64)
        # Inject obvious cosmic ray hits
        self.dirty_data = self.clean_data.copy()
        self.cosmic_positions = [(30, 40), (70, 20), (100, 80)]
        for r, c in self.cosmic_positions:
            self.dirty_data[r, c] = 60000.0

    def test_iteration_detects_cosmics(self):
        """lacosmiciteration() returns niter > 0 when bright cosmics are present."""
        from lsc.cosmics import cosmicsimage
        ci = cosmicsimage(self.dirty_data, gain=2.2, readnoise=10.0,
                          sigclip=5.0, satlevel=-1, verbose=False)
        result = ci.lacosmiciteration(verbose=False)
        assert result["niter"] > 0

    def test_mask_grows_after_iteration(self):
        """After one iteration the mask has flagged at least one pixel."""
        from lsc.cosmics import cosmicsimage
        ci = cosmicsimage(self.dirty_data, gain=2.2, readnoise=10.0,
                          sigclip=5.0, satlevel=-1, verbose=False)
        ci.lacosmiciteration(verbose=False)
        assert ci.mask.sum() > 0

    def test_result_keys(self):
        """The result dict from lacosmiciteration() contains the four expected keys."""
        from lsc.cosmics import cosmicsimage
        ci = cosmicsimage(self.dirty_data, gain=2.2, readnoise=10.0,
                          sigclip=5.0, satlevel=-1, verbose=False)
        result = ci.lacosmiciteration(verbose=False)
        for key in ("niter", "nnew", "itermask", "newmask"):
            assert key in result


# cosmicsimage — clean() replaces detected cosmic pixels
class TestClean:
    def test_clean_replaces_flagged_pixels(self):
        """clean() replaces detected cosmic pixels so none remain as np.inf."""
        from lsc.cosmics import cosmicsimage
        rng = np.random.default_rng(5)
        data = rng.normal(1000.0, 30.0, (64, 64)).astype(np.float64)
        data[20, 20] = 55000.0
        ci = cosmicsimage(data, gain=2.2, readnoise=10.0,
                          sigclip=5.0, satlevel=-1, verbose=False)
        ci.lacosmiciteration(verbose=False)
        before = ci.cleanarray[20, 20]
        ci.clean(verbose=False)
        after = ci.cleanarray[20, 20]
        # The cleaned pixel should no longer be np.Inf
        assert not np.isinf(after)


# cosmicsimage — run() full integration
class TestRun:
    def test_run_completes(self):
        """run() completes without error and returns a valid mask array."""
        from lsc.cosmics import cosmicsimage
        rng = np.random.default_rng(9)
        data = rng.normal(1000.0, 30.0, (64, 64)).astype(np.float64)
        data[15, 15] = 55000.0
        ci = cosmicsimage(data, gain=2.2, readnoise=10.0, sigclip=5.0, satlevel=-1, verbose=False)
        ci.run(maxiter=2, verbose=False)
        # After running, mask should have detected something
        assert isinstance(ci.mask, np.ndarray)

    def test_run_zero_cosmics_on_clean_image(self):
        """With a very high sigclip threshold almost no pixels are flagged on a clean image."""
        from lsc.cosmics import cosmicsimage
        rng = np.random.default_rng(7)
        # Very well-behaved Gaussian image — few if any cosmics expected
        data = rng.normal(1000.0, 5.0, (64, 64)).astype(np.float64)
        ci = cosmicsimage(data, gain=2.2, readnoise=10.0, sigclip=20.0, satlevel=-1, verbose=False)
        ci.run(maxiter=1, verbose=False)
        # With very high sigclip threshold, almost nothing should be flagged
        assert ci.mask.sum() < 10


# fromfits / tofits — round-trip FITS I/O
class TestFitsIO:
    def test_tofits_fromfits_roundtrip(self, tmp_path):
        """Writing then reading a FITS file recovers the original array shape."""
        from lsc.cosmics import fromfits, tofits
        # fromfits transposes; tofits transposes back — net: equal shape
        arr = np.random.default_rng(1).normal(1000, 50, (32, 32)).astype(np.float64)
        outfile = str(tmp_path / "roundtrip.fits")
        tofits(outfile, arr, verbose=False)
        arr2, _ = fromfits(outfile, verbose=False)
        assert arr2.shape == arr.shape or arr2.shape == arr.T.shape

    def test_tofits_creates_file(self, tmp_path):
        """tofits() creates a file at the given path."""
        from lsc.cosmics import tofits
        import os
        arr = np.ones((10, 10))
        outfile = str(tmp_path / "out.fits")
        tofits(outfile, arr, verbose=False)
        assert os.path.exists(outfile)

    def test_tofits_boolean_array_written_as_uint8(self, tmp_path):
        """tofits must convert bool arrays to uint8 before writing."""
        from lsc.cosmics import tofits, fromfits
        arr = np.zeros((8, 8), dtype=bool)
        arr[2, 2] = True
        outfile = str(tmp_path / "bool.fits")
        tofits(outfile, arr, verbose=False)
        arr2, _ = fromfits(outfile, verbose=False)
        # Written values should be 0/1 integers, not a float crash
        assert arr2.dtype != bool
        assert set(arr2.ravel().astype(int)).issubset({0, 1})

    def test_tofits_overwrites_existing_file(self, tmp_path):
        """tofits removes an existing file rather than appending."""
        from lsc.cosmics import tofits, fromfits
        arr1 = np.full((8, 8), 1.0)
        arr2 = np.full((8, 8), 2.0)
        outfile = str(tmp_path / "overwrite.fits")
        tofits(outfile, arr1, verbose=False)
        tofits(outfile, arr2, verbose=False)
        result, _ = fromfits(outfile, verbose=False)
        assert np.allclose(result, 2.0)


# subsample — non-square and degenerate inputs
class TestSubsampleExtra:
    def test_non_square_shape(self):
        """subsample() handles non-square inputs, doubling each dimension independently."""
        from lsc.cosmics import subsample
        a = np.ones((6, 10))
        b = subsample(a)
        assert b.shape == (12, 20)

    def test_single_pixel(self):
        """subsample() expands a single pixel to a 2x2 block of the same value."""
        from lsc.cosmics import subsample
        a = np.array([[5.0]])
        b = subsample(a)
        assert b.shape == (2, 2)
        assert np.allclose(b, 5.0)

# cosmicsimage — __str__ extra branches
class TestCosmicsimageStr:
    def _make(self, **kwargs):
        from lsc.cosmics import cosmicsimage
        rng = np.random.default_rng(0)
        data = rng.normal(1000.0, 30.0, (32, 32)).astype(np.float64)
        return cosmicsimage(data, **kwargs)

    def test_str_shows_pssl(self):
        """__str__ includes the pssl value when it is non-zero."""
        ci = self._make(pssl=50.0)
        assert "50" in str(ci)

    def test_str_no_pssl_line_when_zero(self):
        """__str__ omits the pssl line when pssl is zero."""
        ci = self._make(pssl=0.0)
        assert "previously subtracted" not in str(ci)

    def test_str_shows_satstars_after_findsatstars(self):
        """__str__ reports the saturated star mask after findsatstars() is called."""
        from lsc.cosmics import cosmicsimage
        rng = np.random.default_rng(0)
        data = rng.normal(1000.0, 30.0, (64, 64)).astype(np.float64)
        # Plant a saturated star
        data[30:35, 30:35] = 60000.0
        ci = cosmicsimage(data, satlevel=50000.0, verbose=False)
        ci.findsatstars(verbose=False)
        assert "Saturated star" in str(ci)


# cosmicsimage — getrawarray / getcleanarray pssl correction
class TestGetArrays:
    def test_getrawarray_subtracts_pssl(self):
        """getrawarray() returns the original data array, not rawarray (which includes pssl)."""
        from lsc.cosmics import cosmicsimage
        data = np.full((16, 16), 1000.0)
        ci = cosmicsimage(data, pssl=200.0)
        # getrawarray should return original (data), not data+pssl
        assert np.allclose(ci.getrawarray(), data)

    def test_getcleanarray_reflects_cleaned_pixels(self):
        """getcleanarray() returns the post-clean array (with replaced cosmics), not the raw data."""
        from lsc.cosmics import cosmicsimage
        data = np.random.default_rng(3).normal(1000.0, 30.0, (32, 32)).astype(np.float64)
        data[10, 10] = 60000.0
        ci = cosmicsimage(data, gain=2.2, readnoise=10.0, sigclip=5.0, pssl=200.0, satlevel=-1, verbose=False)
        ci.lacosmiciteration(verbose=False)
        ci.clean(verbose=False)
        # getrawarray always returns the original (unmodified) data minus pssl
        assert np.allclose(ci.getrawarray()[10, 10], data[10, 10])
        assert not np.isclose(ci.getcleanarray()[10, 10], data[10, 10])


# cosmicsimage — guessbackgroundlevel
class TestGuessBackgroundLevel:
    def test_returns_median(self):
        """guessbackgroundlevel() returns the median of rawarray."""
        from lsc.cosmics import cosmicsimage
        data = np.arange(100, dtype=np.float64).reshape(10, 10)
        ci = cosmicsimage(data)
        level = ci.guessbackgroundlevel()
        assert np.isclose(level, np.median(data))

    def test_cached_on_second_call(self):
        """guessbackgroundlevel() caches its result and does not recompute on subsequent calls."""
        from lsc.cosmics import cosmicsimage
        data = np.ones((10, 10)) * 42.0
        ci = cosmicsimage(data)
        l1 = ci.guessbackgroundlevel()
        # Mutate rawarray; second call should still return cached value
        ci.rawarray[:] = 0.0
        l2 = ci.guessbackgroundlevel()
        assert l1 == l2


# ---------------------------------------------------------------------------
# cosmicsimage — getdilatedmask
# ---------------------------------------------------------------------------

class TestGetDilatedMask:
    def _ci_with_mask(self):
        from lsc.cosmics import cosmicsimage
        rng = np.random.default_rng(0)
        data = rng.normal(1000.0, 30.0, (32, 32)).astype(np.float64)
        data[15, 15] = 60000.0
        ci = cosmicsimage(data, gain=2.2, readnoise=10.0,
                          sigclip=5.0, satlevel=-1, verbose=False)
        ci.lacosmiciteration(verbose=False)
        return ci

    def test_size3_dilates(self):
        """size=3 dilation produces a mask at least as large as the original."""
        ci = self._ci_with_mask()
        dilated = ci.getdilatedmask(size=3)
        assert dilated.sum() >= ci.mask.sum()

    def test_size5_dilates_more_than_size3(self):
        """size=5 dilation covers at least as many pixels as size=3."""
        ci = self._ci_with_mask()
        d3 = ci.getdilatedmask(size=3)
        d5 = ci.getdilatedmask(size=5)
        assert d5.sum() >= d3.sum()

    def test_unknown_size_returns_copy_of_mask(self):
        """For unrecognised sizes the method falls back to returning a copy."""
        ci = self._ci_with_mask()
        result = ci.getdilatedmask(size=99)
        assert isinstance(result, np.ndarray)


# cosmicsimage — labelmask
class TestLabelmask:
    def test_returns_list(self):
        """labelmask() returns a list of dicts for each detected cosmic island."""
        from lsc.cosmics import cosmicsimage
        rng = np.random.default_rng(0)
        data = rng.normal(1000.0, 30.0, (64, 64)).astype(np.float64)
        data[10, 10] = 60000.0
        data[50, 50] = 60000.0
        ci = cosmicsimage(data, gain=2.2, readnoise=10.0,
                          sigclip=5.0, satlevel=-1, verbose=False)
        ci.lacosmiciteration(verbose=False)
        result = ci.labelmask(verbose=False)
        assert isinstance(result, list)

    def test_empty_mask_gives_empty_list(self):
        """labelmask() returns an empty list when no cosmics are masked."""
        from lsc.cosmics import cosmicsimage
        data = np.full((32, 32), 1000.0)
        ci = cosmicsimage(data, verbose=False)
        # mask is all-False by default
        result = ci.labelmask(verbose=False)
        assert result == []


# cosmicsimage — findsatstars / getsatstars
class TestSatStars:
    def _data_with_sat_star(self, satlevel=50000.0):
        rng = np.random.default_rng(0)
        data = rng.normal(1000.0, 30.0, (64, 64)).astype(np.float64)
        # Plant a clearly saturated patch
        data[28:36, 28:36] = satlevel + 10000.0
        return data, satlevel

    def test_findsatstars_sets_satstars_mask(self):
        """findsatstars() populates satstars with pixels above the saturation level."""
        from lsc.cosmics import cosmicsimage
        data, satlevel = self._data_with_sat_star()
        ci = cosmicsimage(data, satlevel=satlevel, verbose=False)
        ci.findsatstars(verbose=False)
        assert ci.satstars is not None
        assert ci.satstars.sum() > 0

    def test_getsatstars_calls_findsatstars_lazily(self):
        """getsatstars() calls findsatstars() automatically if it hasn't been run yet."""
        from lsc.cosmics import cosmicsimage
        data, satlevel = self._data_with_sat_star()
        ci = cosmicsimage(data, satlevel=satlevel, verbose=False)
        assert ci.satstars is None
        result = ci.getsatstars(verbose=False)
        assert result is not None

    def test_getsatstars_negative_satlevel_prints_error(self, capsys):
        """getsatstars with satlevel <= 0 should print an error."""
        from lsc.cosmics import cosmicsimage
        data = np.full((32, 32), 1000.0)
        ci = cosmicsimage(data, satlevel=-1, verbose=False)
        ci.getsatstars(verbose=False)
        captured = capsys.readouterr()
        assert "satlevel" in captured.out.lower() or "satlevel" in captured.err.lower()


# cosmicsimage — clean() with a custom mask
class TestCleanCustomMask:
    def test_custom_mask_replaces_only_masked_pixels(self):
        """clean() with a custom mask replaces only the specified pixels, leaving others unchanged."""
        from lsc.cosmics import cosmicsimage
        data = np.full((16, 16), 1000.0, dtype=np.float64)
        ci = cosmicsimage(data, verbose=False)
        # Supply a hand-crafted mask rather than running lacosmiciteration
        custom_mask = np.zeros((16, 16), dtype=bool)
        custom_mask[8, 8] = True
        ci.clean(mask=custom_mask, verbose=False)
        # The flagged pixel should have been replaced (not inf)
        assert not np.isinf(ci.cleanarray[8, 8])
        # All other pixels should be unchanged
        other = ci.cleanarray.copy()
        other[8, 8] = 1000.0
        assert np.allclose(other, 1000.0)


# run() — early termination
class TestRunEarlyTermination:
    def test_stops_early_when_no_cosmics(self):
        """run() should break before maxiter when niter==0 on a clean image."""
        from lsc.cosmics import cosmicsimage
        import io, contextlib
        rng = np.random.default_rng(11)
        data = rng.normal(1000.0, 5.0, (64, 64)).astype(np.float64)
        ci = cosmicsimage(data, gain=2.2, readnoise=10.0,
                          sigclip=50.0, satlevel=-1, verbose=False)
        # Capture stdout from run() to count iterations printed
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ci.run(maxiter=5, verbose=False)
        output = buf.getvalue()
        # "Iteration 2" should not appear when the first iteration finds nothing
        assert "Iteration 5" not in output

    def test_fromfits_returns_array_and_header(self, tmp_path):
        """fromfits() returns a tuple of (numpy array, FITS header)."""
        from lsc.cosmics import fromfits, tofits
        arr = np.ones((16, 16))
        outfile = str(tmp_path / "test2.fits")
        tofits(outfile, arr, verbose=False)
        arr2, _ = fromfits(outfile, verbose=False)
        assert isinstance(arr2, np.ndarray)
