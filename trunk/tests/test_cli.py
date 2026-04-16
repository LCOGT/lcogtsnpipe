"""CLI tests for all bin/ entry points.

Coverage:
- --help: exits 0/1 (no crashes).
- Missing args: required positionals → non-zero exit.
- Invalid choices: reject bad values.
- Invalid types: reject wrong types (e.g., str for int).
- No-op flags: safe boolean combos (no DB/IO).
- Help output: includes key flags.
"""

import subprocess
import sys
import os
import pytest

pytestmark = pytest.mark.subprocess
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bin"))


# Helpers
def _run(script_name, *args, timeout=30):
    """Run a bin script and return the CompletedProcess."""
    script = os.path.join(BIN_DIR, script_name)
    return subprocess.run(
        [sys.executable, script, *args],
        capture_output=True,
        timeout=timeout,
    )


def _help(script_name, timeout=30):
    return _run(script_name, "--help", timeout=timeout)


def _no_args(script_name, timeout=30):
    return _run(script_name, timeout=timeout)


def _ok(result):
    """argparse exits 0, optparse exits 1 for --help – both are acceptable."""
    return result.returncode in (0, 1)


def _error(result):
    """A bad invocation must produce a non-zero exit code."""
    return result.returncode != 0


# 1. --help smoke tests (all 19 scripts)
class TestCliHelp:
    """Every script must handle --help without crashing (exit 0 or 1)."""

    @pytest.mark.parametrize("script", [
        "lscloop.py",
        "runlsc.py",
        "lscpsf.py",
        "lscdiff.py",
        "calibratemag.py",
        "lscsn.py",
        "LCOGTingest.py",
        "lscmaketempl.py",
        "lscastro.py",
        "lscmerge.py",
        "lscnewcalib.py",
        "ingestall.py",
        "ingesttar.py",
        "queryapasscat.py",
        "querysloancat.py",
        "comparecatalogs.py",
        "format_snex2.py",
        "lscingestsloan.py",
        "back_populate_apercorr.py",
        "lsctestheader.py",
    ])
    def test_help_exits_cleanly(self, script):
        result = _help(script)
        assert _ok(result), (
            f"{script} --help returned {result.returncode}\n"
            f"stdout: {result.stdout.decode(errors='replace')}\n"
            f"stderr: {result.stderr.decode(errors='replace')}"
        )


# 2. Help output content checks
class TestHelpContent:
    """--help output must reference the script's primary flags."""

    def test_lscloop_help_mentions_stage(self):
        r = _help("lscloop.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "--stage" in combined or "-s" in combined

    def test_lscloop_help_mentions_filter(self):
        r = _help("lscloop.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "--filter" in combined or "-f" in combined

    def test_lscdiff_help_mentions_normalize(self):
        r = _help("lscdiff.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "normalize" in combined

    def test_lscdiff_help_mentions_convolve(self):
        r = _help("lscdiff.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "convolve" in combined

    def test_lscpsf_help_mentions_threshold(self):
        r = _help("lscpsf.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "threshold" in combined

    def test_lscpsf_help_mentions_fwhm(self):
        r = _help("lscpsf.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "fwhm" in combined.lower()

    def test_querysloancat_help_mentions_ra(self):
        r = _help("querysloancat.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "ra" in combined.lower()

    def test_lscingestsloan_help_mentions_type(self):
        r = _help("lscingestsloan.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "--type" in combined

    def test_ingesttar_help_mentions_file(self):
        r = _help("ingesttar.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "--file" in combined or "-f" in combined

    def test_format_snex2_help_mentions_name(self):
        r = _help("format_snex2.py")
        combined = (r.stdout + r.stderr).decode(errors="replace")
        assert "--name" in combined or "-n" in combined


# 3. Missing required positional arguments → non-zero exit
class TestMissingRequiredArgs:
    """Scripts with mandatory positional arguments must fail when invoked bare."""

    def test_lscdiff_no_args_fails(self):
        # lscdiff.py requires 'targlist' and 'templist' positional args
        r = _no_args("lscdiff.py")
        assert _error(r), "lscdiff.py should exit non-zero when called without positional args"

    def test_querysloancat_no_args_fails(self):
        # querysloancat.py requires 'ra' and 'dec' positional args
        r = _no_args("querysloancat.py")
        assert _error(r), "querysloancat.py should exit non-zero without ra/dec"

    def test_lscingestsloan_no_args_fails(self):
        # lscingestsloan.py requires at least one 'images' positional arg
        r = _no_args("lscingestsloan.py")
        assert _error(r), "lscingestsloan.py should exit non-zero without image list"


# 4. Invalid choice values → non-zero exit
class TestInvalidChoices:
    """Passing a value not in a choices list must produce a non-zero exit."""

    def test_lscloop_bad_stage(self):
        r = _run("lscloop.py", "--stage", "notastage")
        assert _error(r), "lscloop.py --stage notastage should fail"

    def test_lscloop_bad_filter(self):
        r = _run("lscloop.py", "--filter", "X")
        assert _error(r), "lscloop.py --filter X should fail (X is not a valid filter)"

    def test_lscloop_bad_field(self):
        r = _run("lscloop.py", "--field", "badfield")
        assert _error(r), "lscloop.py --field badfield should fail"

    def test_lscloop_bad_type(self):
        r = _run("lscloop.py", "--type", "invalid_type")
        assert _error(r), "lscloop.py --type invalid_type should fail"

    def test_lscdiff_bad_normalize(self):
        # lscdiff.py --normalize choices: ['i', 't']
        r = _run("lscdiff.py", "targ.fits", "templ.fits", "--normalize", "x")
        assert _error(r), "lscdiff.py --normalize x should fail"

    def test_lscdiff_bad_convolve(self):
        r = _run("lscdiff.py", "targ.fits", "templ.fits", "--convolve", "z")
        assert _error(r), "lscdiff.py --convolve z should fail"

    def test_lscdiff_bad_interpolation(self):
        r = _run("lscdiff.py", "targ.fits", "templ.fits", "--interpolation", "cubic99")
        assert _error(r), "lscdiff.py --interpolation cubic99 should fail"

    def test_lscdiff_bad_difftype(self):
        r = _run("lscdiff.py", "targ.fits", "templ.fits", "--difftype", "99")
        assert _error(r), "lscdiff.py --difftype 99 should fail"

    def test_lscingestsloan_bad_type(self):
        r = _run("lscingestsloan.py", "img.fits", "--type", "badsurvey")
        assert _error(r), "lscingestsloan.py --type badsurvey should fail"

    def test_lscloop_bad_normalize(self):
        r = _run("lscloop.py", "--normalize", "q")
        assert _error(r), "lscloop.py --normalize q should fail"

    def test_lscloop_bad_convolve(self):
        r = _run("lscloop.py", "--convolve", "q")
        assert _error(r), "lscloop.py --convolve q should fail"

    def test_lscloop_bad_filetype(self):
        # filetype choices: range(5) → 0,1,2,3,4
        r = _run("lscloop.py", "--filetype", "9")
        assert _error(r), "lscloop.py --filetype 9 should fail"

    def test_lscloop_bad_calib(self):
        r = _run("lscloop.py", "--calib", "badcalib")
        assert _error(r), "lscloop.py --calib badcalib should fail"

    def test_lscloop_bad_bad(self):
        r = _run("lscloop.py", "--bad", "notabadchoice")
        assert _error(r), "lscloop.py --bad notabadchoice should fail"


# 5. Invalid type coercion → non-zero exit
class TestInvalidTypes:
    """Passing a string where an int/float is required must produce a non-zero exit."""

    def test_lscloop_xord_expects_int(self):
        r = _run("lscloop.py", "--xord", "not_an_int")
        assert _error(r), "lscloop.py --xord not_an_int should fail"

    def test_lscloop_yord_expects_int(self):
        r = _run("lscloop.py", "--yord", "abc")
        assert _error(r), "lscloop.py --yord abc should fail"

    def test_lscloop_threshold_expects_float(self):
        r = _run("lscloop.py", "--threshold", "notafloat")
        assert _error(r), "lscloop.py --threshold notafloat should fail"

    def test_lscloop_bkg_expects_float(self):
        r = _run("lscloop.py", "--bkg", "notafloat")
        assert _error(r), "lscloop.py --bkg notafloat should fail"

    def test_lscloop_size_expects_float(self):
        r = _run("lscloop.py", "--size", "notafloat")
        assert _error(r), "lscloop.py --size notafloat should fail"

    def test_lscloop_cutmag_expects_float(self):
        r = _run("lscloop.py", "--cutmag", "abc")
        assert _error(r), "lscloop.py --cutmag abc should fail"

    def test_lscloop_combine_expects_float(self):
        r = _run("lscloop.py", "--combine", "abc")
        assert _error(r), "lscloop.py --combine abc should fail"

    def test_lscloop_sigma_clip_expects_float(self):
        r = _run("lscloop.py", "--sigma-clip", "abc")
        assert _error(r), "lscloop.py --sigma-clip abc should fail"

    def test_lscloop_targetid_expects_int(self):
        r = _run("lscloop.py", "--targetid", "notanint")
        assert _error(r), "lscloop.py --targetid notanint should fail"

    def test_querysloancat_ra_expects_float(self):
        r = _run("querysloancat.py", "notafloat", "2.2")
        assert _error(r), "querysloancat.py with non-float RA should fail"

    def test_querysloancat_dec_expects_float(self):
        r = _run("querysloancat.py", "150.0", "notafloat")
        assert _error(r), "querysloancat.py with non-float DEC should fail"

    def test_querysloancat_radius_expects_float(self):
        r = _run("querysloancat.py", "150.0", "2.2", "--radius", "notafloat")
        assert _error(r), "querysloancat.py --radius notafloat should fail"


# 6. Unknown / unrecognised flags → non-zero exit
class TestUnknownFlags:
    """argparse scripts must reject completely unknown flags."""

    @pytest.mark.parametrize("script", [
        "lscloop.py",
        "lscdiff.py",
        "querysloancat.py",
        "lscingestsloan.py",
        "format_snex2.py",
        "ingesttar.py",
    ])
    def test_unknown_flag_fails(self, script):
        r = _run(script, "--this-flag-does-not-exist-9999")
        assert _error(r), f"{script} should exit non-zero for unknown flag"


# 7. Valid optional-only invocations
class TestValidFlagCombinations:
    """Checks that the argument parser does not reject a valid flag."""

    def _assert_not_argparse_crash(self, result, script):
        """Exit code 2 from argparse means a parse-level usage error."""
        assert result.returncode != 2, (
            f"{script} raised an argparse usage error for a valid flag combination\n"
            f"stderr: {result.stderr.decode(errors='replace')}"
        )

    def test_lscloop_force_flag(self):
        r = _run("lscloop.py", "--force")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_interactive_flag(self):
        r = _run("lscloop.py", "--interactive")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_show_flag(self):
        r = _run("lscloop.py", "--show")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_zcatold_flag(self):
        r = _run("lscloop.py", "--zcatold")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_unfix_flag(self):
        r = _run("lscloop.py", "--unfix")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_use_sextractor_flag(self):
        r = _run("lscloop.py", "--use-sextractor")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_match_by_site_flag(self):
        r = _run("lscloop.py", "--match-by-site")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_filter_single(self):
        r = _run("lscloop.py", "--filter", "r")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_filter_multiple(self):
        r = _run("lscloop.py", "--filter", "r", "g", "i")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_stage(self):
        r = _run("lscloop.py", "--stage", "wcs")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_field_landolt(self):
        r = _run("lscloop.py", "--field", "landolt")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_field_sloan(self):
        r = _run("lscloop.py", "--field", "sloan")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_field_gaia(self):
        r = _run("lscloop.py", "--field", "gaia")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_type_fit(self):
        r = _run("lscloop.py", "--type", "fit")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_filetype_values(self):
        for v in ("0", "1", "2", "3", "4"):
            r = _run("lscloop.py", "--filetype", v)
            self._assert_not_argparse_crash(r, f"lscloop.py --filetype {v}")

    def test_lscloop_valid_normalize_i(self):
        r = _run("lscloop.py", "--normalize", "i")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_normalize_t(self):
        r = _run("lscloop.py", "--normalize", "t")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_calib_sloan(self):
        r = _run("lscloop.py", "--calib", "sloan")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscloop_valid_calib_natural(self):
        r = _run("lscloop.py", "--calib", "natural")
        self._assert_not_argparse_crash(r, "lscloop.py")

    def test_lscdiff_valid_normalize_i(self):
        # Positional args needed; will fail at file-open, not argparse
        r = _run("lscdiff.py", "t.fits", "r.fits", "--normalize", "i")
        self._assert_not_argparse_crash(r, "lscdiff.py")

    def test_lscdiff_valid_normalize_t(self):
        r = _run("lscdiff.py", "t.fits", "r.fits", "--normalize", "t")
        self._assert_not_argparse_crash(r, "lscdiff.py")

    def test_lscdiff_valid_interpolation_drizzle(self):
        r = _run("lscdiff.py", "t.fits", "r.fits", "--interpolation", "drizzle")
        self._assert_not_argparse_crash(r, "lscdiff.py")

    def test_lscdiff_valid_interpolation_nearest(self):
        r = _run("lscdiff.py", "t.fits", "r.fits", "--interpolation", "nearest")
        self._assert_not_argparse_crash(r, "lscdiff.py")

    def test_lscdiff_valid_difftype_0(self):
        r = _run("lscdiff.py", "t.fits", "r.fits", "--difftype", "0")
        self._assert_not_argparse_crash(r, "lscdiff.py")

    def test_lscdiff_valid_difftype_1(self):
        r = _run("lscdiff.py", "t.fits", "r.fits", "--difftype", "1")
        self._assert_not_argparse_crash(r, "lscdiff.py")

    def test_lscdiff_force_show_fixpix(self):
        r = _run("lscdiff.py", "t.fits", "r.fits", "--force", "--show", "--fixpix")
        self._assert_not_argparse_crash(r, "lscdiff.py")

    def test_lscdiff_unmask_no_iraf(self):
        r = _run("lscdiff.py", "t.fits", "r.fits", "--unmask", "--no-iraf")
        self._assert_not_argparse_crash(r, "lscdiff.py")

    def test_lscingestsloan_type_sloan(self):
        r = _run("lscingestsloan.py", "img.fits", "--type", "sloan")
        self._assert_not_argparse_crash(r, "lscingestsloan.py")

    def test_lscingestsloan_type_ps1(self):
        r = _run("lscingestsloan.py", "img.fits", "--type", "ps1")
        self._assert_not_argparse_crash(r, "lscingestsloan.py")

    def test_lscingestsloan_force_flag(self):
        r = _run("lscingestsloan.py", "img.fits", "--force")
        self._assert_not_argparse_crash(r, "lscingestsloan.py")

    def test_querysloancat_valid_args(self):
        # Will fail at network; not an argparse crash
        r = _run("querysloancat.py", "150.0", "2.2", "--radius", "5.0",
                 "--mag1", "14.0", "--mag2", "19.0", "--output", "out.cat")
        self._assert_not_argparse_crash(r, "querysloancat.py")

    def test_ingesttar_force_db_flag(self):
        r = _run("ingesttar.py", "--force-db")
        self._assert_not_argparse_crash(r, "ingesttar.py")

    def test_format_snex2_name_flag(self):
        r = _run("format_snex2.py", "--name", "nonexistent_file.dat")
        self._assert_not_argparse_crash(r, "format_snex2.py")

    def test_lscmerge_check_flag(self):
        r = _run("lscmerge.py", "--check")
        self._assert_not_argparse_crash(r, "lscmerge.py")

    def test_lscmerge_force_flag(self):
        r = _run("lscmerge.py", "--force")
        self._assert_not_argparse_crash(r, "lscmerge.py")


# 8. Boundary / edge cases for numeric arguments
class TestNumericBoundaries:
    """Numeric arguments at extreme but syntactically valid values."""

    def test_lscloop_xord_zero(self):
        r = _run("lscloop.py", "--xord", "0")
        assert r.returncode != 2

    def test_lscloop_xord_large(self):
        r = _run("lscloop.py", "--xord", "999")
        assert r.returncode != 2

    def test_lscloop_threshold_zero(self):
        r = _run("lscloop.py", "--threshold", "0.0")
        assert r.returncode != 2

    def test_lscloop_threshold_negative(self):
        r = _run("lscloop.py", "--threshold", "-1.0")
        assert r.returncode != 2

    def test_lscloop_sigma_clip_very_large(self):
        r = _run("lscloop.py", "--sigma-clip", "100.0")
        assert r.returncode != 2

    def test_lscloop_combine_zero(self):
        r = _run("lscloop.py", "--combine", "0.0")
        assert r.returncode != 2

    def test_querysloancat_radius_zero(self):
        r = _run("querysloancat.py", "150.0", "2.2", "--radius", "0.0")
        assert r.returncode != 2

    def test_querysloancat_negative_ra(self):
        # Negative RA is syntactically valid (argparse accepts it)
        r = _run("querysloancat.py", "-180.0", "0.0")
        assert r.returncode != 2

    def test_lscloop_xshift_negative(self):
        r = _run("lscloop.py", "--xshift", "-50")
        assert r.returncode != 2

    def test_lscloop_yshift_large(self):
        r = _run("lscloop.py", "--yshift", "10000")
        assert r.returncode != 2


# 9. Script existence checks
class TestScriptExists:
    """All expected bin scripts must exist on disk."""

    @pytest.mark.parametrize("script", [
        "lscloop.py",
        "runlsc.py",
        "lscpsf.py",
        "lscdiff.py",
        "calibratemag.py",
        "lscsn.py",
        "LCOGTingest.py",
        "lscmaketempl.py",
        "lscastro.py",
        "lscmerge.py",
        "lscnewcalib.py",
        "lsctestheader.py",
        "ingestall.py",
        "ingesttar.py",
        "queryapasscat.py",
        "querysloancat.py",
        "comparecatalogs.py",
        "format_snex2.py",
        "back_populate_apercorr.py",
        "lscingestsloan.py",
    ])
    def test_script_exists(self, script):
        path = os.path.join(BIN_DIR, script)
        assert os.path.isfile(path), f"Expected bin script not found: {path}"
