"""
Tests for lsc.sites — pure data structures and chosecolor() helper.
"""
import pytest

pytestmark = pytest.mark.unit


def test_extinction_keys():
    """Verify all expected observatory sites are present in the extinction table."""
    from lsc.sites import extinction
    expected_sites = {"lsc", "coj", "ogg", "elp", "cpt", "tfn", None, "PS1", "SDSS"}
    assert expected_sites.issubset(set(extinction.keys()))


def test_extinction_values_are_floats():
    """Ensure every per-band extinction coefficient is a plain float."""
    from lsc.sites import extinction
    for site, bands in extinction.items():
        for band, val in bands.items():
            assert isinstance(val, float), f"{site}/{band} = {val!r} is not float"


def test_extinction_positive():
    """Extinction coefficients must be non-negative for all sites and bands."""
    from lsc.sites import extinction
    for site, bands in extinction.items():
        for band, val in bands.items():
            assert val >= 0.0, f"Negative extinction at {site}/{band}: {val}"


def test_filterst_canonical_bands():
    """Check that all standard photometric bands are registered in filterst."""
    from lsc.sites import filterst
    for band in "UBVRIugrizw":
        assert band in filterst, f"Band '{band}' missing from filterst"


def test_filterst1_reverse_of_filterst():
    """filterst1 must be the inverse mapping of filterst."""
    from lsc.sites import filterst, filterst1
    for canonical, aliases in filterst.items():
        if isinstance(canonical, str) and len(canonical) == 1:
            for alias in aliases:
                assert filterst1[alias] == canonical, (
                    f"filterst1['{alias}'] = '{filterst1[alias]}', expected '{canonical}'"
                )


def test_chosecolor_basic():
    """Verify adjacent color pairs (BV, VR) are assigned correctly for Landolt filters."""
    from lsc.sites import chosecolor
    color = chosecolor(["B", "V", "R"])
    # BV should be available for both B and V
    assert "BV" in color["B"]
    assert "BV" in color["V"]
    # VR should be available for both V and R
    assert "VR" in color["V"]
    assert "VR" in color["R"]
    # RI not possible — I not in allfilter
    assert "RI" not in color.get("R", [])


def test_chosecolor_sloan():
    """Verify adjacent color pairs (gr, ri, iz) are assigned correctly for Sloan filters."""
    from lsc.sites import chosecolor
    color = chosecolor(["u", "g", "r", "i", "z"])
    assert "gr" in color["g"]
    assert "ri" in color["r"]
    assert "iz" in color["i"]


def test_chosecolor_usegood_selects_preferred():
    """With usegood=True, each filter should resolve to a single preferred color pair."""
    from lsc.sites import chosecolor
    color = chosecolor(["B", "V", "R", "I"], usegood=True)
    # With usegood, V prefers VR, B prefers BV
    assert color["V"] == ["VR"]
    assert color["B"] == ["BV"]


def test_chosecolor_missing_pair():
    """Only one filter — no color pair is possible."""
    from lsc.sites import chosecolor
    color = chosecolor(["g"])
    assert color["g"] == []


def test_filterst_sloan_group():
    """All aliases for Sloan bands (ugriz) must appear in the 'sloan' group."""
    from lsc.sites import filterst
    for band in "ugriz":
        for alias in filterst[band]:
            assert alias in filterst["sloan"]


def test_filterst_landolt_group():
    """All aliases for Landolt bands (UBVRI) must appear in the 'landolt' group."""
    from lsc.sites import filterst
    for band in "UBVRI":
        for alias in filterst[band]:
            assert alias in filterst["landolt"]

def test_ps1_sdss_aliases_match_source_sites():
    """PS1 and SDSS extinction entries must be the same objects as ogg and elp respectively."""
    from lsc.sites import extinction
    assert extinction["PS1"] is extinction["ogg"]
    assert extinction["SDSS"] is extinction["elp"]


def test_filterst1_no_orphan_keys():
    """Every key in filterst1 must trace back to a filterst alias."""
    from lsc.sites import filterst, filterst1
    all_aliases = {alias for k, v in filterst.items()
                   if isinstance(k, str) and len(k) == 1
                   for alias in v}
    for alias in filterst1:
        assert alias in all_aliases, f"filterst1 has orphan key '{alias}'"


def test_chosecolor_usegood_fallback_when_preferred_unavailable():
    """With usegood=True and only one filter available, result should be empty not an error."""
    from lsc.sites import chosecolor
    # V alone — VR can't form, should be empty not an error
    color = chosecolor(["V"], usegood=True)
    assert color["V"] == []


def test_filterst_apass_group():
    """All aliases for APASS bands (BVgri) must appear in the 'apass' group."""
    from lsc.sites import filterst
    for band in "BVgri":
        for alias in filterst[band]:
            assert alias in filterst["apass"]


def test_filterst_gaia_contains_all():
    """The 'gaia' group must contain aliases for every supported photometric band."""
    from lsc.sites import filterst
    for band in "UBVRIugrizw":
        for alias in filterst[band]:
            assert alias in filterst["gaia"]


def test_filterst_empty_string_group():
    """The empty-string group must equal the union of the landolt and sloan groups."""
    from lsc.sites import filterst
    combined = set(filterst["landolt"]) | set(filterst["sloan"])
    assert set(filterst[""]) == combined


def test_filterst_band_lists_nonempty():
    """Every canonical band must have at least one alias registered in filterst."""
    from lsc.sites import filterst
    for band in "UBVRIugrizw":
        assert len(filterst[band]) > 0, f"filterst['{band}'] is empty"


def test_chosecolor_multiple_pairs_usegood_collapses():
    """With usegood=True and multiple possible pairs, each filter collapses to one preferred pair."""
    from lsc.sites import chosecolor
    color = chosecolor(["U", "B", "V"], usegood=True)
    assert color["B"] == ["BV"]   # prefers BV over UB
    assert color["U"] == ["UB"]
    assert color["V"] == ["BV"]


def test_all_sites_have_same_band_keys():
    """Every site in the extinction table must define exactly the same set of photometric bands."""
    from lsc.sites import extinction
    reference_bands = set(extinction[None].keys())
    for site, bands in extinction.items():
        assert set(bands.keys()) == reference_bands, (
            f"Site '{site}' has different bands: {set(bands.keys()) ^ reference_bands}"
        )

def test_filterst1_no_collisions():
    """Ensure no alias points to two different canonical bands."""
    from lsc.sites import filterst
    all_aliases = []
    for canonical, aliases in filterst.items():
        if isinstance(canonical, str) and len(canonical) == 1:
            all_aliases.extend(aliases)
    assert len(all_aliases) == len(set(all_aliases)), "Duplicate alias found across different bands"

def test_extinction_default_not_empty():
    """Ensure the 'None' (default) site has the required base data."""
    from lsc.sites import extinction
    assert len(extinction[None]) > 0
    # Example: Check that the default V-band extinction is roughly standard
    assert 0.0 <= extinction[None]['V'] <= 1.0 

def test_chosecolor_invalid_input():
    """Ensure the helper doesn't crash with nonsense input."""
    from lsc.sites import chosecolor
    # Should return empty dict or handle gracefully, not raise KeyError
    assert chosecolor(["NOT_A_BAND"]) == {"NOT_A_BAND": []}