import pytest

from lsc.instruments import (
    get_instrument,
    InstrumentProfile,
    SinistroProfile,
    SophiaProfile,
    MuscatProfile,
    SbigProfile,
    SpectralProfile,
    EmccdProfile,
    QhyProfile,
    GenericProfile,
)


def test_profile_must_define_every_field():
    with pytest.raises(TypeError, match='db_type'):
        class IncompleteProfile(InstrumentProfile):
            name = 'incomplete'


def test_resolve_by_name():
    assert isinstance(get_instrument('fa01'), SinistroProfile)
    assert isinstance(get_instrument('fl05'), SinistroProfile)
    assert isinstance(get_instrument('sinistro'), SinistroProfile)
    assert isinstance(get_instrument('kb12'), SbigProfile)
    assert isinstance(get_instrument('sbig'), SbigProfile)
    assert isinstance(get_instrument('fs02'), SpectralProfile)
    assert isinstance(get_instrument('spectral'), SpectralProfile)
    assert isinstance(get_instrument('em01'), EmccdProfile)
    assert isinstance(get_instrument('emccd'), EmccdProfile)
    assert isinstance(get_instrument('sq01'), QhyProfile)
    assert isinstance(get_instrument('qhy'), QhyProfile)
    assert isinstance(get_instrument('muscat'), MuscatProfile)
    assert isinstance(get_instrument('sophia'), SophiaProfile)
    assert isinstance(get_instrument('unknown'), GenericProfile)
    assert isinstance(get_instrument(''), GenericProfile)


def test_resolve_ep_sophia_vs_muscat():
    # ep60 and above is SOPHIA (1m)
    assert isinstance(get_instrument('ep60'), SophiaProfile)
    assert isinstance(get_instrument('ep61'), SophiaProfile)
    assert isinstance(get_instrument('ep70'), SophiaProfile)

    # ep below 60 is MuSCAT (2m)
    assert isinstance(get_instrument('ep01'), MuscatProfile)
    assert isinstance(get_instrument('ep05'), MuscatProfile)
    assert isinstance(get_instrument('ep59'), MuscatProfile)

    # Bare 'ep' defaults to MuSCAT for backward compatibility
    assert isinstance(get_instrument('ep'), MuscatProfile)


def test_external_instruments_are_generic():
    for name in ('FLI', 'KeplerCam', 'DECam', '1m0', 'all'):
        assert isinstance(get_instrument(name), GenericProfile)


def test_resolve_from_header():
    hdr_sinistro = {'INSTRUME': 'fa01', 'TELESCOP': '1m0-01'}
    assert isinstance(get_instrument(hdr_sinistro), SinistroProfile)

    hdr_sophia = {'INSTRUME': 'ep60', 'TELESCOP': '1m0-01'}
    assert isinstance(get_instrument(hdr_sophia), SophiaProfile)

    hdr_muscat = {'INSTRUME': 'ep05', 'TELESCOP': '2m0-01'}
    assert isinstance(get_instrument(hdr_muscat), MuscatProfile)

    hdr_spectral = {'INSTRUME': 'fs02', 'TELESCOP': '2m0-02'}
    assert isinstance(get_instrument(hdr_spectral), SpectralProfile)


def test_db_type():
    assert get_instrument('fa01').db_type == 'Sinistro'
    assert get_instrument('fl05').db_type == 'Sinistro'
    assert get_instrument('ep60').db_type == 'SOPHIA'
    assert get_instrument('ep05').db_type == 'MUSCAT'
    assert get_instrument('kb12').db_type == 'SBIG'
    assert get_instrument('fs02').db_type == 'Spectral'
    assert get_instrument('sq01').db_type == 'QHY'
    assert get_instrument('em01').db_type is None
    assert get_instrument('generic').db_type is None


def test_pixel_scale_and_astrometry():
    sinistro = get_instrument('fa01')
    assert sinistro.pixel_scale == 0.389
    assert sinistro.astrometry_search_size == 40

    sophia = get_instrument('ep60')
    assert sophia.pixel_scale == 0.389
    assert sophia.astrometry_search_size == 40

    assert get_instrument('unknown').pixel_scale is None

    muscat = get_instrument('ep05')
    assert muscat.pixel_scale == 0.27
    assert muscat.astrometry_search_size == 20

    sbig = get_instrument('kb12')
    assert sbig.pixel_scale == 0.464
    assert sbig.astrometry_search_size == 20

    spectral = get_instrument('fs02')
    assert spectral.pixel_scale == 0.30
    assert spectral.astrometry_search_size == 20

    emccd = get_instrument('em01')
    assert emccd.pixel_scale == 0.278
    assert emccd.astrometry_search_size == 20

    qhy = get_instrument('sq01')
    assert qhy.pixel_scale == 0.734
    assert qhy.astrometry_search_size == 20


def test_saturation_limits():
    sophia = get_instrument('ep60')
    assert sophia.saturation_level == 236000.0
    assert sophia.sextractor_saturation_cap == 236000.0
    assert sophia.banzai_datamax == 236000.0

    sinistro = get_instrument('fa01')
    assert sinistro.saturation_level == 45000.0
    assert sinistro.sextractor_saturation_cap == 55000.0
    assert sinistro.banzai_datamax == 75000.0


def test_header_keys_spectral_legacy():
    hdr_old = {'INSTRUME': 'fs02', 'DATE-OBS': '2013-05-01T12:00:00', 'RDNOISE': 8.5}
    inst = get_instrument(hdr_old)
    keys = inst.get_header_keys(hdr_old)
    assert keys['datamax'] == 60000.0
    assert keys['wcserr'] == 'WCS_ERR'
    assert keys['JD'] == 'MJD'
    assert keys['mjd'] == 'MJD'
    assert keys['telescop'] == 'TELID'
    assert keys['ron'] == 'RDNOISE'
    assert keys['pixscale'] == 0.30104

    hdr_new = {'INSTRUME': 'fs02', 'DATE-OBS': '2015-05-01T12:00:00'}
    keys_new = inst.get_header_keys(hdr_new)
    assert keys_new['datamax'] == 'SATURATE'
    assert keys_new['wcserr'] == 'WCSERR'
    assert keys_new['JD'] == 'MJD-OBS'
    assert keys_new['telescop'] == 'TELESCOP'


def test_color_terms():
    sinistro_terms = get_instrument('fa01').get_color_terms()
    assert sinistro_terms['ggr'] == 0.109

    # SOPHIA color terms are not measured yet
    sophia_terms = get_instrument('ep60').get_color_terms()
    assert all(val == 0.0 for val in sophia_terms.values())

    # MuSCAT has its own color terms
    muscat_terms = get_instrument('ep05').get_color_terms()
    assert muscat_terms['ggr'] == 0.0087

    # Natural calibration gives 0 for all
    nat_terms = get_instrument('fa01').get_color_terms('natural')
    assert all(val == 0.0 for val in nat_terms.values())

    # Sloanprime terms
    sinistro_sp = get_instrument('fa01').get_color_terms('sloanprime')
    assert sinistro_sp['UUB'] == 0.059

    spectral_sp = get_instrument('fs02').get_color_terms('sloanprime')
    assert spectral_sp['iri'] == 0.027


def test_readkey3_integration():
    from lsc.util import readkey3

    hdr_sinistro = {
        'INSTRUME': 'fa01',
        'TELESCOP': '1m0-01',
        'DATE-OBS': '2023-01-15T04:30:00',
        'DAY-OBS': '20230115',
        'OBJECT': 'SN2023test',
        'MJD-OBS': 59959.1875,
        'EXPTIME': 120.0,
        'FILTER': 'rp',
        'GAIN': 3.1,
        'RDNOISE': 10.0,
        'RA': '12:34:56.7',
        'DEC': '-20:30:40.0',
        'WCSERR': 0,
        'SATURATE': 65000,
    }
    assert readkey3(hdr_sinistro, 'instrume') == 'fa01'
    assert readkey3(hdr_sinistro, 'telescop') == '1m0-01'
    assert readkey3(hdr_sinistro, 'date-obs') == '20230115'
    assert readkey3(hdr_sinistro, 'ut') == '04:30:00'
    assert readkey3(hdr_sinistro, 'object') == 'SN2023test'
    assert readkey3(hdr_sinistro, 'filter') == 'rp'
    assert readkey3(hdr_sinistro, 'mjd') == 59959.1875
    assert readkey3(hdr_sinistro, 'datamin') == -100.0

    hdr_sophia = {
        'INSTRUME': 'ep60',
        'TELESCOP': '1m0-01',
        'DATE-OBS': '2026-10-15T02:00:00',
        'DAY-OBS': '20261015',
        'OBJECT': 'SN2026test',
        'MJD-OBS': 60000.0,
        'FILTER': 'zs',
        'GAIN': 1.1,
        'RDNOISE': 8.0,
        'SATURATE': 260000,
    }
    assert readkey3(hdr_sophia, 'instrume') == 'ep60'
    assert readkey3(hdr_sophia, 'filter') == 'zs'
    assert readkey3(hdr_sophia, 'datamin') == -100.0

    hdr_old_spectral = {
        'INSTRUME': 'fs02',
        'DATE-OBS': '2013-05-01T12:00:00',
        'DAY-OBS': '20130501',
        'MJD': 56413.5,
        'TELID': '2m0-01',
        'WCS_ERR': 0,
        'RDNOISE': 8.0,
    }
    assert readkey3(hdr_old_spectral, 'mjd') == 56413.5
    assert readkey3(hdr_old_spectral, 'wcserr') == 0
    assert readkey3(hdr_old_spectral, 'datamax') == 60000.0
    assert readkey3(hdr_old_spectral, 'telescop') == '2m0-01'


def test_guess_instrument_type_function():
    from lsc.mysqldef import guess_instrument_type

    assert guess_instrument_type('fa01') == 'Sinistro'
    assert guess_instrument_type('fl05') == 'Sinistro'
    assert guess_instrument_type('ep60') == 'SOPHIA'
    assert guess_instrument_type('ep05') == 'MUSCAT'
    assert guess_instrument_type('kb12') == 'SBIG'
    assert guess_instrument_type('fs02') == 'Spectral'
    assert guess_instrument_type('sq01') == 'QHY'
    assert guess_instrument_type('em01') is None
    assert guess_instrument_type('unknown') is None

