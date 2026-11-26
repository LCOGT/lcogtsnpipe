import re 

BASE_HEADER_MAP = {
    'object': 'OBJECT',
    'date-obs': 'DATE-OBS',
    'ut': 'DATE-OBS',
    'date-night': 'DAY-OBS',
    'RA': 'RA',
    'DEC': 'DEC',
    'CAT-RA': 'CAT-RA',
    'CAT-DEC': 'CAT-DEC',
    'datamin': -100.0,
    'datamax': 'SATURATE',
    'observer': 'OBSERVER',
    'exptime': 'EXPTIME',
    'wcserr': 'WCSERR',
    'instrume': 'INSTRUME',
    'JD': 'MJD-OBS',
    'mjd': 'MJD-OBS',
    'filter': 'FILTER',
    'gain': 'GAIN',
    'ron': 'RDNOISE',
    'airmass': 'AIRMASS',
    'type': 'OBSTYPE',
    'propid': 'PROPID',
    'userid': 'USERID',
    'telescop': 'TELESCOP',
}

NATURAL_COLOR_TERMS = {
    'UUB': 0.0, 'uug': 0.0, 'BUB': 0.0, 'BBV': 0.0, 'VBV': 0.0, 'VVR': 0.0,
    'gug': 0.0, 'ggr': 0.0, 'RVR': 0.0, 'RRI': 0.0, 'rrz': 0.0, 'zrz': 0.0,
    'rgr': 0.0, 'rri': 0.0, 'iri': 0.0, 'iiz': 0.0, 'IRI': 0.0, 'ziz': 0.0
}

SLOANPRIME_STANDARD_COLOR_TERMS = {
    'UUB': 0.059, 'uug': 0.0, 'BUB': -0.095, 'BBV': 0.06, 'VBV': 0.03, 'VVR': -0.059,
    'gug': 0.13, 'ggr': 0.054, 'RVR': -0.028, 'RRI': -0.033, 'rrz': 0.0, 'zrz': 0.0,
    'ggi': 0.0, 'igi': 0.0, 'rgr': 0.003, 'rri': -0.007, 'iri': 0.028, 'iiz': 0.110,
    'IRI': 0.013, 'ziz': -0.16
}

ZERO_COLOR_TERMS = {
    'uug': 0.0, 'gug': 0.0, 'ggr': 0.0, 'rgr': 0.0, 'rri': 0.0, 'iri': 0.0,
    'iiz': 0.0, 'ziz': 0.0, 'UUB': 0.0, 'BUB': 0.0, 'BBV': 0.0, 'VBV': 0.0,
    'VVR': 0.0, 'RVR': 0.0, 'RRI': 0.0, 'IRI': 0.0
}

# BVgri color terms from Valenti et al. 2016, MNRAS, 459, 3939
SINISTRO_COLOR_TERMS = {
    'uug': 0.0, 'ggr': 0.109, 'rri': 0.027, 'iri': 0.036, 'BBV': -0.024, 'VBV': -0.014,
    'UUB': 0.059, 'BUB': -0.095, 'VVR': -0.059, 'RVR': -0.028, 'RRI': -0.033, 'IRI': 0.013,
    'ziz': -0.04
}

MUSCAT_COLOR_TERMS = {
    'uug': 0.0, 'ggr': 0.0087, 'rri': 0.0166, 'iri': 0.0217, 'BBV': 0.0, 'VBV': 0.0,
    'UUB': 0.0, 'BUB': 0.0, 'VVR': 0.0, 'RVR': 0.0, 'RRI': 0.0, 'IRI': 0.0,
    'ziz': 0.0152
}

SPECTRAL_COLOR_TERMS = {
    'UUB': 0.0, 'uug': 0.0, 'BUB': 0.0, 'BBV': 0.0, 'VBV': 0.0, 'VVR': 0.0,
    'gug': 0.0, 'ggr': 0.105, 'RVR': 0.0, 'RRI': 0.0, 'rrz': 0.0, 'zrz': 0.0,
    'rgr': 0.013, 'rri': 0.029, 'iri': 0.0874, 'iiz': 0.0, 'IRI': 0.0, 'ziz': -0.15
}

SPECTRAL_SLOANPRIME_COLOR_TERMS = {
    'UUB': 0.0, 'uug': 0.0, 'BUB': 0.0, 'BBV': 0.0, 'VBV': 0.0, 'VVR': 0.0,
    'gug': 0.0, 'ggr': 0.0, 'RVR': 0.0, 'RRI': 0.0, 'rrz': 0.0, 'zrz': 0.0,
    'rgr': 0.0, 'rri': 0.0, 'iri': 0.027, 'iiz': 0.0, 'IRI': 0.0, 'ziz': 0.0
}

LEGACY_SPECTRAL_PIXEL_SCALES = {'fs01': 0.304, 'fs02': 0.30104, 'fs03': 0.304}


class InstrumentProfile:
    """Shared instrument behavior. Subclasses must define the required fields."""

    # Required: every camera must define these
    name: str
    db_type: str | None
    pixel_scale: float | None
    template_name: str | None
    template_pixel_scale: float | None
    template_imagesize: int | None
    color_terms: dict | None

    # Shared defaults: override only when a camera differs
    astrometry_search_size = 20
    saturation_level = 45000.0
    sextractor_saturation_cap = 55000.0
    banzai_datamax = 75000.0
    sloanprime_color_terms = SLOANPRIME_STANDARD_COLOR_TERMS

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        missing = [field for field in InstrumentProfile.__annotations__ if not hasattr(cls, field)]
        if missing:
            raise TypeError(f"{cls.__name__} must define: {', '.join(missing)}")

    def get_header_keys(self, hdr):
        """Standard BANZAI / modern LCO header keyword mapping."""
        return BASE_HEADER_MAP.copy()

    def get_color_terms(self, calib=None):
        """Return calibration color terms for this instrument."""
        if calib == 'natural':
            return NATURAL_COLOR_TERMS.copy()
        if calib == 'sloanprime':
            return self.sloanprime_color_terms.copy()
        if self.color_terms is not None:
            return self.color_terms.copy()
        print('No color terms exist for telescope/instrument set up. No color term applied, '
              'use --unfix to calculate a color term from the field stars in the image')
        return ZERO_COLOR_TERMS.copy()


class SinistroProfile(InstrumentProfile):
    name = 'sinistro'
    db_type = 'Sinistro'
    pixel_scale = 0.389
    template_name = 'sinistro'
    template_pixel_scale = 0.387
    template_imagesize = 4020
    color_terms = SINISTRO_COLOR_TERMS
    astrometry_search_size = 40


class SophiaProfile(SinistroProfile):
    """New camera in the Sinistro instrument package, with the same pixels and optics."""
    name = 'sophia'
    db_type = 'SOPHIA'
    color_terms = None  # not measured yet; use --unfix
    saturation_level = 236000.0  # ~260 ke- full well at 1.1 e-/ADU
    sextractor_saturation_cap = saturation_level
    banzai_datamax = saturation_level


class MuscatProfile(InstrumentProfile):
    name = 'muscat'
    db_type = 'MUSCAT'
    pixel_scale = 0.27
    template_name = 'muscat'
    template_pixel_scale = 0.27
    template_imagesize = 2020
    color_terms = MUSCAT_COLOR_TERMS


class SbigProfile(InstrumentProfile):
    name = 'sbig'
    db_type = 'SBIG'
    pixel_scale = 0.464
    template_name = 'sbig'
    template_pixel_scale = 0.467
    template_imagesize = 2030
    color_terms = None


class QhyProfile(InstrumentProfile):
    name = 'qhy'
    db_type = 'QHY'
    pixel_scale = 0.734
    template_name = 'qhy'
    template_pixel_scale = 0.74
    template_imagesize = 162
    color_terms = None


class SpectralProfile(InstrumentProfile):
    name = 'spectral'
    db_type = 'Spectral'
    pixel_scale = 0.30
    template_name = 'spectral'
    template_pixel_scale = 0.30104
    template_imagesize = 2020
    color_terms = SPECTRAL_COLOR_TERMS
    sloanprime_color_terms = SPECTRAL_SLOANPRIME_COLOR_TERMS

    def get_header_keys(self, hdr):
        date_obs = str(hdr.get('DATE-OBS') or '')
        if not date_obs or date_obs >= '2014-04-01':
            return super().get_header_keys(hdr)

        keys = {
            **BASE_HEADER_MAP,
            'datamax': 60000.0,
            'wcserr': 'WCS_ERR',
            'JD': 'MJD',
            'mjd': 'MJD',
            'ron': next((key for key in ('RDNOISE', 'READNOIS') if key in hdr), 'ron'),
            'telescop': 'TELID',
        }
        pixel_scale = LEGACY_SPECTRAL_PIXEL_SCALES.get(str(hdr.get('INSTRUME') or '').lower())
        if pixel_scale:
            keys['pixscale'] = pixel_scale
        return keys


class EmccdProfile(SpectralProfile):
    name = 'emccd'
    db_type = None
    pixel_scale = 0.278


class GenericProfile(InstrumentProfile):
    name = 'generic'
    db_type = None
    pixel_scale = None
    template_name = None
    template_pixel_scale = None
    template_imagesize = None
    color_terms = None
    saturation_level = 32000.0

    def get_header_keys(self, hdr):
        return {
            'object': 'OBJECT',
            'RA': 'RA',
            'DEC': 'DEC',
            'CAT-RA': 'RA',
            'CAT-DEC': 'DEC',
            'ron': 'RDNOISE',
            'mjd': 'MJD-OBS',
            'ut': 'DATE-OBS',
            'date-obs': 'DATE-OBS',
            'date-night': 'DAY-OBS',
            'datamax': 'SATURATE',
        }


# ---------------------------------------------------------------------------
# Instrument identification
# ---------------------------------------------------------------------------

PROFILES = {
    profile.name: profile
    for profile in (
        SinistroProfile(),
        SophiaProfile(),
        MuscatProfile(),
        SbigProfile(),
        SpectralProfile(),
        EmccdProfile(),
        QhyProfile(),
        GenericProfile(),
    )
}

CAMERA_PREFIXES = {
    'fa': 'sinistro',
    'fl': 'sinistro',
    'kb': 'sbig',
    'sq': 'qhy',
    'fs': 'spectral',
    'em': 'emccd',
}

FIRST_SOPHIA_CAMERA = 60


def get_instrument(hdr_or_name):
    """Return the InstrumentProfile for a FITS header or instrument name."""
    name = hdr_or_name.get('INSTRUME') if hasattr(hdr_or_name, 'get') else hdr_or_name
    name = str(name or '').strip().lower()
    if name in PROFILES:
        return PROFILES[name]

    camera = re.fullmatch(r'([a-z]{2})(\d*)', name)
    if not camera:
        return PROFILES['generic']
    prefix, number = camera.groups()
    if prefix == 'ep':  # MuSCAT and SOPHIA share this prefix
        is_sophia = number and int(number) >= FIRST_SOPHIA_CAMERA
        return PROFILES['sophia' if is_sophia else 'muscat']
    return PROFILES[CAMERA_PREFIXES.get(prefix, 'generic')]
