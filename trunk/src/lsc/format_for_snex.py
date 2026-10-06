import getpass
import json
import os
import sys

import requests
from astropy.time import Time

import lsc

DEFAULT_GROUPIDCODE = 1703768065789
_conn = None
_auth = None
_warned = False


def _db():
    global _conn
    if _conn is None:
        _conn = lsc.mysqldef.dbConnect(*lsc.mysqldef.getconnection('lcogt2'))
    return _conn


def _query(sql):
    return lsc.mysqldef.query([sql], _db())


def configured():
    return bool(os.environ.get('SNEX_URL', '') and os.environ.get('SNEX_API_TOKEN', ''))


def _credentials():
    global _auth
    if _auth is None:
        _auth = (lsc.util.userinput('SNEx username: '), getpass.getpass('SNEx password: '))
    return _auth


class SNExAuthError(RuntimeError):
    pass


def _api(method, path, **kwargs):
    url = os.environ['SNEX_URL'].rstrip('/')
    headers = {'Accept': 'application/json'}
    if os.environ.get('SNEX_API_TOKEN', ''):
        headers['Authorization'] = 'Token ' + os.environ['SNEX_API_TOKEN']
        response = requests.request(method, url + path, headers=headers, allow_redirects=False, **kwargs)
    else:
        response = requests.request(method, url + path, auth=_credentials(), headers=headers, allow_redirects=False, **kwargs)
    if response.status_code in (301, 302, 401) or (response.status_code == 403 and (
            'Invalid token' in response.text or 'credentials were not provided' in response.text)):
        raise SNExAuthError('SNEx at {} answered {} {}; check SNEX_URL and SNEX_API_TOKEN'.format(
            url, response.status_code, response.headers.get('Location', '')))
    if not response.ok:
        raise RuntimeError('{} {} returned {}: {}'.format(method, path, response.status_code, response.text[:300]))
    return response.json() if response.content else None


def groups_for(groupidcode):
    groupidcode = int(groupidcode or DEFAULT_GROUPIDCODE)
    return [{'name': g['name']} for g in _query('select name, idcode from `groups`')
            if g['idcode'] and int(g['idcode']) & groupidcode]


def _with_groups(payload, groups, encode=False):
    if groups is None:
        return payload
    return dict(payload, groups=json.dumps(groups) if encode else groups)


def _post_with_groups(groups, post):
    default_groups = groups_for(None) or None
    try:
        return post(groups or default_groups)
    except RuntimeError as e:
        if 'does not exist' not in str(e):
            raise
        print('{}; using the default groups'.format(e))
        return post(default_groups)


def photometry_groups(row):
    groupidcode = row['groupidcode']
    if not groupidcode:
        raw = _query('select groupidcode from photlcoraw where filename like "{}%"'.format(row['filename'].split('.')[0]))
        groupidcode = raw[0]['groupidcode'] if raw else None
    return groups_for(groupidcode)


def is_standard(targetid):
    return bool(_query('select targets.id from targets join classifications on classifications.id=targets.classificationid '
                       'where targets.id={:d} and classifications.name="Standard"'.format(int(targetid))))


def get_target(targetid, cache=None):
    if cache is not None and targetid in cache:
        return cache[targetid]
    row = _query('select * from targets where id={:d}'.format(int(targetid)))[0]
    standard = is_standard(targetid)
    names = {}
    for n in _query('select name from targetnames where targetid={:d} order by id'.format(int(targetid))):
        names.setdefault(n['name'].replace(' ', '').lower(), n['name'])
    names = list(names.values())
    params = {'name': names[0], 'type': 'SIDEREAL', 'permissions': 'PRIVATE',
              'ra': float(row['ra0']), 'dec': float(row['dec0']),
              'standard': standard, 'aliases': [{'name': n} for n in names[1:]],
              'comment': 'Target created by the pipeline'}
    if standard:
        target = _api('POST', '/api/targets/', json=params)
    else:
        target = _post_with_groups(groups_for(row['groupidcode']), lambda groups: _api(
            'POST', '/api/targets/', json=_with_groups(params, groups)))
    target_id = target['id']
    standard = target.get('standard', standard)
    if cache is not None:
        cache[targetid] = (target_id, standard)
    return target_id, standard


def photometry_value(row):
    value = {
        'magnitude': row['mag'],
        'error': row['dmag'],
        'filter': row['filter'],
        'telescope': row['telescope'],
        'instrument': row['instrument'],
        'basename': row['filename'].split('.')[0],
        'exptime': row['exptime'],
        'fwhm': row['fwhm'],
        'wcs': None if row['wcs'] is None else int(row['wcs']),
        'psfx': row['psfx'],
        'psfy': row['psfy'],
        'psfmag': row['psfmag'],
        'psfdmag': row['psfdmag'],
        'apmag': row['apmag'],
        'dapmag': row['dapmag'],
        'background_subtracted': row['filetype'] == 3,
    }
    if value['background_subtracted']:
        value['reduction_type'] = 'manual'
        value['template_source'] = next((source for source in ('SDSS', 'PS1') if source in row['filename']), 'LCO')
        algorithm = {0: 'Hotpants', 1: 'PyZOGY'}.get(row['difftype'])
        if algorithm:
            value['subtraction_algorithm'] = algorithm
        value['difference_image'] = row['filepath'] + row['filename']
        template = _query('select photlco.filepath, photlco.filename from photpairing join photlco '
                          'on photlco.filename=photpairing.nametemplate where photpairing.nameout="{}"'.format(row['filename']))
        if template:
            value['template_image'] = template[0]['filepath'] + template[0]['filename']
    return value


def upload_photometry(filenames, final=True):
    if not os.environ.get('SNEX_URL', ''):
        raise RuntimeError('SNEX_URL must be set to upload to SNEx')
    if not configured() and not sys.stdin.isatty():
        raise RuntimeError('SNEX_API_TOKEN must be set to upload to SNEx without a terminal')
    targets, failed = {}, []
    results = {'created': 0, 'updated': 0, 'unchanged': 0, 'kept_final': 0}
    for filename in filenames:
        try:
            rows = _query('select * from photlco where filename="{}"'.format(filename))
            if not rows or rows[0]['filetype'] not in (1, 3):
                continue
            row = rows[0]
            target_id, standard = get_target(row['targetid'], targets)
            value = photometry_value(row)
            if final:
                value['final_reduction'] = True
            datum = {
                'target': target_id,
                'timestamp': Time(row['mjd'], format='mjd').isot + 'Z',
                'source_name': 'LCO',
                'source_location': value['basename'],
                'value': value,
            }
            if standard:
                posted = _api('POST', '/api/photometry/', json=dict(datum, groups=[]))
            else:
                posted = _post_with_groups(photometry_groups(row), lambda groups: _api(
                    'POST', '/api/photometry/', json=_with_groups(datum, groups)))
            results[posted['result']] += 1
        except requests.ConnectionError as e:
            raise RuntimeError('Could not connect to SNEx at {}: {}'.format(os.environ['SNEX_URL'], e))
        except SNExAuthError:
            raise
        except Exception as e:
            failed.append(filename)
            print('SNEx upload failed for {}: {}'.format(filename, e))
    print('SNEx photometry: {created} created, {updated} updated, {unchanged} unchanged, {kept_final} kept as final, '
          '{failed} failed'.format(failed=len(failed), **results))
    return failed


def _post_image(target_id, path, data_product_type, groups):
    if not os.path.isfile(path):
        path += '.fz'
    product_id = os.path.basename(path).split('.fits')[0]

    def post(groups):
        with open(path, 'rb') as f:
            return _api('POST', '/api/dataproducts/', files={'file': (os.path.basename(path), f, 'image/fits')},
                        data=_with_groups({'target': target_id, 'data_product_type': data_product_type,
                                           'product_id': product_id}, groups, encode=True))
    return _post_with_groups(groups, post)


def upload_difference_images(filenames):
    if not configured():
        raise RuntimeError('SNEX_URL and SNEX_API_TOKEN must be set to upload images to SNEx')
    targets, posted, failed = {}, 0, []
    for filename in filenames:
        try:
            rows = _query('select * from photlco where filename="{}"'.format(filename))
            if not rows or rows[0]['filetype'] != 3:
                continue
            row = rows[0]
            target_id, standard = get_target(row['targetid'], targets)
            value = photometry_value(row)
            for key, data_product_type in (('difference_image', 'difference_image'), ('template_image', 'template_image')):
                if value.get(key) and not _post_image(target_id, value[key], data_product_type, photometry_groups(row)).get('already_posted'):
                    posted += 1
        except (requests.ConnectionError, SNExAuthError):
            raise
        except Exception as e:
            failed.append(filename)
            print('SNEx image upload failed for {}: {}'.format(filename, e))
    print('SNEx difference images: {} files posted, {} frames failed'.format(posted, len(failed)))
    return failed


def upload_raw_spectrum(filepath, filename):
    global _warned
    posted_marker = os.path.join(filepath, filename.replace('.fits', '.posted'))
    if os.path.exists(posted_marker):
        return
    if not configured():
        if not _warned:
            print('ERROR: SNEX_URL and SNEX_API_TOKEN must be set to post raw spectra to the SNEx FLOYDS inbox')
            _warned = True
        return
    try:
        rows = _query('select * from speclcoraw where filename="{}"'.format(filename))
        if not rows:
            return
        row = rows[0]
        if row['type'] in ('LAMPFLAT', 'ARC') or 'bad' in (row['filepath'] or ''):
            return
        names = [n['name'] for n in _query('select name from targetnames where targetid={:d}'.format(int(row['targetid'])))]
        if any('test_' in n for n in names):
            return
        if row['propid'] == 'FLOYDS standards' or is_standard(row['targetid']):
            return
        target_id, standard = get_target(row['targetid'])
        if standard:
            return
        params = {
            'target': target_id,
            'data_product_type': 'raw_spectrum',
            'product_id': filename.split('.')[0],
            'extra_data': json.dumps({'path': row['filepath'], 'filename': filename,
                                      'dateobs': str(row['dateobs']), 'propid': row['propid']}),
        }
        def post(groups):
            with open(os.path.join(filepath, filename.replace('.fits', '.png')), 'rb') as f:
                return _api('POST', '/api/dataproducts/', data=_with_groups(params, groups, encode=True), files={'thumbnail': f})
        if not _post_with_groups(groups_for(row['groupidcode']), post).get('already_posted'):
            print('posted {} to the SNEx FLOYDS inbox'.format(filename))
        open(posted_marker, 'w').close()
    except Exception as e:
        print('SNEx FLOYDS inbox post failed for {}: {}'.format(filename, e))
