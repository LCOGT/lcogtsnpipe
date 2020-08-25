#!/usr/bin/env python
description = ">> make catalogue from table"

import os
from argparse import ArgumentParser
import lsc
import numpy.ma as np
from numpy import pi, newaxis
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import matplotlib.pyplot as plt
from datetime import datetime
import sys

def crossmatch(cat0, cat1, threshold=1., racol0='ra', deccol0='dec', racol1='ra', deccol1='dec', right_join=False):
    dra = cat0[racol0] - cat1[racol1][:, newaxis]
    ddec = cat0[deccol0] - cat1[deccol1][:, newaxis]
    sep = np.sqrt(dra**2 + ddec**2) * 3600.
    matches = np.min(sep, axis=1) < threshold
    inds = np.argmin(sep, axis=1)
    out = Table(cat0[inds], masked=True)
    if right_join:
        for col in out.colnames:
            out[col].mask = ~matches
        cat1 = cat1.copy()
    else:
        out = out[matches]
        cat1 = cat1[matches]
    return out, cat1
    
def get_image_data(lista, magcol=None, errcol=None, refcat=None):
    filename_equals = ['filename="{}"'.format(os.path.basename(fn).replace('.sn2.fits', '.fits')) for fn in lista]
    t = Table(lsc.mysqldef.query(['''select filter, filepath, filename, airmass, shortname, dayobs, instrument, targetid,
                                     zcol1, z1, c1, dz1, dc1, zcol2, z2, c2, dz2, dc2, psfmag, psfdmag, apmag, dapmag
                                     from photlco left join telescopes on photlco.telescopeid=telescopes.id where ''' + 
                                     ' or '.join(filename_equals)], lsc.conn), masked=True)
    t['filter'] = [lsc.sites.filterst1[filt] for filt in t['filter']]
    if magcol in t.colnames and errcol in t.colnames:
        t.rename_column(magcol, 'instmag')
        t.rename_column(errcol, 'dinstmag')
    elif magcol is not None and errcol is not None:
        print 'Cross-matching {} catalogs. This may take a while...'.format(len(lista))
        catalogs = []
        badrows = []
        for i, image in enumerate(t):
            sn2file = image['filepath']+image['filename'].replace('.fits', '.sn2.fits')
            if os.path.isfile(sn2file):
                cat = Table.read(sn2file)
                if refcat is None:
                    refcat = cat
                    racol1 = 'ra0'
                    deccol1 = 'dec0'
                else:
                    racol1 = 'ra'
                    deccol1 = 'dec'
                cat_match, _ = crossmatch(cat, refcat, racol0='ra0', deccol0='dec0',
                                          racol1=racol1, deccol1=deccol1, right_join=True)
                catalogs.append(cat_match)
            else:
                badrows.append(i)
        t.remove_rows(badrows)
        t['ra'] = np.array([cat['ra'] for cat in catalogs])
        t['dec'] = np.array([cat['dec'] for cat in catalogs])
        t['instmag'] = np.array([cat[magcol] for cat in catalogs])
        t['dinstmag'] = np.array([cat[errcol] for cat in catalogs])
    for col in t.colnames:
        if t[col].dtype.kind == 'f':
            t[col].mask |= t[col] >= 9999.
    return t

def average_in_flux(mag, dmag, axis=None):
    flux = 10**(mag / -2.5)
    dflux = np.log(10) / 2.5 * flux * dmag
    avg_dflux = np.power(np.sum(np.power(dflux, -2), axis), -0.5)
    avg_flux = np.sum(flux * np.power(dflux, -2), axis) * avg_dflux**2
    avg_mag = -2.5 * np.log10(avg_flux)
    avg_dmag = 2.5 / np.log(10) * np.divide(avg_dflux, avg_flux)
    return avg_mag, avg_dmag

def combine_nights(combined_catalog, filterlist, refcat):
    header = ['BEGIN CATALOG HEADER',
              'nfields 13',
              '    ra    1 0 d degrees %10.6f',
              '    dec   2 0 d degrees %10.6f',
              '    id    3 0 c INDEF %3d']
    for filt in filterlist:
        header.append('    {}    {:2d} 0 r INDEF %6.3f'.format(filt, len(header) - 1))
        header.append('    {}err {:2d} 0 r INDEF %6.3f'.format(filt, len(header) - 1))
    header += ['END CATALOG HEADER', '']
    catalog = Table([refcat['ra'], refcat['dec'], refcat['id']], meta={'comments': header}, masked=True)
    for filt in filterlist:
        mags = combined_catalog['mag'][combined_catalog['filter'] == filt]
        median = np.median(mags, axis=0)
        absdev_mag = mags - median
        mad = np.median(np.abs(absdev_mag), axis=0) * np.sqrt(pi / 2)
        mags.mask |= np.abs(absdev_mag) > 5 * mad
        catalog[filt] = np.median(mags, axis=0)
        catalog[filt+'err'] = np.median(np.abs(mags - catalog[filt]), axis=0) * np.sqrt(pi / 2)
    return catalog

if __name__ == "__main__":
    parser = ArgumentParser(description=description)
    parser.add_argument("imglist", help="file containing a list of images on which to run")
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("-F", "--force", action="store_true", help="overwrite existing catalogs")
    parser.add_argument("-e", "--exzp", help='filename for external zero point from different field')
    parser.add_argument("-s", "--stage", default='abscat', choices=['mag', 'abscat', 'local'],
                        help='calibrate the local sequence (abscat) or the supernova (mag)?')
    parser.add_argument("-t", "--typemag", default='fit', choices=['fit', 'ph'], 
                        help='PSF photometry (fit) or aperture photometry (ph)?')
    parser.add_argument("-f", "--field", choices=['landolt', 'sloan', 'apass'],
                        help='Landolt (UBVRI), Sloan (ugriz), or APASS (BVgri) filters?')
    parser.add_argument("-c", "--catalog", help="use only stars that match this reference catalog")
    parser.add_argument("--minstars", default=0, type=int, help="minimum number of catalog matches for inclusion")
    parser.add_argument("-o", "--output", default='{SN}_{field}_{datenow}.cat', help='output filename')
    args = parser.parse_args()
    
    with open(args.imglist) as f:
        lista = f.read().splitlines()
    if not lista:
        sys.exit('calibratemag.py: ' + args.imglist + ' is empty')
    if args.stage == 'local' and args.exzp is None:
        sys.exit('No file with external zeropoints was passed to calibratemag.py. It may be that no standards were found by myloopdef.get_standards')
    if args.stage in ['abscat', 'local'] and args.catalog is not None:
        try:
            refcat = Table.read(args.catalog, format='ascii', fill_values=[('9999.000', '0')])
            colnames = [row.split()[0] for row in refcat.meta['comments'] if len(row.split()) == 6]
            for old, new in zip(refcat.colnames, colnames):
                refcat.rename_column(old, new)
        except ascii.core.InconsistentTableError: # real Landolt catalogs are different
            refcat = Table.read(args.catalog, format='ascii', names=['id', 'ra', 'dec', 'U', 'B', 'V', 'R', 'I',
                                'vary', 'Uerr', 'Berr', 'Verr', 'Rerr', 'Ierr', 'col13', 'col14', 'col15', 'col16', 'col17'],
                                fill_values=[('99.999', '0'), ('9.9999', '0')])

    if args.typemag == 'fit' and args.stage == 'mag': # PSF photometry for supernova
        targets = get_image_data(lista, 'psfmag', 'psfdmag')
    elif args.typemag == 'ph' and args.stage == 'mag': # aperture photometry for supernova
        targets = get_image_data(lista, 'apmag', 'dapmag')
    elif args.typemag == 'fit': # PSF photometry for local sequence
        targets = get_image_data(lista, 'smagf', 'smagerrf', refcat)
    elif args.typemag == 'ph': # aperture photometry for local sequence
        targets = get_image_data(lista, 'magp3', 'merrp3', refcat)
        
    color_to_use = lsc.sites.chosecolor(targets['filter'], True)
    colors_to_calculate = set(sum(color_to_use.values(), []))

    # copy average zero points & color terms from the standards to the science images
    if args.exzp:
        with open(args.exzp) as f:
            lista2 = f.read().splitlines()
        standards = get_image_data(lista2)
        standards = standards.group_by(['dayobs', 'shortname', 'instrument', 'filter', 'zcol1', 'zcol2'])
        for icol in ['zcol1', 'z1', 'dz1', 'c1', 'dc1', 'zcol2', 'z2', 'dz2', 'c2', 'dc2']:
            targets[icol].mask = True
        for group in standards.groups:
            matches_in_targets = ((targets['dayobs'] == group['dayobs'][0]) & (targets['shortname'] == group['shortname'][0])
                                   & (targets['instrument'] == group['instrument'][0]) & (targets['filter'] == group['filter'][0]))
            if not np.any(matches_in_targets):
                continue
            targets['zcol1'][matches_in_targets] = group['zcol1'][0]
            targets['zcol2'][matches_in_targets] = group['zcol2'][0]
            targets['z1'][matches_in_targets], targets['dz1'][matches_in_targets] = average_in_flux(group['z1'], group['dz1'])
            targets['z2'][matches_in_targets], targets['dz2'][matches_in_targets] = average_in_flux(group['z2'], group['dz2'])
            if np.all(group['dc1']):
                dc1 = np.sum(group['dc1']**-2)**-0.5
                targets['c1'][matches_in_targets] = np.sum(group['c1'] * group['dc1']**-2) * dc1**2
                targets['dc1'][matches_in_targets] = dc1
            else:
                targets['c1'][matches_in_targets] = np.mean(group['c1'])
                targets['dc1'] = 0.
            if np.all(group['dc2']):
                dc2 = np.sum(group['dc2']**-2)**-0.5
                targets['c2'][matches_in_targets] = np.sum(group['c2'] * group['dc2']**-2) * dc2**2
                targets['dc2'][matches_in_targets] = dc2
            else:
                targets['c2'][matches_in_targets] = np.mean(group['c2'])
                targets['dc2'] = 0.

    # generate average colors for each night at each site
    targets['site'] = [row['shortname'].split()[0].lower() if row['shortname'] is not None else None for row in targets]
    extinction = [lsc.sites.extinction[row['site']][row['filter']] for row in targets]
    targets['instmag_amcorr'] = (targets['instmag'].T - extinction * targets['airmass']).T
    targets = targets.group_by(['dayobs', 'shortname', 'instrument'])
    for filters in colors_to_calculate:
        colors, dcolors = [], []
        for group in targets.groups:
            f0 = group['filter'] == filters[0]
            f1 = group['filter'] == filters[1]
            m0, dm0 = average_in_flux(group['instmag_amcorr'][f0], group['dinstmag'][f0], axis=0)
            m1, dm1 = average_in_flux(group['instmag_amcorr'][f1], group['dinstmag'][f1], axis=0)
            z0, dz0 = average_in_flux(group['z1'][f0], group['dz1'][f0])
            z1, dz1 = average_in_flux(group['z2'][f1], group['dz2'][f1])
            if np.all(group['dc1'][f0]):
                dc0 = np.sum(np.power(group['dc1'][f0], -2))**-0.5
                c0 = np.sum(group['c1'][f0] * np.power(group['dc1'][f0], -2)) * dc0**2
            else:
                dc0 = 0.
                c0 = np.mean(group['c1'][f0])
            if np.all(group['dc2'][f1]):
                dc1 = np.sum(np.power(group['dc2'][f1], -2))**-0.5
                c1 = np.sum(group['c2'][f1] * np.power(group['dc2'][f1], -2)) * dc1**2
            else:
                dc1 = 0.
                c1 = np.mean(group['c2'][f1])
            color = np.divide(m0 - m1 + z0 - z1, 1 - c0 + c1)
            dcolor = np.abs(color) * np.sqrt(
                        np.divide(dm0**2 + dm1**2 + dz0**2 + dz1**2, (m0 - m1 + z0 - z1)**2)
                        + np.divide(dc0**2 + dc1**2, (1 - c0 + c1)**2)
                                            )
            for row in group:
                colors.append(color)
                dcolors.append(dcolor)
        targets[filters] = np.array(colors)
        targets['d'+filters] = np.array(dcolors)

    # calibrate all the instrumental magnitudes
    zcol = [color_to_use[row['filter']][0] if color_to_use[row['filter']] else row['filter']*2 for row in targets]
    zeropoint = np.choose(zcol == targets['zcol1'], [targets['z2'], targets['z1']])
    dzeropoint = np.choose(zcol == targets['zcol1'], [targets['dz2'], targets['dz1']])
    colorterm = np.choose(zcol == targets['zcol1'], [targets['c2'], targets['c1']])
    dcolorterm = np.choose(zcol == targets['zcol1'], [targets['dc2'], targets['dc1']])
    uzcol, izcol = np.unique(zcol, return_inverse=True)
    color_used = np.choose(izcol, [targets[col].T if col in targets.colnames else 0. for col in uzcol]).filled(0.) # if no other filter, skip color correction
    dcolor_used = np.choose(izcol, [targets['d'+col].T if col in targets.colnames else 0. for col in uzcol]).filled(0.)
    targets['mag'] = (targets['instmag_amcorr'].T + zeropoint + colorterm * color_used).T
    targets['dmag'] = np.sqrt(targets['dinstmag'].T**2 + dzeropoint**2 + dcolorterm**2 * color_used**2 + colorterm**2 * dcolor_used**2).T

    if args.stage == 'mag':
        # write mag & dmag to database
        targets['dmag'].mask = targets['mag'].mask
        query = 'INSERT INTO photlco (filename, targetid, mag, dmag) VALUES\n'
        query += ',\n'.join(['("{}", {}, {}, {})'.format(row['filename'], row['targetid'], row['mag'], row['dmag']) for row in targets.filled(9999.)])
        query += '\nON DUPLICATE KEY UPDATE mag=VALUES(mag), dmag=VALUES(dmag)'
        print query
        lsc.mysqldef.query([query], lsc.conn)
    elif args.stage == 'abscat': 
       # write all the catalogs to files & put filename in database
        for row in targets:
            good = ~row['mag'].mask
            if not np.any(good):
                print 'no good magnitudes for', row['filename']
                lsc.mysqldef.updatevalue('photlco', 'abscat', 'X', row['filename'])
                continue
            outtab = Table([row['ra'][good].T, row['dec'][good].T, row['mag'][good].T, row['dmag'][good].T],
                           meta={'comments': ['daophot+standardfield', '         ra           dec       {0}     d{0}'.format(row['filter'])]})
            outtab['col2'].format = '%6.3f'
            outtab['col3'].format = '%5.3f'
            outfile = row['filename'].replace('.fits', '.cat')
            try:
                outtab.write(row['filepath'] + outfile, format='ascii.fixed_width_no_header', delimiter='',
                             overwrite=args.force, fill_values=[(ascii.masked, '9999.0')])
                lsc.mysqldef.updatevalue('photlco', 'abscat', outfile, row['filename'])
            except IOError as e:
                print e, '-- use -F to overwrite'
    elif args.stage == 'local':
        if args.field == 'landolt':
            filterlist = ['U', 'B', 'V', 'R', 'I']
        elif args.field == 'sloan':
            filterlist = ['u', 'g', 'r', 'i', 'z']
        elif args.field == 'apass':
            filterlist = ['B', 'V', 'g', 'r', 'i']
        else:
            raise Exception('Need to give --field (landolt, sloan, apass) for -s local')
        # make master catalog and write to file
        catalog = combine_nights(targets, filterlist, refcat)
        catfile = os.path.basename(args.catalog)
        if args.interactive:
            plt.ion()
            fig = plt.figure(1, figsize=(11, 8.5))
            for filt in filterlist:
                nightly_by_filter = targets[(targets['filter'] == filt) & (np.sum(~targets['mag'].mask, axis=1) > args.minstars)]
                if not nightly_by_filter:
                    print 'no calibrated stars in', filt
                    continue
                fig.clear()
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
                ax1.set_title('Filter: ' + filt)
                lgd_handle = ax1.plot(nightly_by_filter['mag'].data.T, color='k', alpha=0.5, marker='_', ls='none')[:1]
                lgd_handle.append(ax1.errorbar(range(len(catalog)), catalog[filt], yerr=catalog[filt+'err'], marker='o', ls='none'))
                lgd_label = ['individual images', 'output (median of images)']
                if filt in refcat.colnames:
                    lgd_handle += ax1.plot(refcat[filt], marker='o', mfc='none', ls='none', zorder=10)
                    lgd_label.append(catfile)
                ax1.invert_yaxis()
                ax1.set_xlabel('Star ID in {} (line number)'.format(catfile))
                ax1.set_ylabel('Apparent Magnitude')
                ax1.legend(lgd_handle, lgd_label, loc='best')

                nightly_by_filter['offset'] = nightly_by_filter['mag'] - catalog[filt].T
                if filt in refcat.colnames:
                    nightly_by_filter['offset from refcat'] = nightly_by_filter['mag'] - refcat[filt].T
                ax2.axhline(0., color='k', lw=1)
                lgd_handle = ax2.plot(nightly_by_filter['offset'], color='k', alpha=0.5, marker='_', ls='none')[:1]
                lgd_handle += ax2.plot(np.median(nightly_by_filter['offset'], axis=1), marker='o', ls='none')
                lgd_label = ['individual stars', 'median offset']
                if filt in refcat.colnames:
                    lgd_handle += ax2.plot(np.median(nightly_by_filter['offset from refcat'], axis=1), marker='o', mfc='none', ls='none')
                    lgd_label.append('offset from ' + catfile)
                ax2.set_xticks(range(len(nightly_by_filter)))
                ax2.set_xticklabels(nightly_by_filter['filename'], rotation='vertical', size='xx-small')
                ax2.set_ylabel('Offset from Median (mag)')
                ax2.legend(lgd_handle, lgd_label, loc='best')
                fig.set_tight_layout(True)
                fig.subplots_adjust(bottom=0.28, hspace=0.33)
                fig.canvas.draw_idle()
                
                if filt in refcat.colnames:
                    plt.figure(2)
                    plt.clf()
                    ax3 = plt.subplot(111)
                    ax3.axhline(0., color='k', lw=1)
                    diffs = (catalog[filt] - refcat[filt]).data
                    median_diff = np.median(diffs)
                    ax3.plot(catalog[filt], diffs, label='individual stars', marker='o', ls='none')
                    ax3.axhline(median_diff, label='median: {:.3f} mag'.format(median_diff), ls='--', lw=1)
                    ax3.set_title('Filter: ' + filt)
                    ax3.set_xlabel('Apparent Magnitude in Output Catalog')
                    ax3.set_ylabel('Output - {} (mag)'.format(catfile))
                    ax3.legend(loc='best')
                    plt.draw()
                
                raw_input('Press enter to continue.')

        snname = os.path.basename(catfile).split('_')[0] if args.catalog else 'catalog'
        filename = args.output.format(SN=snname, field=args.field,
                                      datenow=datetime.now().strftime('%Y%m%d_%H_%M'))
        catalog['ra'].format = '%10.6f'
        catalog['dec'].format = '%10.6f'
        for col in catalog.colnames[3:]:
            catalog[col].format = '%6.3f'
        catalog.write(filename, format='ascii.fixed_width_no_header', delimiter='',
                      fill_values=[(ascii.masked, '9999.0')])
        print 'catalog written to', filename
