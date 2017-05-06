#!/usr/bin/env python

description = ">> Calibrate local reference stars"

import numpy as np
import os
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import lsc
from datetime import datetime
from astropy.table import Table, vstack, join
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u

def read_all_catalogs(lista, refcoords=None):
    print 'Cross-matching {} catalogs. This may take a while...'.format(len(lista))
    combined = Table()
    for filename in lista:
        cat = Table.read(filename, format='ascii.commented_header', delimiter='\s', header_start=1)
        if not len(cat):
            continue
        filt = cat.colnames[2]
        cat['filter'] = cat.colnames[2]
        cat['filename'] = os.path.basename(filename)
        cat.rename_column(filt, 'mag')
        cat.rename_column('d'+filt, 'dmag')
        coords = SkyCoord(cat['ra'], cat['dec'], unit=(u.hourangle, u.deg))
        if refcoords is None:
            refcoords = coords
        i0, sep, _ = coords.match_to_catalog_sky(refcoords)
        cat['id'] = i0 + 1 # to match ID in reference catalog (1-indexed)
        combined = vstack([combined, cat[sep.arcsec < 3.5]])
    return combined

if __name__ == "__main__":
    parser = ArgumentParser(description=description)
    parser.add_argument("catlist", help="file containing a list of catalog files on which to run")
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("-f", "--field", default='landolt', choices=['landolt', 'sloan', 'apass'],
                        help='Landolt (UBVRI), Sloan (ugriz), or APASS (BVgri) filters?')
    parser.add_argument("-c", "--catalog", help="use only stars that match this reference catalog")
    parser.add_argument("-o", "--output",  default='{SN}_{field}_{datenow}.cat', help='output filename')
    args = parser.parse_args()
    
    with open(args.catlist) as ff:
        lista = ff.read().splitlines()
    if args.field == 'landolt':
        filterlist = ['U', 'B', 'V', 'R', 'I']
    elif args.field == 'sloan':
        filterlist = ['u', 'g', 'r', 'i', 'z']
    elif args.field == 'apass':
        filterlist = ['B', 'V', 'g', 'r', 'i']
    
    # read reference catalog and find matching stars in our catalogs
    refcat = Table.read(args.catalog, format='ascii', names=['ra', 'dec', 'id', 'B', 'Berr',
                        'V', 'Verr', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr'],
                        fill_values=[('9999.000', '0')])
    refcoords = SkyCoord(refcat['ra'], refcat['dec'], unit=u.deg)
    combined_catalog = read_all_catalogs(lista, refcoords)
    
    # initialize catalog with IRAF astrometry file header
    header = ['BEGIN CATALOG HEADER',
              'nfields 13',
              '    ra    1 0 d degrees %10.6f',
              '    dec   2 0 d degrees %10.6f',
              '    id    3 0 c INDEF %3d']
    for filt in filterlist:
        header.append('    {}    {:2d} 0 r INDEF %6.3f'.format(filt, len(header) - 1))
        header.append('    {}err {:2d} 0 r INDEF %6.3f'.format(filt, len(header) - 1))
    header += ['END CATALOG HEADER', '']
    catalog = Table(names=['ra', 'dec', 'id']+sum([[f, f+'err'] for f in filterlist], []),
                    dtype=[float, float, int, float, float, float, float, float, float,
                    float, float, float, float], meta={'comments': header}, masked=True)

    # combine nights to get a single mag per star per filter
    grouped_by_star = combined_catalog.group_by('id')
    grouped_by_star.remove_columns(['ra', 'dec', 'filename'])
    for star_group in grouped_by_star.groups:
        starid = star_group['id'][0]
        row = {'ra': refcoords[starid-1].ra, 'dec': refcoords[starid-1].dec, 'id': starid}
        for filt in filterlist:
            if filt in star_group['filter']:
                mags = star_group['mag'][star_group['filter'] == filt]
                row[filt] = np.median(mags)
                row[filt+'err'] = np.median(np.abs(mags - row[filt])) * np.sqrt(np.pi / 2)
        catalog.add_row(row)
    
    if args.interactive:
        plt.ion()
        fig = plt.figure(figsize=(11, 8.5))
        for filt in filterlist:
            if filt not in combined_catalog['filter']:
                continue
            nightly_by_filter = join(combined_catalog[combined_catalog['filter'] == filt], catalog[['id', filt]])
            fig.clear()
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            ax1.set_title('Filter: ' + filt)
            ax1.plot(nightly_by_filter['id'], nightly_by_filter['mag'], label='individual images', color='k', alpha=0.5, marker='_', ls='none')
            ax1.errorbar(catalog['id'], catalog[filt], catalog[filt+'err'], label='output (median of images)', marker='o', ls='none')
            if filt in refcat.colnames:
                ax1.plot(refcat['id'], refcat[filt], label='APASS', marker='o', mfc='none', ls='none', zorder=10)
            ax1.invert_yaxis()
            ax1.set_xlabel('Star ID in APASS Catalog')
            ax1.set_ylabel('Apparent Magnitude')
            ax1.legend(loc='best')

            nightly_by_filter['offset'] = nightly_by_filter['mag'] - nightly_by_filter[filt]
            if filt in refcat.colnames:
                nightly_by_filter['offset from APASS'] = nightly_by_filter['mag'] - refcat[filt][nightly_by_filter['id'] - 1]
            nightly_by_filter.remove_columns(['ra', 'dec', 'mag', 'dmag', 'filter', 'id'])
            offsets_by_file = nightly_by_filter.group_by('filename')
            filestats = offsets_by_file.groups.aggregate(np.median)
            filenames, filenums = np.unique(nightly_by_filter['filename'], return_inverse=True)
            ax2.axhline(0., color='k')
            ax2.plot(filenums, nightly_by_filter['offset'], label='individual stars', color='k', alpha=0.5, marker='_', ls='none')
            ax2.plot(filestats['offset'], label='median offset', marker='o', ls='none')
            if filt in refcat.colnames:
                ax2.plot(filestats['offset from APASS'], marker='o', mfc='none', ls='none')
            ax2.set_xticks(range(len(filenames)))
            ax2.set_xticklabels(filenames, rotation='vertical', size='xx-small')
            ax2.set_ylabel('Offset from Median (mag)')
            ax2.legend(loc='best')
            fig.set_tight_layout(True)
            fig.subplots_adjust(bottom=0.28, hspace=0.33)
            fig.canvas.draw_idle()
            
            if filt in refcat.colnames:
                plt.figure(2)
                plt.clf()
                ax3 = plt.subplot(111)
                ax3.axhline(0., color='k')
                diffs = (catalog[filt] - refcat[filt][catalog['id'] - 1]).data
                median_diff = np.median(diffs)
                ax3.plot(catalog[filt], diffs, label='individual stars', marker='o', ls='none')
                ax3.axhline(median_diff, label='median: {:.2f} mag'.format(median_diff), ls='--')
                ax3.set_title('Filter: ' + filt)
                ax3.set_xlabel('Apparent Magnitude in Output Catalog')
                ax3.set_ylabel('Offset from APASS (mag)')
                ax3.legend(loc='best')
                plt.draw()
            
            raw_input('Press enter to continue.')

    # write catalog to file
    snname = os.path.basename(args.catalog).split('_')[0] if args.catalog else 'catalog'
    filename = args.output.format(SN=snname, field=args.field,
                                  datenow=datetime.now().strftime('%Y%m%d_%H_%M'))
    catalog['ra'].format = '%10.6f'
    catalog['dec'].format = '%10.6f'
    catalog['id'].format = '%3d'
    for col in catalog.colnames[3:]:
        catalog[col].format = '%6.3f'
    catalog.write(filename, format='ascii.fixed_width_no_header', delimiter='',
                  fill_values=[(ascii.masked, '9999.0')])
    print 'catalog written to', filename
