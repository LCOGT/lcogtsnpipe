#!/usr/bin/env python
description = ">> make catalogue from table"

import os
from argparse import ArgumentParser
import lsc
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import ascii

def makecatalog(lista, sn2columns=None, ref_coords=None):
    filename_equals = ['filename="{}"'.format(os.path.basename(fn).replace('.sn2.fits', '.fits')) for fn in lista]
    t = Table(lsc.mysqldef.query(['''select filter, filepath, filename, airmass, shortname, dayobs,
                                     zcol1, z1, c1, dz1, dc1, zcol2, z2, c2, dz2, dc2
                                     from photlco, telescopes where photlco.telescopeid=telescopes.id and (''' + 
                                     ' or '.join(filename_equals) + ')'], lsc.conn), masked=True)
    t['filter'] = [lsc.sites.filterst1[filt] for filt in t['filter']]
    if sn2columns is not None:
        print 'Cross-matching {} catalogs. This may take a while...'.format(len(lista))
        catalogs = []
        for image in t:
            cat = Table.read(image['filepath']+image['filename'].replace('.fits', '.sn2.fits'))        
            coords = SkyCoord(cat['ra'], cat['dec'], unit=(u.hourangle, u.deg))
            if ref_coords is None:
                ref_coords = coords
            i0, _, _ = ref_coords.match_to_catalog_sky(coords)
            catalogs.append(cat[i0])
        t['ra'] = [cat['ra'] for cat in catalogs]
        t['dec'] = [cat['dec'] for cat in catalogs]
        t['instmag'] = [cat[sn2columns[0]] for cat in catalogs]
        t['dinstmag'] = [cat[sn2columns[1]] for cat in catalogs]
    for col in t.colnames:
        if np.all(np.isreal(t[col])): # any numerical columns
            t[col].mask = t[col] == 9999.
    return t

def average_in_flux(mag, dmag, axis=None):
    flux = 10**(mag / -2.5)
    dflux = np.log(10) / 2.5 * flux * dmag
    avg_dflux = np.sum(dflux**-2, axis)**-0.5
    avg_flux = np.sum(flux * dflux**-2, axis) * avg_dflux**2
    avg_mag = -2.5 * np.log10(avg_flux)
    avg_dmag = 2.5 / np.log(10) * avg_dflux / avg_flux
    return avg_mag, avg_dmag

if __name__ == "__main__":
    parser = ArgumentParser(description=description)
    parser.add_argument("imglist", help="file containing a list of images on which to run")
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("-e", "--exzp", help='filename for external zero point from different field')
    parser.add_argument("-t", "--typemag", default='fit', choices=['fit', 'ph'], help='type of magnitude')
    parser.add_argument("-F", "--force", action="store_true", help="overwrite existing catalogs")
    args = parser.parse_args()

    with open(args.imglist) as f:
        lista = f.read().splitlines()

    apass_path = lsc.util.getcatalog(lista[0].replace('.sn2.fits', '.fits'), 'apass')
    apass = Table.read(apass_path, format='ascii', names=['ra', 'dec', 'id', 'B', 'Berr',
                       'V', 'Verr', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr'],
                       fill_values=[('9999.000', '0')])
    ref_coords = SkyCoord(apass['ra'], apass['dec'], unit=u.deg)

    if args.typemag == 'fit': # PSF photometry
        targets = makecatalog(lista, ('smagf', 'smagerrf'), ref_coords)
    elif args.typemag == 'ph': # aperture photometry
        targets = makecatalog(lista, ('magp3', 'merrp3'), ref_coords)

    color_to_use = lsc.myloopdef.chosecolor(targets['filter'], True)
    colors_to_calculate = set(sum(color_to_use.values(), []))

    # copy average zero points & color terms from the standards to the science images
    if args.exzp:
        with open(args.exzp) as f:
            lista2 = f.read().splitlines()
        standards = makecatalog(lista2)
        standards = standards.group_by(['dayobs', 'shortname', 'filter', 'zcol1', 'zcol2'])
        for group in standards.groups:
            matches_in_targets = (targets['dayobs'] == group['dayobs'][0]) & (targets['shortname'] == group['shortname'][0]) & (targets['filter'] == group['filter'][0])
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
    extinction = np.array([lsc.sites.extinction[row['shortname'].split()[0].lower()][row['filter']] for row in targets])
    targets['instmag_amcorr'] = targets['instmag'] - np.atleast_2d(extinction * targets['airmass']).T
    targets = targets.group_by(['dayobs', 'shortname'])
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
                dc0 = np.sum(group['dc1'][f0]**-2)**-0.5
                c0 = np.sum(group['c1'][f0] * group['dc1'][f0]**-2) * dc0**2
            else:
                dc0 = 0.
                c0 = np.mean(group['c1'][f0])
            if np.all(group['dc2'][f1]):
                dc1 = np.sum(group['dc2'][f1]**-2)**-0.5
                c1 = np.sum(group['c2'][f1] * group['dc2'][f1]**-2) * dc1**2
            else:
                dc1 = 0.
                c1 = np.mean(group['c2'][f1])
            color = (m0 - m1 + z0 - z1) / (1 - c0 + c1)
            dcolor = np.abs(color) * np.sqrt((dm0**2 + dm1**2 + dz0**2 + dz1**2) / (m0 - m1 + z0 - z1)**2 +
                                             (dc0**2 + dc1**2) / (1 - c0 + c1)**2)
            colors.append(np.tile(color.data, (len(group), 1)))
            dcolors.append(np.tile(dcolor.data, (len(group), 1)))
        targets[filters] = np.vstack(colors)
        targets['d'+filters] = np.vstack(dcolors)

    # calibrate all the catalogs
    zeropoint = np.array([row['z1'] if row['zcol1'] == color_to_use[row['filter']][0] else row['z2'] for row in targets])
    dzeropoint = np.array([row['dz1'] if row['zcol1'] == color_to_use[row['filter']][0] else row['dz2'] for row in targets])
    colorterm = np.array([row['c1'] if row['zcol1'] == color_to_use[row['filter']][0] else row['c2'] for row in targets])
    dcolorterm = np.array([row['dc1'] if row['zcol1'] == color_to_use[row['filter']][0] else row['dc2'] for row in targets])
    color_used = np.array([row[color_to_use[row['filter']][0]] for row in targets]).T
    dcolor_used = np.array([row['d'+color_to_use[row['filter']][0]] for row in targets]).T
    targets['mag'] = (targets['instmag_amcorr'].T + zeropoint + colorterm * color_used).T
    targets['dmag'] = np.sqrt(targets['dinstmag'].T**2 + dzeropoint**2 + dcolorterm**2 * color_used**2 + colorterm**2 + dcolor_used**2).T
    
    # write all the catalogs to files
    for row in targets:
        good = ~row['mag'].mask
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
