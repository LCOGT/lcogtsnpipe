#!/usr/bin/env python

from os.path import splitext
from astropy.io import ascii
import argparse

parser = argparse.ArgumentParser(description='Reformat photometry files to upload to snex2')
parser.add_argument('-n', '--name', type=str, help='filename to reformat')
args = parser.parse_args()

phot = ascii.read(args.name)

lines = ['MJD,filter,mag,error\n']

for row in phot:
    mjd = row[1] - 2400000.5
    filt = row[-2] if isinstance(row[-1], int) else row[-1]
    mag  = [val for val in list(row)[2:-2:2] if val != 9999.][0]
    dmag = [val for val in list(row)[3:-2:2] if val != 0.][0]

    lines.append('{mjd},{filt},{mag},{dmag}\n'.format(mjd=mjd, filt=filt, 
                                                      mag=mag, dmag=dmag))

[filename_root, ext] = splitext(args.name)
filename_snex2 = filename_root + '_snex2' + '.csv'
print('Saving new file as {name}'.format(name=filename_snex2))
with open(filename_snex2, 'w') as f:
    for line in lines:
        f.write(line)
