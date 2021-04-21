Start with the automatically reduced photometry. I will use NAME as the
name of the object, TARGDATE-TARGDATE as the range of dates for the
science images (images where the target was visible), and TEMPDATE as
the epoch(s) of the reference images. ID refers to the ID number of a
specific file. IN refers to the instrument prefix (kb for SBIG, fl for
Sinistro, fs for Spectral). TEMPTEL is the instrument prefix for the
template image (kb, fl, fs, or SDSS). If you don't give this option, it
will select the same instrument given with -T.

Prepare the Reference Images
============================

LCOGT References
----------------

1.  Look through your reference images and choose the best one for each
    camera--filter combination. Make a note of their ID numbers.

-   ds9 &\
    lscloop.py -n NAME -e TEMPDATE -s checkwcs

    If wcs is off by a little bit, rerun the wcs stage with the
    following command, then run checkwcs:

    lscloop.py -n NAME -e TEMPDATE -s wcs --mode astrometry -F

2.  Mark your chosen images as references.

-   lscloop.py -n NAME -e TEMPDATE -d ID -s template

3.  Run cosmic ray rejection on the reference images.

-   lscloop.py -n NAME -e TEMPDATE --filetype 4 -s cosmic

4.  Generate PSFs for the reference images. I recommend doing this
    interactively for each image (using --show), since these PSFs will
    affect all of your subtractions.

-   lscloop.py -n NAME -e TEMPDATE --filetype 4 -s psf --show

SDSS References
---------------

1.  Choose a set of science images that includes one image with each
    camera--filter combination used. Then run the following command once
    for each of those images, using the ID numbers to choose individual
    frames. This will take a while.

-   lscloop.py -n NAME -e TARGDATE -d ID -s ingestsloan

    Make a note of the TEMPDATE for each SDSS frame you download.

2.  Generate PSFs for the SDSS images. I usually get good results using
    a FWHM of $\sim 5$ pixels for SDSS.

-   ds9 &\
    lscloop.py -n NAME -e TEMPDATE --filetype 4 -s psf --show --fwhm 5
    --use-sextractor

3.  Copy the variance images to the right place for use as cosmic ray
    masks.

-   lscloop.py -n NAME -e TEMPDATE --filetype 4 -s cosmic

Do the Image Subtraction
========================

1.  Run cosmic ray rejection on all the science images. This will take a
    while.

-   lscloop.py -n NAME -e TARGDATE-TARGDATE -s cosmic

2.  In the mean time, check to make sure all the science images have
    good PSFs from the automatic pipeline.

-   lscloop.py -n NAME -e TARGDATE-TARGDATE -s checkpsf

    Generate a new PSF for any images you marked as bad. I usually get
    good results using a FWHM of $\sim 7$ pixels, but this depends on
    the seeing.

    lscloop.py -n NAME -e TARGDATE-TARGDATE -b psf -s psf --show --fwhm
    7

    If the images are saturated, adjust the Dmax value. For my images
    75000 worked, but depends on the gain (saturation is 65535$*$gain):

    lscloop.py -n NAME -e TARGDATE-TARGDATE -b psf -s psf --show --fwhm
    7 --datamax 75000

3.  Once all the cosmic ray rejection is done (for science and reference
    images), run the subtraction. This will take a while. By default,
    --tempdate=19990101-20080101 (useful for SDSS), --temptel=IN,
    --fixpix=False, and --difftype=0 (0 = HOTPANTS, 1 = Optimal).

-   lscloop.py -n NAME -e TARGDATE-TARGDATE --normalize t -T IN
    \[--tempdate TEMPDATE\] \[--temptel TEMPTEL\] \[--fixpix\]
    \[--difftype 1\] -s diff

    If you want, look over the results. Make sure to choose Frame $>$
    Tile on DS9.

    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s checkdiff

Extract the Photometry
======================

By default, --filetype 3 selects both HOTPANTS and Optimally subtracted
images. To select just one or the other, use --difftype 0 for HOTPANTS
or --difftype 1 for Optimal.

1.  First copy the convolved PSFs over.

-   lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s psf

2.  Then calculate the instrumental magnitudes. Use a flat background
    for subtracted images.

-   lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -x 1 -y 1 -s
    psfmag

    If the apmags are all bad, then there was probably a bad pixel near
    the supernova messing it up. To fix this, use --datamax and
    --datamin.

    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s psfmag
    --datamax 100000 --datamin -1000 -F

    If you want, look over the results and redo any bad fits.

    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s checkmag

3.  Calculate the zero point and color term for each set of filters. (If
    the field is not in SDSS, use --field apass for the Sloan filters
    too.)

-   lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -f landolt
    --field apass -s zcat\
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -f sloan -s
    zcat

4.  Calculate the apparent magnitudes. Make sure you don't separate the
    filters if you want to include color corrections.

-   lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s mag

5.  If the lightcurve has outliers or faint data, run the getmag stage
    to interactively check out suspicious points. Can do --type fit to
    check uncalibrated magnitudes.

-   lscloop.py -n NAME -e TARGDATE-TARGDATE -s getmag --type mag --show
    \[-f FILTER\] \[--filetype 3\]
