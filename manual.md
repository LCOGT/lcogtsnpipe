# Table of Contents
* [Ingest Data](#Ingest-data)
* [Create apass and sloan catalogs](#Create-apass-and-sloan-catalogs-for-new-objects)
* [Running the pipeline](#Running-the-pipeline)
* [Cookbook](#Cookbook)
* [Creating a Landolt Catalog](#Creating-a-Landolt-Catalog)
* [Difference Imaging](#Difference-Imaging)
* [Definitions](#Definitions)
    * [Telescopes](##Telescopes)
    * [Filetype](##Filetype)
    * [Status](##Status)
    * [Database](##Database)

# Ingest data
* Download observations from SNEx and run ingesttar.py. This data will have already been reduced with BANZAI.
* Download observations from archive.lco.global using LCOGTingest.py. It is recommended that you download data reduced with BANZAI. This will only work for public data (e.g., standard stars images and observations from over 1 year ago) unless you have an account for the observing portal. For example (all arguments optional; see ```LCOGTingest.py -h``` for options):
```
LCOGTingest.py -n NAME -s YYYY-MM-DD -e YYYY-MM-DD -t EXPOSE -r reduced --public
```

# Create apass and sloan catalogs for new objects
* run `comparecatalogs.py` to generate new catalogs
* Note if you are trying to reduce U band, you need to generate a local catalog. See [Creating an Landolt Catalog](#Creating-a-Landolt-Catalog) for details.

# Cookbook
## Basic reduction
This is a description of the stream-lined steps that are recommended for processing most data. This cookbook assumes that you have downloaded reduced data and have created catalogs for each object as described above.

1. Verify your PSF catalog. It is recommended that the `psf` stage be run with the catalog option pointing to one of the catalogs you downloaded (usually apass or sloan) so that the same stars are used for build the PSF for every image. You should run the `psf` stage using the `--show` option on a few images in each filter to verify that the star selected from the catalog are good (e.g. not double stars, near a hot pixel, or saturated). In addition to showing you the model PSF, if a DS9 window is open prior to running, the image will be shown with the stars using to build the PSF circled in cyan. Note that these ID number represent the line number (after the header) in the catalog file and not the catalog ID numbers. Remove any stars that are being selected and are not good from the catalog in `$LCOSNPIPE/trunk/src/lsc/standards/cat/...`. Additionally, you copy the file to the installation directory. The easiest way to do this is to re-install the pipelines `python setup.py install`. Repeat this stage until you are satisfied with the stars selected for the catalog.  
**Example:**
    ```
    ds9&
    lscloop.py -n 2016cok -e 20160528-20160530 -s psf --show --catalog $LCOSNDIR/standard/cat/apass/AT2016cok_apass.cat
    ```

2. Run the `psf` stage on all of the data without the `--show` option  
**Example**
    ```
    lscloop.py -n 2016cok -e 20160528-20180104 -s psf --catalog $LCOSNDIR/standard/cat/apass/AT2016cok_apass.cat
    ```

3. Check the PSF and image quality. Run the `checkpsf` stage with the `--show`. I like to do this by filter so I can develop a good baseline of what a field should look like in a given filter. Inspect each PSF and each image. Mark the PSFs that you want to redo with `n` and the images that you never want to see again as `b`. At this stage the only PSFs that I've redone are ones with exceptional seeing, where stars that normally wouldn't were saturated. See https://www.overleaf.com/read/sccbqgnhwyfh for a description of different PSF shapes you may see.  
**Example for checking r-band images**
    ```
    ds9&
    lscloop.py -n 2016cok -e 20160528-20180104 -s checkpsf --show -f B
    ```
4. Redo the PSF for any images you marked with `n`. These are selected using the `-b psf` option.  
**Example** 
lscloop.py -n 2016cok -e 20160528-20180104 -s psf -b psf --catalog $LCOSNDIR/standard/cat/apass/AT2016cok_apass.cat

5. Calculate instrumental magnitudes by running the `psfmag` stage. This will derive aperture and psf photometry.   
**Example**
```
lscloop.py -n 2016cok -e 20160528-20180104 -s psfmag
```

6. Calculate the photometric zeropoint for each image. The `zcat` stage should be run for at least two filters at a time so that the color term can be evaluated.  
**Example**
```
lscloop.py -n 2016cok -e 20160528-20180104 -s zcat
```

7. Calculate the apparent magnitude for an image using the instrumental magnitude derived from PSF photometry using the `--type fit` option. This is another stage that should be run on more than one filter so the color term can be used.  
**Example**
```
lscloop.py -n 2016cok -e 20160528-20180104 -s mag --type fit
```

8. Visually inspect your light curve using `getmag` stage with the `--show` option. This will bring up a plot with your light curve (I like to do this one filter at a time). You can click on individual points to bring up a second window which shows a cut out of the supernova on the left and the residual after the PSF subtraction on the right. You are given the option to remove bad points from the light curve. It is recommended that you use the `n` option at this stage, allowing you to easily redo these observations from any stage using the `-b mag` option. In general, however, this is used to eliminate points with a lot of scatter due to poor observing conditions.


# Running the pipeline:
In general, the pipeline is called as a series of stages. Each stage should be called with the general format:
```
lscloop.py -e YYYYMMDD-YYYYMMDD -n NAME -s STAGE
```
Where:
* ```-e YYYYMMDD-YYYYMMDD``` is the date range of the observations you want to process
* ```-n NAME``` is the name of the object you want to process
* ```-s STAGE``` is the stage which is detailed below, in recommended order

## Broadly used options
* ```-F``` force a step to be run again, even if it succeeded last time
* ```-f FILTER``` run only on observations from one filter or set of filters (U,u,B,g,V,r,R,i,I,z,w,landolt, apass, sloan)
* ```--id``` run only on a specific image specified by a 3 digit number in the filename. For example you would use ```--id 046``` to run on the file elp1m008-fl05-20180302-0046-s91.
* ```-T``` run only on observations from one telescope. Valid options: 1m0, 2m0, or 0m4
* ```-I``` run only on observations from one instrument. Valid options: kb, fl, fs, sinistro, sbig
* ```-b STAGE``` run only on observations marked at bad at a given stage (where stage is quality, wcs, psf, psfmag, zcat)

## Steps:
### Cosmic ray rejection:
* **Call**: ```-s cosmic```
* **Description**: Runs lacosmic to clean the images. This stage only needs to be run prior to the `-s diff` stage as that is the only step that uses the cosmic ray mask
* **Recommended Options**:
* **Other Options**: 

### WCS solution
* **Call**: ```-s wcs```
* **Description**: Generate a WCS solution for a given observation. In general, this stage does not need to be run when starting analysis from observations that have been reduced with Banzai (which is recommended) as Banzai solves the WCS solution.
* **Recommended Options**: 

* **Other Options**: 
    * ```--catalog=CATALOG```: where CATALOG is either usnoa2, usnob1, 2mass, or apass. If this option is not specified then the following catalogs are tried: ['2mass', 'usnoa2', 'usnob1']
* **How to tell if this step worked**:
    * If the pipeline thinks it was successful, then the column wcs in the photlco table of the database will have a 0 in it, if it wasn't able to find a solution the value in this column will be 9999
    * Additionally, when the pipeline fails on this step it sets the quality from 127 to 1. 
    * To inspect images for which the wcs step failed:
        * open DS9
        * run the standard lscloop.py with stage ```-s checkquality -b quality --show```
        * You will then be asked if the image is good or bad. If you say its good, the quality will be set to 127, if bad, the quality will stay 1
    * To inspect the WCS solution for images for which the wcs step succeeded:
        * open DS9
        * run the standard lscloop.py with the stage ```-s checkwcs -b wcs --show```
        * This displays image for which the wcs value is 9999 in DS9 with the catalog of stars on top of it. You should see the stars in the image at the locations marked by the catalog.
        * You will then be asked if this image is good with the option to answer: yes, no, or bad. Here if you answer bad, the quality will be set from 127 to 1, indicating that the image should not be used. If you mark no, then the wcs will be reset, indicating that you would like to try to calibrate the wcs again.
        * if you run this step without ```-b wcs``` then you can inspect all observations
* **How to fix this step if the image looks ok but the wcs failed**
    * rerun lscloop.py with the step ```-s wcs -b wcs``` to run only on those files for which the wcs failed
    * try running with the --catalog option
    * try using ```--use_sextractor``` flag TODO: describe this

### Generate a PSF model
* **Call**: ```-s psf```
* **Description**: creates a model psf for each image
* **Recommended Options**:
    * ```--catalog=CATALOG``` this uses the stars found at the location of stars in the catalog to generate a model PSF (TODO: check that this is true)
* **Other Options**: 
* **How to tell if this step worked**:
    * This step fills in the psf column with a filename (for the psf model) if it succeeds and an X if it fails.
    * To inspect images for which the pipeline failed:
        * open DS9
        * run the standard lscloop.py with the stage ```-s checkpsf -b psf --show``` TODO: check if you need --show or -i
        * This displays image for which the psf value is X in DS9 and in a separate window shows the psf model. The psf model should be a single source (e.g. not a blended gaussian). TODO: is it ok if its elongated?
        * You will then be asked if this image is good with the option to answer: yes, no, or bad. Here if you answer bad, the quality will be set from 127 to 1, indicating that the image should not be used. If you mark no, then the psf will be reset, indicating that you would like to try to calibrate the psf again.
        * if you run this step without ```-b psf``` then you can inspect all observations
        * BETA: ```-s checkpsf``` has the option to use matplotlib to display the psf rather than iraf with the option ```--no_iraf```. In this case, you can manipulate the figure (rotate and pan) to better see the psf
* **How to fix this step if the image looks ok but the psf step failed**
    * Try running without the  --catalog option
    * Try running with the ```--fwhm``` option. This is especially true for the 0.4m telescopes. Possible values to try: 5 and 7
    * 

### Generate instrumental magnitudes using psf photometry
* **Call**: ```-s psfmag```
* **Description**: Generate instrumental magnitudes using psf photometry
* **Recommended Options**:
* **Other Options**: 
* **How to tell if this step worked**:
    * This step sets the psfmag column to the instrumental magnitude if this step works and to 9999 if this step fails.
* **How to fix this step if the image looks ok but the psfmag step failed**

### Calculate the zeropoint
* **Call**: ```-s zcat```
* **Description**: Calculate the zeropoint of an image using the apparent magnitude of stars in the image as found in the catalog provided. This will be used to convert instrumental magnitudes to apparent magnitudes
* TODO: What does the pipeline default to when no catalog is provided?
* **Recommended Options**:
    * ```--catalog=CATALOG``` where CATALOG should be a catalog that corresponds to the filter you are trying to calibrate (e.g. apass, sloan, or landolt). The pipeline will use the location and apparent magnitude of the stars in this catalog to calculate the conversion between instrumental and apparent magnitude.
* **Other Options**: 
* **How to tell if this step worked**:
    * This set fills in the zcat column in the photlco database with a file name if it succeeded and an X if it failed
    * The most likely failure mode is that there are not enough stars between the catalog and the image. You can visually check this by:
        * Openning DS9
        * running ```-s zcat -i```
        * if not enough stars are found, then "no calibration: FILTER" is printed to the screen and only the last file is displayed in DS9. You can run this file by file using the --id option described under [Broadly used options](#Broadly-used-options)
        * if enough stars are found, this also displays a plot calibrated_mag - instrumental_mag vs color.  Each point is a star being used to calculate the zero point. Green are accepted, red are excluded by some automatic sigma-clipping. You can toggle each point by clicking on it if you don’t like the automatic rejection.
* **How to fix this step if the image looks ok but the zcat step failed**

### Generate instrumental magnitudes using psf photometry
* **Call**: ```-s mag```
* **Description**: Generate instrumental magnitudes using psf photometry
* **Recommended Options**:
    * ```--type fit```: this generates an apparent magnitude from the instrumental magnitude derived from psf photometry.
* **Other Options**: 
    * ```--type ph```: generate apparent magnitudes from aperture photometry
    * ```--type mag```: TODO: what does this mean?
* **How to tell if this step worked**:
    * If this step worked, an apparent magnitude will be put in the mag column of the photlco table. If this step failed this field will be set to 9999
* **How to fix this step if the image looks ok but the psfmag step failed**

### Evaluating Image quality from lightcurve outliers:
* **Call**: ```-s getmag``` TODO: check if you need --show or -i
* **Description**: view whole light curve and inspect cutouts of object and residuals after the psf subtraction. 
* **Recommended Options**:
* **Other Options**:
    * ```-f FILTER``` where FILTER is U,B,g,V,r,R,i,I, landolt, apass, or sloan. This displays light curve for only this filter or set of filters
    * If a point looks funny, then you can click on it. This brings up an image with a cutout of the image, and a cut out of the residual after the psf has been subtracted. It will ask you if you want to make this image with yes, no, or bad (TODO: verify this and fill in what happens if you answer no vs bad)


#
* Check whether there are enough stars in common between a field and the catalog being used in the zcat step ```-s zcat -i```


# Creating a Landolt Catalog:
### Download landolt standard star catalogs
* These need to be obtained from LCO (Jamie gave me a zip file with them)
* Put these files in the directory $LCOSNDIR/standard/cat/landolt
* Reinstall the pipeline ```python setup.py install```

### Download and Calibrate the Standard Star Observations
* Find standard stars that were observed during the epoch you need in the filter you need (U). Ideally, you also obtain at least 1 other filter (for color correction) - recommended B and V
* Search archive.lco.global for obstype=STANDARD, set date range, filter, if applicable telescope and site (e.g. fo 2018zd I set site=elp and telescope=fl05)
* You will be using standards from the same telescope and observed on the same night as your observations
* Download the reduced observations
* run python ingestzip.py -f downloaded_directory_or_zip_filename
* Use the standard star image to find the nightly zeropoint for the nights on which you observed your target and the standard (run cosmic, wcs, psf, psfmag, zcat)
    * When running the zcat step, use ```--catalog=$LCOSNDIR/standard/cat/landolt/<standard catalog name>```
* When you run -s local on the SN (in the next step), it queries the database for any standards in the same filter-site-night as the SN observations. You can check that these are identified correctly with ```lscloop.py -n 'SN 2018zd' -e 20180302-20180330 -f landolt -T 1m0 --standard STANDARD``` where STANDARD is the name of your standard (e.g. L95)

### Create a catalog of stars (local sequence) in SN field for Landolt filters
* This calculates the apparent magnitude for the stars in a given filter in your SN field and creates a catalog that will be used in the zcat step to generate the zeropoint for each observation. If a single class of telescopes seems to be an outlier, then you should limit the telescopes used to calculate the apparent magnitude to a single class of telescopes with the ```-T 1m0``` flag
```
lscloop.py -n 'SN 2018zd' -e 20180302-20180330 -f landolt -s local -i --standard STANDARD --catalog=CATALOG
```
Where STANDARD is the standard you calibrated previously and CATALOG is the name of the apass or sloan catalog created for your supernova (e.g. SN2018zd_apass.cat). You should use the full epoch so that you calculate the apparent magnitudes for a variety of night (ideally) making your values more robust to outliers.

* Check plots for:
    - if there is a big difference between the BV filters of Landolt and APASS for a large number of stars (a few stars will be variable)
    - I guess the main thing is that some of the standard field observations will be bad because of clouds or whatever, so make sure they are not pulling the median zero point away too much
* This outputs a catalog (the name of which is printed to the screen)

### Move the catalog to the catalog directory:
* Move the catalog created in the previous step to the directory $LCOSNPIPE/trunk/src/lsc/standard/cat/landolt
* Reinstall the pipeline ```python setup.py install```

### Manually update the database
* update the targets table of the database (manually) to know about this catalog ```update targets set landolt_cat=CATALOG where id=ID``` where CATALOG is the catalog you moved in the previous step and ID is the id of your supernova
### Run the photometry as usual 
* e.g. ```lscloop.py -n 'SN 2018zd' -e 20180302-20190419 -f landolt -s zcat```

# Difference Imaging
* An outline of the steps to perform difference imaging can be found here: https://www.authorea.com/users/75900/articles/96044-image-subtraction-with-lcogtsnpipe
* When running the `zcat` stage, `--field` is used to select the catalog to which you will calibrate the magnitudes of those images. If you do not give `--field`, it uses whatever you put in `-f` by default (which is correct for Sloan, for example).
* When the `-s diff --difftype 1` fails with the error message `ERROR:root:Too few stars in common at 5-sigma; lower and try again` this could be the result of bad WCS, bad background subtraction, bad cosmic ray mask, bad quality image. Regardless, most of these issues are solved by adding the `--unmask` flag to ignore the bad pixel mask when looking for stars
* An experimental feature has been added under the `--no_iraf` flag. This will use the Python package `reproject` (specifically the `reproject_interp` function) to align the images instead of using IRAF's `gregister`. We added this after experiencing problems where `gregister` does not align images correctly, but we have not looked deeply into the root cause of that problem. We think that `reproject_interp` does not use the same algorithm as `gregister`, and in fact it may use a worse algorithm, so this option should be used with caution. Note that IRAF is still used for `--fixpix` regardless of the `--no_iraf` flag.

# Definitions
## Telescopes
| Short Name   | Long Name                         | ?             | lscloop keyword|
|--------------|-----------------------------------|---------------|----------------|
| OGG 2m       | Haleakala Observatory - 2m        | ogg2m001-fs02 | 2m0            |
| COJ 2m       | Siding Springs Observatory - 2m   | coj2m002-fs03 | 2m0            |
| COJ 1m       | Siding Springs Observatory - 1m   | coj1m003-kb71 | 1m0            |
| LSC 1m       | CTIO - Region IV                  | lsc1m004-kb77 | 1m0            |
| LSC 1m       | CTIO - Region IV                  | lsc1m005-kb78 | 1m0            |
| ELP 1m       | McDonald Observatory - 1m         | elp1m008-kb74 | 1m0            |
| LSC 1m       | CTIO - Region IV                  | lsc1m009-fl03 | 1m0            |
| CPT 1m       | SAAO - Sutherland Facilities - 1m | cpt1m010      | 1m0            |
| COJ 1m       | Siding Springs Observatory - 1m   | coj1m011-kb05 | 1m0            |
| CPT 1m       | SAAO - Sutherland Facilities - 1m | cpt1m012-kb75 | 1m0            |
| CPT 1m       | SAAO - Sutherland Facilities - 1m | cpt1m013-kb76 | 1m0            |

## Filetype:
| Number | Meaning         |
|--------|-----------------|
|1       | Original image  |
|2       | Merged image    |
|3       | Difference image|
|4       | Reference image |

sn2.fits is a table with aperture and psf photometry measurements created by the psf stage
## Status
| Status | Meaning |
|--------|----------|
| -4     | bad quality |
| -3     | image not ingested in the data table |
| -2     | image not in working directory |
| -1     | sn2.fits (fits table with all psf and aperture measurements created by the psfmag step) not in working directory |
| 0      | not done and previous stage not done |
| 1      | not done and possible since previous stage is done |
| 2      | done and possible to do again |
| 3      | local sequence catalog available |

## Database
### apercorr column in photlco
The aperture correction column as was added to the photlco table with [PR 52](https://github.com/LCOGT/lcogtsnpipe/pull/52); see also [PR 59](https://github.com/LCOGT/lcogtsnpipe/pull/59) which fixes a number of bugs related to the aperture correction . Note that despite its name, this correction has a different definition from what is usually referred to as an aperture correction. During the PSF stage, both aperture and PSF photometry is performed on the catalog stars. The aperture photometry calculated with an aperture that is int(fwhm * 3 + 0.5) is compared to the PSF photometry for all of the catalog stars. The assumption is that the aperture contains all of the stellar flux and the aperture photometry therefore represents the true photometry. On the otherhand, small differences between the model and true PSF may lead to errors in the PSF photometry. For this reason, the aperture correction is calculated as the sigma-clipped mean difference between the aperture and PSF photometry and the error is the sigma-clipped standard deviation for the catalog stars. If the aperture correction exceeds a pre-defined threshold, the psf stage will fail under the assumption that a large aperture correction indicates a poor PSF model.

The aperture correction is applied to the PSF photometry in the sn2 catalog file (which will be used to calculate the zeropoint during the zcat stage) and to the supernova photometry during the psfmag stage. For this reason, for standard PSF photometry on filetype 1 objects when the PSF photometry is used to calculate the zeropoint, the apparent magnitudes is independent of the aperture correction as the aperture correction that is applied during the psfmag stage is undone by the zeropoint. 

This situation is more complicated for difference imaging as for PyZOGY the PSF is a model that neither represents the original nor the template PSF and for HOTPANTS the PSF can either be from the template or the original image. For this reason and the fact that aperture photometry should be sufficient for difference imaging, the ability to run PSF photometry on difference imaging has been disabled. However, should it ever be enabled again, the following decisions have been made. When the PSF stage is run for PyZOGY, an "image" of the model PSF is created and the PSF stage is run as if that is the image with a catalog of one star centered on the middle of the frame. In this case, the aperture correction is meaningless and the pipeline sets the aperture correction to 0. For HOTPANTS, the CONVOL00 keyword is read from the header to determine which PSF was used for difference imaging. The aperture correction is then read from the sn2 file corresponding to the original image or the template image. This is because the sn2 file from the difference image is copied from the original or the template based on the normalization keyword and not the PSF keyword.