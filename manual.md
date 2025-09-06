# Table of Contents
* [Ingest Data](#Ingest-data)
* [Create apass and sloan catalogs](#Create-apass-and-sloan-catalogs-for-new-objects)
* [Cookbook](#Cookbook)
* [Running the pipeline](#Running-the-pipeline)
* [Creating a Landolt Catalog](#Creating-a-Landolt-Catalog)
* [Difference Imaging](#Difference-Imaging)
* [Tips and Tricks](#Tips-and-Tricks)
* [Common Errors](#Common-Errors)
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
If you wish to download data that is not public yet, request an API key in the GSP slack. Someone at LCO will send it to you privately and you can use it in the following manner (replacing `sdf087asdf987asdf0786987` with your API key):
```
export LCO_API_KEY=sdf087asdf987asdf0786987
LCOGTingest.py -n NAME -s YYYY-MM-DD -e YYYY-MM-DD -t EXPOSE -r reduced
```

# Create gaia, apass, and sloan catalogs for new objects
* run `comparecatalogs.py` to generate new catalogs
* Note if you are trying to reduce U band, you need to generate a local catalog. See [Creating an Landolt Catalog](#Creating-a-Landolt-Catalog) for details.
* As an experimental feature, you can now download Pan-STARRS1 3π catalogs instead of SDSS catalogs for calibrating _griz_ images.
  To use this feature run `comparecatalogs.py -p`. Note that you will not be able to calibrate _u_-band images. This will only work when `sloan_cat=NULL` in the database (for the initial run). If the pipeline has already determined that a sloan catalog doesn't exist you will have to reset this field manually in the database.

# Cookbook
## Basic reduction
This is a description of the stream-lined steps that are recommended for processing most data. This cookbook assumes that you have downloaded reduced data and have created catalogs for each object as described above.

1. Run the `psf` stage on all of the data without the `--show` option. The default option is to use the Gaia catalog unless another catalog is specified with the `--catalog` keyword.
**Example**
    ```
    lscloop.py -n 2016cok -e 20160528-20180104 -s psf 
    ```

3. Check the PSF and image quality. Run the `checkpsf` stage with the `--show`. I like to do this by filter so I can develop a good baseline of what a field should look like in a given filter. Inspect each PSF and each image, the stars used to build are marked in blue if you use the `--no_iraf` flag. Mark the PSFs that you want to redo with `n` and the images that you never want to see again as `b` (this sets the quality to 1, you have to run the `checkquality -b quality` stage to recover these files). See https://www.overleaf.com/read/sccbqgnhwyfh for a description of different PSF shapes you may see. In addition to the PSF, you should look at which stars are identified in the image as these will be used for the flux calibration later.
**Example for checking r-band images**
    ```
    ds9&
    lscloop.py -n 2016cok -e 20160528-20180104 -s checkpsf --show -f B --no_iraf
    ```
4. Redo the PSF for any images you marked with `n`. These are selected using the `-b psf` option. Suggestions for improving the PSFs: change the fwhm (`--fwhm` to be bigger or smaller (depending on which stars you want selected), increase the `--nstars` parameter (e.g. `--nstars 12`) to increase the number of stars that get averaged together to make the PSF; adjust the `--datamax` and `--datamin` parameters to exclude bright stars or cosmic rays (e.g. `--datamax 60000 --datamin 0`); use a different catalog (as described above). You can play around with different FWHM values using the `--show` option and saying the PSF is bad, you will then be prompted to enter a new PSF.
**Example** 
    ```
    lscloop.py -n 2016cok -e 20160528-20180104 -s psf -b psf 
    ```

5. Calculate instrumental magnitudes by running the `psfmag` stage. This will derive aperture and psf photometry. If the aperture photometry is failing often, you may need to set `--datamax` and `--datamin`.  
**Example**
    ```
    lscloop.py -n 2016cok -e 20160528-20180104 -s psfmag
    ```

6. Calculate the photometric zeropoint for each image. The `zcat` stage should be run for at least two filters at a time so that the color term can be evaluated. You can specify which catalog of field stars to use as calibration with the `--field` argument.
**Example**
    ```
    lscloop.py -n 2016cok -e 20160528-20180104 -s zcat -f landolt --field apass
    lscloop.py -n 2016cok -e 20160528-20180104 -s zcat -f sloan --field sloan
    ```

7. Calculate the apparent magnitude for an image using the instrumental magnitude derived from PSF photometry using the `--type fit` option. This is another stage that should be run on more than one filter so the color term can be used.  
**Example**
    ```
    lscloop.py -n 2016cok -e 20160528-20180104 -s mag --type fit
    ```

8. Visually inspect your light curve using `getmag` stage with the `--show` option. This will bring up a plot with your light curve (I like to do this one filter at a time). You can click on individual points to bring up a second window which shows a cut out of the supernova on the left and the residual after the PSF subtraction on the right. You are given the option to remove bad points from the light curve. It is recommended that you use the `d` option at this stage, allowing you to easily redo these observations from any stage using the `-b mag` option. In general, however, this is used to eliminate points with a lot of scatter due to poor observing conditions. By supplying `--output output_filename.csv` at this stage, the output can be written to a file. When photometry is final, use the `--uploadtosnex2` flag to send your light curve to SNEx 2.
**Example**
    ```
    lscloop.py -n 2016cok -e 20160528-20180104 -s getmag --type mag --show
    lscloop.py -n 2016cok -e 20160528-20180104 -s getmag --type mag --output output_filename.csv
    ```


# Running the pipeline:
In general, the pipeline is called as a series of stages. Each stage should be called with the general format:
```
lscloop.py -e YYYYMMDD-YYYYMMDD -n NAME -s STAGE
```
Where:
* ```-e YYYYMMDD-YYYYMMDD``` is the date range of the observations you want to process, this can be a single date if so desired
* ```-n NAME``` is the name of the object you want to process
* ```-s STAGE``` is the stage which is detailed below, in recommended order

## Broadly used options
* ```-F``` force a step to be run again, even if it succeeded last time
* ```-f FILTER``` run only on observations from one filter or set of filters (U,u,B,g,V,r,R,i,I,z,w,landolt, apass, sloan)
* ```--id, -d``` run only on a specific image specified by a 3 digit number in the filename. For example you would use ```--id 046``` to run on the file elp1m008-fl05-20180302-0046-s91.
* ```-T``` run only on observations from one telescope. Valid options: 1m0, 2m0, or 0m4. Note: because of the implementation, there is a bug/feature to `-T` that you can search for any substring in the filename
* ```-I``` run only on observations from one instrument. Valid options: kb, fl, fs, sinistro, sbig, muscat, ep, sq, qhy
* ```-b STAGE``` run only on observations marked at bad at a given stage (where stage is quality, wcs, psf, psfmag, zcat, mag)
* in many of the quality check stages the user is asked whether a file is good or bad and there are two options for bad, `b` and `n`. In general the `b` option should only be used for unusable images as it removes the observation completely from future processing and hides it. The only way to get this observation back is to run `-s checkquality -b quality`. If you answer `n` that stage of the pipeline will be reset for you to run again (as if the stage failed).

## Steps:
### Cosmic ray rejection:
* **Call**: ```-s cosmic```
* **Description**: Runs lacosmic to clean the images. This stage only needs to be run prior to the `-s diff` stage as that is the only step that uses the cosmic ray mask
* **Recommended Options**:
* **Other Options**: 

### WCS solution
* **Call**: ```-s wcs```
* **Description**: Generate a WCS solution for a given observation. In general, this stage does not need to be run when starting analysis from observations that have been reduced with Banzai (which is recommended) as Banzai solves the WCS solution (files processed by Banzai end in `.e91.fits`). If Banzai fails to solve the WCS solution, then the WCS column will be populated with 4 - this means you may need to run the WCS stage manually. If the file was marked with `n` in the `checkwcs` stage or previously marked with `b` and resurrected, the value will be 9999.0. The easiest way to fix this is the run the `wcs` stage with `-s wcs -b wcs`.
* **Recommended Options**: 

* **Other Options**: 
    * ```--catalog=CATALOG```: where CATALOG is either usnoa2, usnob1, 2mass, or apass. If this option is not specified then the following catalogs are tried: ['2mass', 'usnoa2', 'usnob1']
* **How to tell if this step worked**:
    * If the pipeline thinks it was successful, then the column wcs in the photlco table of the database will have a 0 in it, if it wasn't able to find a solution the value in this column will be 9999
    * Alternately, run ```lscloop.py -n <snname> -e <epoch> -b wcs``` and any other criteria you'd like to filter on to get a list of files that failed the WCS stage
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
    * try using ```--mode astrometry``` to use astrometry.net to solve the WCS. Note this only works if astrometry.net has been installed. astrometry.net is not installed in the docker version of the pipeline
    * try using ```--xshift 1 --yshift 1``` to start the solution from a slightly different location
    * If you get the error: "ERROR: catalog empty 2MASS" this means `lscastrodef.querycatalogue` is failing. It could mean one of two things:
        - the vizier query in `lscastrodef.vizq` is timing out. If this is the issue, you will either have to wait or manually solve the WCS with astrometry.net and then update the WCS column in the photlco table
        - the `CRPIX` value in the header of the file you downloaded from SNEx was corrupt. You can check this by looking at the `CRPIX1` and `CRPIX2` header keywords. These should not be very large numbers (like 1E42). If the file is corrupt, you can download it directly from the LCO archive or if you don't have permission, request in the pipeline channel that someone download the files for you (or for large downloads, give you the proper archive permissions).
        - For either case, you should also check that the apass, sloan, and gaia catalogs you downloaded for the target are not empty

### Generate a PSF model
* **Call**: ```-s psf```
* **Description**: creates a model psf for each image
* **Recommended Options**:
* **Other Options**: 
    * ```--catalog=CATALOG``` this uses the stars found at the location of stars in the catalog to generate a model PSF
    * ```--max_apercorr``` if the median difference between the psf and aperture magnitudes for catalog stars is more than this value, the PSF stage will fail. In practice, if this value is large its an indication that the model PSF is not characterizing the actual PSF well.
* **How to tell if this step worked**:
    * This step fills in the psf column with a filename (for the psf model) if it succeeds and an X if it fails.
    * Alternately, run ```lscloop.py -n <snname> -e <epoch> -b psf``` and any other criteria you'd like to filter on to get a list of files that failed the PSF stage
    * To inspect all images and their PSFs:
        * run the standard lscloop.py with the stage ```-s checkpsf --show --no_iraf``` (if running without the `--no_iraf` flag you should open DS9 first)
        * This displays an image with the stars in the catalog marked and the stars used to build the PSF marked in a different color. In a separate window the psf model is shown. The psf model should be a single source although a non-gaussian PSF may still be an accurate representation of the PSF for the image. The stars chosen to build the PSF should not be so bright that they are saturated.
        * In addition to inspecting the PSF, you should look at which stars are identified in the image as these will be used for the flux calibration later. The stars should be evenly distributed across the image and the bright stars should be selected.
        * You will then be asked if this image is good with the option to answer: yes, no, or bad. Here if you answer bad, the quality will be set from 127 to 1, indicating that the image should not be used and it will not be included in any future processing. If you mark no, then the psf will be reset, indicating that you would like to try to calibrate the psf again.

* **How to fix this step if the image looks ok but the psf step failed**
    * Try running with the `--fwhm` option. This is especially true for the 0.4m telescopes. Possible values to try: 5 and 7
    * Try adjusting `--datamin` or `--datamax` to select different stars. You can use the datamin and datamax output during the PSF stage to change which stars are selected
    * Try running with the  `--field` option. For some sources, the pipeline may have trouble picking suitable stars using the default 'gaia' catalog. So, try running `-s psf --field apass` or `-s psf --field sloan`.
    * Try running with `--use-sextractor` argument.
    * `--field` and `--use-sextractor` commands can be combined with `--datamin` and `--datamax` to select suitable stars.

### Generate instrumental magnitudes using psf and aperture photometry
* **Call**: ```-s psfmag```
* **Description**: Generate instrumental magnitudes using psf and aperture photometry
* **Recommended Options**:
* **Other Options**: 
    * -i run this stage interactively
    * `--datamin` and `--datamax` manually set the datamin and datamax keywords used by iraf to select the PSF stars. These may be necessary to use if the aperture photometry is failing and you're interested in it
* **How to tell if this step worked**:
    * This step sets the psfmag and apmag column to the instrumental magnitude if this step works and to 9999 if this step fails.
    * Alternately, run ```lscloop.py -n <snname> -e <epoch> -b psfmag``` (or ```-b apmag```) and any other criteria you'd like to filter on to get a list of files that failed the PSF stage
* **How to fix this step if the image looks ok but the psfmag step failed**
    * Try setting `--datamin` and/or `--datamax`

### Calculate the zeropoint
* **Call**: ```-s zcat```
* **Description**: Calculate the zeropoint of an image using the apparent magnitude of stars in the image as found in the catalog provided. This will be used to convert instrumental magnitudes to apparent magnitudes
* **Recommended Options**:
    * ```-f sloan```, ```-f apass```, or ```-f landolt```: more than one filter should be provided so that the zeropoint color term can be calculated
* **Other Options**: 
    * ```--catalog=CATALOG``` where CATALOG should be a catalog that corresponds to the filter you are trying to calibrate (e.g. apass, sloan, or landolt). The pipeline will use the location and apparent magnitude of the stars in this catalog to calculate the conversion between instrumental and apparent magnitude.
    * --unfix allow the slope of the color term to vary (use with extreme caution)
* **How to tell if this step worked**:
    * This set fills in the zcat column in the photlco database with a file name if it succeeded and an X if it failed
    * Alternately, run ```lscloop.py -n <snname> -e <epoch> -b zcat``` and any other criteria you'd like to filter on to get a list of files that failed the zcat stage
    * The most likely failure mode is that there are not enough stars between the catalog and the image. You can visually check this by:
        * Openning DS9
        * running ```-s zcat -i```
        * if not enough stars are found, then "no calibration: FILTER" is printed to the screen and only the last file is displayed in DS9. You can run this file by file using the --id option described under [Broadly used options](#Broadly-used-options)
        * if enough stars are found, this also displays a plot calibrated_mag - instrumental_mag vs color.  Each point is a star being used to calculate the zero point. Green are accepted, red are excluded by some automatic sigma-clipping. You can toggle each point by clicking on it if you don’t like the automatic rejection.
* **How to fix this step if the image looks ok but the zcat step failed**

### Generate apparent magnitudes using from instrumental magnitudes
* **Call**: ```-s mag```
* **Description**: Generate apparent magnitudes using instrumental magnitudes
* **Recommended Options**:
   *  ```-f sloan```, ```-f apass```, or ```-f landolt```: more than one filter should be provided to use the zeropoint calculated with a color term. If a single filter is given the mag column will be populated with 9999.
* **Other Options**: 
    * The default for non-difference images is to use PSF instrumental magnitudes for this stage (```--type fit```), however, for difference imaging,  aperture instrumental magnitudes (```--type ph```) are used for consistency between methods (because the aperture correction is applied to HOTPANTS subtractions but not to PyZOGY subtractions). This may cause additional scatter in your light curve.
    * ```--type ph```: generate apparent magnitudes from aperture photometry
    * ```--type fit```: generate apparent magnitudes from PSF photometry
* **How to tell if this step worked**:
    * If this step worked, an apparent magnitude will be put in the mag column of the photlco table. If this step failed this field will be set to 9999 (if every field is 9999 check that you asked for more than one filter).
    * Alternately, run ```lscloop.py -n <snname> -e <epoch> -b mag``` and any other criteria you'd like to filter on to get a list of files that failed the PSF stage
* **How to fix this step if the image looks ok but the mag step failed**
    * If all mag values are 9999, it is likely that you ran this stage with a single filter. More than one filter should be provided to use the zeropoint calculated with a color term. Try running with one of the recommended options for ```-f```.
    * Since this is just a combination of the output of the psfmag stage and the zcat stage, inspect the outputs of these stages and redo as needed

### Evaluating Image quality from lightcurve outliers:
* **Call**: ```-s getmag --show```
* **Description**: view whole light curve and inspect cutouts of object and residuals after the psf subtraction. 
* **Recommended Options**:
* **Other Options**:
    * ```-f FILTER``` where FILTER is U,B,g,V,r,R,i,I, landolt, apass, or sloan. This displays light curve for only this filter or set of filters
    * ```-o, --output <filename>``` write the output to a csv file
    * ```--updatetosnex2``` upload output to SNEx 2
    * If a point looks funny, then you can click on it. This brings up an image with a cutout of the image, and a cut out of the residual after the psf has been subtracted. It will ask you if you want to make this image with yes, delete (d), upper limit (u), or bad. In theory, yes should keep the observation and all derived quantities, delete (d) removes everything back to the psf stage, upper limit (u) set the magtype to -1, and bad (b) sets the quality to 1. It is recommended that you use delete to remove an observation rather than b. Currently, there is a bug which sets images marked as bad to upper limits (i.e. b is aliased to u; see issue #56).
    
# Creating a Landolt Catalog:
### Download landolt standard star catalogs
* Landolt standard star catalogs (L*.cat) are now part of the pipeline repository and can be found in $LCOSNDIR/standard/cat/landolt. 
* For SA (selected area) standards, you may need to create your own catalog. First check Table 2 in [Landolt, 2009](https://ui.adsabs.harvard.edu/abs/2009AJ....137.4186L). This has proper coordinates although, as with all of the catalogs, the colors will need to be converted into magnitudes to match the other landolt catalogs. If the field is not in that paper, you can search [Vizier](https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=II/183A/table2&-out.max=1000&-out.all=2&-out.all=2&-out.form=HTML%20Table&-oc.form=sexa). **Warning**: The coordinates in the original (Landolt, 1992) published paper have no decimal places, which is not precise enough for our needs. You must manually fix this! You may be able to look up more precise coordinates in [Simbad](https://simbad.cds.unistra.fr/simbad/sim-fid) by searching by identifier, e.g., "SA 109*" (don't forget the wildcard to get all the stars in the field). However, this table may not contain uncertainties on the magnitudes, so you might have to manually assemble the final catalog. We then encourage you to then contribute the corrected catalog back to the pipeline for others to use.
   * Put these files in the directory $LCOSNDIR/standard/cat/landolt
   * Reinstall the pipeline ```python setup.py install```

### Download and Calibrate the Standard Star Observations
* Find standard stars that were observed during the epoch you need in the filter you need (U). The standard reduction practice is to create a landolt catalog for U, B, and V filters and calibrate your data with these. At a minimum, you also obtain B band for the color correction. [SNEx 2](https://test.supernova.exchange) now has a link to a list all standard star observations that match your supernova observations on the photometry tab between the list of photometry and the plot of photometry.
* Search archive.lco.global for observation type "STANDARD" and/or proposal "standard" (see note below), set date range, filter, site, telescope, and instrument. If you can't find a standard on the same night, the next best option is to all for different instruments (but same site and night). If you do this, you will need to use the `--match-by-site` option when you run the `local` stage.
* You will be using standards from the same telescope and observed on the same night as your observations
* Use LCOGTingest.py to download and ingest observations. For example ``` LCOGTingest.py -S lsc -T 1m0a -f U -I fl05 -s 2019-09-26 -e 2019-10-08  -r reduced -t STANDARD --public``` Where `-S` is the site, `-T` is the telescope, `-f` is the filter, -I is the instrument, `-s` is the start date, `-e` is the end date, `-r reduced` downloads files that have already been reduced with BANZAI, and `--public` is always required. Note: the date format has dashes and the search is exclusive of the end data. Also note: for many years in the late 2010s, LCO scheduled some (but not all) standard field observations with observation type EXPOSE instead of STANDARD. You can still identify them by their proposal ID "standard" (i.e., in the command above, replace `-t STANDARD` with `-P standard`).
 Also also note: you cannot search the LCO archive by telescope number (e.g., 1m0-08), but you can search by instrument number (e.g., fl05), which is equivalent for a given night (at least). The image filenames contain both the telescope number and the instrument number.
* If you are downloading a standard for the first time, manually update the targets table so that the standards you downloaded have `classificationid=1` and to set the catalog to the one you downloaded earlier.
    * ```mysql -u supernova -D supernova -p``` (if you are using docker-compose you also need ```-h supernovadb```)
    * Get the targetid and landolt catalog of the standard you want to update: ```SELECT targetnames.targetid, name, classificationid, landolt_cat FROM targets JOIN targetnames ON targets.id=targetnames.targetid WHERE name="L107"``` (you should replace L107 with the name of your standard)
    * Check that you're selecting the right row: ```SELECT targetnames.targetid, name, classificationid, landolt_cat FROM targets JOIN targetnames ON targets.id=targetnames.targetid WHERE targetid=55``` (you should replace 55 with the targetid you found in the last step)
    * If `classificationid` is not 1, update `classificationid` for the row selected: ```UPDATE targets SET classificationid=1 WHERE id=<TARGETID>``` (where <TARGETID> is replaced with the targetid of your standard)
   * If `landolt_cat` is empty, you should update it with the name of the catalog you downloaded (e.g. L107_landolt_stetson.cat). ``` UPDATE targets SET landolt_cat=<LANDOLT_FILENAME> WHERE targetid=<TARGETID>``` (where <TARGETID> is replaced with the targetid of your standard and <LANDOLT_FILENAME> is the name of your downloaded landolt catalog)
* When you run a command with `--standard all`, it queries the database for any standards in the same filter-telescope-night as the SN observations and runs the stage on those images instead. You can check that the standard images are identified correctly with ```lscloop.py -n 'SN 2018zd' -e 20180302-20180330 -f landolt -T 1m0 --standard all```. If you want to calibrate to a single standard field, you can also give, e.g., `--standard L94`. If a photometric night with a supernova observation and standard star taken on the same telescope with the same instrument does not exist. You should request observations of your SN field with the landolt filters until such a night exists. There is no requirement that your SN be visible at this time. If you get the error `NameError: name 'refcat' is not defined` then you need `comparecatalogs.py` to download catalogs for your supernova.
* Create photometry catalogs for the standard star images. You need these to calculate the zero points.
```
lscloop.py -n 'SN 2018zd' -e 20180302-20180330 -f landolt --standard all -s psf
```
* Calculate the zero points for the standard star images. You should run this interactively to make sure there are enough stars identified in the image.
```
lscloop.py -n 'SN 2018zd' -e 20180302-20180330 -f landolt --standard all -s zcat -i
```

### Create a catalog of stars (local sequence) in SN field for Landolt filters
* This calculates the apparent magnitude for the stars in a given filter in your SN field and creates a catalog that will be used in the zcat step to generate the zeropoint for each observation. If a single class of telescopes seems to be an outlier, then you should limit the telescopes used to calculate the apparent magnitude to a single class of telescopes with the ```-T 1m0``` flag
```
lscloop.py -n 'SN 2018zd' -e 20180302-20180330 -f landolt -s local -i --standard all
```
Where STANDARD is the standard you calibrated previously. You should use the full epoch so that you calculate the apparent magnitudes for a variety of night (ideally) making your values more robust to outliers.

* Check plots for:
    - if there is a big difference between the BV filters of Landolt and APASS for a large number of stars (a few stars will be variable)
    - I guess the main thing is that some of the standard field observations will be bad because of clouds or whatever, so make sure they are not pulling the median zero point away too much
* This outputs a catalog (the name of which is printed to the screen)
* Your catalog needs at least B-band magnitudes, in addition to U, to perform color corrections. If there were no B-band standards taken for your dates/site, you'll need to manually add the B-band magnitudes of the stars into your new catalog. The id's should match those in the APASS catalog for your SN, so you can add the appropriate values to the B and Berr columns in your new landolt catalog before proceeding. 

### Move the catalog to the catalog directory:
* Move the catalog created in the previous step to the directory $LCOSNPIPE/trunk/src/lsc/standard/cat/landolt
* Reinstall the pipeline ```python setup.py install```

### Manually update the database
* update the targets table of the database (manually) to know about this catalog ```update targets set landolt_cat=CATALOG where id=ID``` where CATALOG is the catalog you moved in the previous step and ID is the id of your supernova
### Run the photometry as usual 
* e.g. ```lscloop.py -n 'SN 2018zd' -e 20180302-20190419 -f landolt -s zcat```
* note: because you are calculating the conversion from instrumental to apparent magnitudes, you do not need to run any stage before the zcat stage

# Difference Imaging
Start with the automatically reduced photometry. I will use NAME as the
name of the object, TARGDATE-TARGDATE as the range of dates for the
science images (images where the target was visible), and TEMPDATE as
the epoch(s) of the reference images. ID refers to the ID number of a
specific file. IN refers to the instrument prefix (kb for SBIG, fl for
Sinistro, fs for Spectral). TEMPTEL is the instrument prefix for the
template image (kb, fl, fs, fa, or SDSS). If you don't give this option, it
will select the same instrument given with -T.

## Prepare the Reference Images

### LCOGT References

1.  Look through your reference images and choose the best one for each camera--filter combination. Make a note of their ID numbers.
    ```
    ds9 &
    ```
    ```
    lscloop.py -n NAME -e TEMPDATE -s checkwcs
    ```
    If wcs is off by a little bit, rerun the wcs stage with the following command, then run `checkwcs`:
    ```
    lscloop.py -n NAME -e TEMPDATE -s wcs --mode astrometry -F
    ```

2.  Mark your chosen images as references.
    ```
    lscloop.py -n NAME -e TEMPDATE -d ID -s template
    ```

3.  Run cosmic ray rejection on the reference images.
    ```
    lscloop.py -n NAME -e TEMPDATE --filetype 4 -s cosmic
    ```

4.  Generate PSFs for the reference images. I recommend doing this interactively for each image (using `--show`), since these PSFs will affect all of your subtractions.
    ```
    lscloop.py -n NAME -e TEMPDATE --filetype 4 -s psf --show
    ```

### SDSS References


1.  Choose a set of science images that includes one image with each camera--filter combination used. Then run the following command once for each of those images, using the ID numbers to choose individual frames. This will take a while.
    ```
    lscloop.py -n NAME -e TARGDATE -d ID -s ingestsloan
    ````

    Make a note of the TEMPDATE for each SDSS frame you download.

2.  Generate PSFs for the SDSS images. I usually get good results using a FWHM of $\sim 5$ pixels for SDSS.
    ```
    ds9 &\
    ```
    ```
    lscloop.py -n NAME -e TEMPDATE --filetype 4 -s psf --show --fwhm 5 --use-sextractor
    ```

3.  Copy the variance images to the right place for use as cosmic ray masks.
    ```
    lscloop.py -n NAME -e TEMPDATE --filetype 4 -s cosmic
    ```

### PAN-STARRS (PS1) References


1.  Choose a set of science images that includes one image with each camera--filter combination used. Then run the following command once for each of those images, using the ID numbers to choose individual frames. This will take a while.
    ```
    lscloop.py -n NAME -e TARGDATE -d ID -s ingestps1
    ````

    Make a note of the TEMPDATE for each PS1 frame you download.

2.  Generate PSFs for the PS1 images.
    ```
    ds9 &\
    ```
    ```
    lscloop.py -n NAME -e TEMPDATE --filetype 4 -s psf --show --fwhm 5 --use-sextractor
    ```

3.  Copy the variance images to the right place for use as cosmic ray masks.
    ```
    lscloop.py -n NAME -e TEMPDATE --filetype 4 -s cosmic
    ```

## Do the Image Subtraction

1.  Run cosmic ray rejection on all the science images. This will take a while.
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE -s cosmic
    ```

2.  In the mean time, check to make sure all the science images have good PSFs from the automatic pipeline.
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE -s checkpsf
    ```

    Generate a new PSF for any images you marked as bad. I usually get good results using a FWHM of $\sim 7$ pixels, but this depends on the seeing.
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE -b psf -s psf --show --fwhm 7
    ```

    If the images are saturated, adjust the Dmax value. For my images 75000 worked, but depends on the gain (saturation is 65535$*$gain):
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE -b psf -s psf --show --fwhm 7 --datamax 75000
    ```

3.  Once all the cosmic ray rejection is done (for science and reference images), run the subtraction. This will take a while. By default, `--tempdate=19990101-20080101` (useful for SDSS, PS1, and any other archival template image), `--temptel=IN`, `--fixpix=False`, and `--difftype=0` (0 = HOTPANTS, 1 = Optimal). For archival templates, the telescope name (`-T`) must be `kb, fs, fl, fa, ep,` or `sq` and `--temptel` is `PS1` or `SDSS`
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --normalize t -T IN [--tempdate TEMPDATE] [--temptel TEMPTEL] [--fixpix] [--difftype 1] -s diff
    ```

    If you want, look over the results. Make sure to choose Frame $>$Tile on DS9.
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s checkdiff
    ```

    * When the `-s diff --difftype 1` fails with the error message `ERROR:root:Too few stars in common at 5-sigma; lower and try again` this could be the result of bad WCS, bad background subtraction, bad cosmic ray mask, bad quality image. Regardless, most of these issues are solved by adding the `--unmask` flag to ignore the bad pixel mask when looking for stars

    * An experimental feature has been added under the `--no_iraf` flag. This will use the Python package `reproject` (specifically the `reproject_interp` function) to align the images instead of using IRAF's `gregister`. We added this after experiencing problems where `gregister` does not align images correctly, but we have not looked deeply into the root cause of that problem. We think that `reproject_interp` does not use the same algorithm as `gregister`, and in fact it may use a worse algorithm, so this option should be used with caution. Note that IRAF is still used for `--fixpix` regardless of the `--no_iraf` flag.

## Extract the Photometry


By default, `--filetype 3` selects both HOTPANTS and Optimally subtracted images. To select just one or the other, use `--difftype 0` for HOTPANTS or `--difftype 1` for Optimal (PyZOGY).

1.  First copy the convolved PSFs over.
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s psf
    ```

2.  Then calculate the instrumental magnitudes. Use a flat background for subtracted images.
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -x 1 -y 1 -s psfmag
    ```

    If the `apmags` are all bad, then there was probably a bad pixel near the supernova messing it up. To fix this, use `--datamax` and `--datamin`.
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s psfmag --datamax 100000 --datamin -1000 -F
    ```

    If you want, look over the results and redo any bad fits.

    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s checkmag
    ```

3.  Calculate the zero point and color term for each set of filters. (If the field is not in SDSS, use `--field apass` for the Sloan filters too.)
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -f landolt --field apass -s zcat
    ```
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -f sloan -s zcat
    ```
    When running the `zcat` stage, `--field` is used to select the catalog to which you will calibrate the magnitudes of those images. If you do not give `--field`, it uses whatever you put in `-f` by default (which is correct for Sloan, for example).

4.  Calculate the apparent magnitudes. Make sure you don't separate the filters if you want to include color corrections.
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE --filetype 3 -s mag
    ```

5.  If the lightcurve has outliers or faint data, run the getmag stage to interactively check out suspicious points. Can do `--type ph` or `--type fit` to check uncalibrated magnitudes (aperture and psf, respectively).
    ```
    lscloop.py -n NAME -e TARGDATE-TARGDATE -s getmag --type mag --show \[-f FILTER\] \[--filetype 3\]
    ```











# Tips and Tricks
* Cleanly stopping a pipeline process: use ctrl+Z to suspend the process, and then kill %1 (or whatever the job number is). Sometimes this works better than ctrl+C when there are multiple processes
   
* To check which standard stars you have downloaded for a given SN use the command `lscloop.py -n SN_NAME -e SN_EPOCH --standard all`

# Common Errors and how to fix them
* `sn2.fits file not found`: this usually indicates that the PSF stage was not run or ran and failed. You can check this by running lscloop.py with your target and epoch but no stage and checking the PSF column (X means it failed, filename means it succeeded). The most commonly missed error in the PSF stage is "The difference between the aperture and psf magnitudes exceeds the minimum allowed".

* `The difference between the aperture and psf magnitudes exceeds the minimum allowed  ## > 0.1` : The aperture correction compares the median aperture and psf photometry of field stars in your image. It was decided that if this differs by more than ~10% (0.1 mag) then something is likely wrong with your PSF (or possibly your aperture) and this deserves a closer look. First, run the psf stage with `--show` or the checkpsf stage with `--no_iraf` and check your PSF stars are good (not double, not near the edge, not saturated, not under diffraction spikes, etc). If your PSF stars look good, you can usually fix this with one of the following:
    - increase the `fwhm`. This value should be printed to the screen during the PSF stage. You can suggest a FWHM to lscloop with the keyword `--fwhm` (e.g. `--fwhm 7`)
    - Remove the brightest star by setting the `--datamax` keyword to a value lower than the first PSF star datamax, which can be read from the output
    - Increase the number of PSF stars using the `--nstars` flag (e.g. `--nstars 12`)
    - Try using different `--field` for selecting psf stars e.g. `-s psf --field apass` or `-s psf --field sloan`.
    - Try running with `--use-sextractor` argument.
    - `--field` and `--use-sextractor` commands can be combined with `--datamin` and `--datamax` to select suitable stars.
   
* I accidentally marked a file as `b` how do I bring it back? -OR- One of my files disappeared, how do I get it back?
   Files marked with `b` are considered failed observations that no amount of finessing will bring back (e.g. no signal). For this reason, the pipeline hides them from all future processing. The only way to get this observation back is to run `-s checkquality -b quality`. If you answer `n` that stage of the pipeline will be reset for you to run again (as if the stage failed).
   
* My instrumental aperture magnitudes for different imaging are all 9999. Try running the `psfmag` stage with `--datamin -1000 --datamax 100000`
   
* All of my apparent magnitudes are 9999 after I run the mag stage. This occurs when you run the mag stage with only one filter because the color correction can't be applied. Try running with `-f sloan` or `-f landolt`
   
# Definitions
## Telescopes
Current facilities can be found at: https://lco.global/observatory/sites/  
Telescope numbers are unique within a given mirror size

| Short Name   | Long Name                         | Telescope id (-T)  | Telescope Number | Site id (-S)|
|--------------|-----------------------------------|----------------|------------------|----------------|
| OGG 2m       | Haleakala Observatory - 2m        | 2m0            | 01               | ogg            |
| COJ 2m       | Siding Springs Observatory - 2m   | 2m0            | 02               | coj            |
|--------------|-----------------------------------|----------------|------------------|----------------|
| COJ 1m       | Siding Springs Observatory - 1m   | 1m0            | 03               | coj            |
| LSC 1m       | CTIO - Region IV                  | 1m0            | 04               | lsc            |
| LSC 1m       | CTIO - Region IV                  | 1m0            | 05               | lsc            |
| ELP 1m       | McDonald Observatory - 1m         | 1m0            | 06               | elp            |
| ELP 1m       | McDonald Observatory - 1m         | 1m0            | 08               | elp            |
| LSC 1m       | CTIO - Region IV                  | 1m0            | 09               | lsc            |
| CPT 1m       | SAAO - Sutherland Facilities - 1m | 1m0            | 10               | cpt            |
| COJ 1m       | Siding Springs Observatory - 1m   | 1m0            | 11               | coj            |
| CPT 1m       | SAAO - Sutherland Facilities - 1m | 1m0            | 12               | cpt            |
| CPT 1m       | SAAO - Sutherland Facilities - 1m | 1m0            | 13               | cpt            |
|--------------|-----------------------------------|----------------|------------------|----------------|
| COJ 0.4m     | Siding Springs Observatory - 0.4m | 0m4            | 03               | coj            |
| OGG 0.4m     | Haleakala Observatory - 0.4m      | 0m4            | 04               | ogg            |
| COJ 0.4m     | Siding Springs Observatory - 0.4m | 0m4            | 05               | coj            |
| OGG 0.4m     | Haleakala Observatory - 0.4m      | 0m4            | 06               | ogg            |
| CPT 0.4m     | SAAO - Sutherland Facilities -0.4m| 0m4            | 07               | cpt            |
| LSC 0.4m     | CTIO - Region IV                  | 0m4            | 09               | lsc            |
| TFN 0.4m     | Teide Observatory 0.4m            | 0m4            | 10               | tfn            |
| ELP 0.4m     | McDonald Observatory - 0.4m       | 0m4            | 11               | elp            |
| LSC 0.4m     | CTIO - Region IV                  | 0m4            | 12               | lsc            |
| TFN 0.4m     | Teide Observatory 0.4m            | 0m4            | 14               | tfn            |

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
