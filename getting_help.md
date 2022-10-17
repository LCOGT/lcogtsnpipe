# How to get LCO SN Photometry Pipeline Help
`lcogtsnpipe` is not currently maintained by anyone and any improvements and help offered to the GSP is done on a volunteer basis. We are really excited to see the reduction of GSP data occurring outside of LCO and will work to enable this for you. To ensure the efficient use of everyone's time and better assist you and the broader community, we kindly ask that you follow these guidelines:

* Review the LCOGT documentation:
    * Check the Cookbook section of the [manual](https://github.com/LCOGT/lcogtsnpipe/blob/master/manual.md)
    * Check the Running the Pipeline, Creating a Landolt Catalog, and Difference Imaging sections of the [manual](https://github.com/LCOGT/lcogtsnpipe/blob/master/manual.md) (depending on what you're doing)

* Check out the [tutorial](https://docs.google.com/document/d/14ADvdbS-19flwtU7TRJ1lyx8IjgBK6D1KubZnQkQciA/edit) for an example of a full reduction

* Confirm you are using the most recent version of the pipeline
    * you can check this by going to your pipeline directory and running `git pull origin master`. If you are up to date you will get a message that ends with `Already up to date`
    * If you do not get this message, run step 8 of the [installation instructions](https://github.com/LCOGT/lcogtsnpipe#readme) (Build Docker), then follow the `Use` instructions (also part of the installation instructions)

* Confirm that if you updated the code your stages are actually running (either you removed entries from the database or you are using the -F flag to force reprocessing)

* Read through all of the output (I know, its a lot). 
    * Make sure no errors or warnings anywhere in the output. The pipeline is made to run in batch mode and not stop all processing when a single image fails somewhere in the reduction. This means that it may appear to finish correctly but may have actually errored earlier.
    * Check that the error you are looking at is the first error. Sometimes IRAF issues an error and then the pipeline wrapper issues a secondary error because it didn't get the output it expected, which is not as informative as the original error. 
    * The first table contains a tremendous amount of information about what stages have and have not been run. Make sure that the status of the file is what you expect it to be. 

* Use the Slack search bar to search the pipeline channel for similar issues

* Quality check your results. If the code is running but the results aren't what you expect, use the check stages to ensure that your observation is what you expect it to be (i.e. its not blank, the PSF is well characterized, the difference images are well subtracted, etc.)

* Check with your research group and advisor. Even if they are unfamiliar with the pipeline - a fresh set of eyes, especially ones that have done photometry before can be really helpful at spotting reduction issues

If the above does not solve your issue:
* Request assistance on Slack:
    * Post in the most relevant #pipeline channel
    * Confirm that you have followed the above guidelines
    * Describe the documentation you were following (as applicable)
    * Share the full command you are running
    * If there was a crash, paste in a substantial portion of your trace (this goes back to the multiple errors/errors farther up issue)

* Request assistance via office hours:
    * office hours will be held weekly at 11am pacific time on Tuesdays. Zoom link in the GSP Pipeline Slack channel description.
    * be prepared to screenshare

* `lcogtsnpipe` is currently developed and maintained by just a handful of GSP volunteers at present. Your contributions to this open source project would be greatly appreciated. You can contribute by: 
    * forking the code and improving the documentation. 
    * forking the code and improving the code/algorithms.  
    * Your changes will be reviewed and added to the code base via a github pull request. 


*These guidelines inspired and drawn from the PypeIt project*
