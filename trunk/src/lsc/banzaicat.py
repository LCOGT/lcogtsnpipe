import lsc
from astropy.io import fits
from scipy.stats import sigmaclip
from operator import itemgetter

def make_cat(filename,datamax=75000,banzai_sig=3,banzai_crlim=3):

	if datamax == None: datamax = 75000

	hdul = fits.open(filename)
	banzai_cat = hdul[1].data

	print "Total number of sources in BANZAI catalog: {0}".format(len(banzai_cat))

	ellipticities = [x['ELLIPTICITY'] for x in banzai_cat]
	backgrounds = [x['BACKGROUND'] for x in banzai_cat]
	fwhms = [x['FWHM'] for x in banzai_cat]

	filtered_el, lo, hi = sigmaclip(ellipticities, low=banzai_sig, high=banzai_sig)
	filtered_bg, lo, hi = sigmaclip(backgrounds, low=banzai_sig, high=banzai_sig)
	filtered_fwhm, lo, hi = sigmaclip(fwhms, low=banzai_sig, high=banzai_sig)

	id_num = 0
	sources = []

	for source in banzai_cat:
		if (source['FLAG'] == 0 
				and source['PEAK'] <= datamax
				and source['ELLIPTICITY'] in filtered_el 
				and source['BACKGROUND'] in filtered_bg
				and source['FWHM'] in filtered_fwhm 
				and source['FWHM'] > banzai_crlim):
			id_num += 1
			
			StN = source['PEAK']/source['BACKGROUND']	

			sources.append([source['RA'],source['DEC'],StN,id_num])

	print ("Number of sources in BANZAI catalog after filtering: "
		"{0}".format(len(sources)))
	print ("({0}-sigma clipping on source ellipticity, "
		"background level, and FWHM.)".format(banzai_sig))

	#Sort by S/N	
	sources = sorted(sources, key=itemgetter(2), reverse=True)

	header = "# BEGIN CATALOG HEADER\n"
	header += "# nfields 13\n"
	header += "#     ra     1  0 d degrees %10.5f\n"
	header += "#     dec    2  0 d degrees %10.5f\n"
	header += "#     id     3  0 c INDEF %15s\n"
	header += "# END CATALOG HEADER\n"
	header += "#\n"

	with open('banzai.cat','w') as banzai_cat_file:

		banzai_cat_file.write(header)
		for source in sources:
			line = "{0:10.5f}\t{1:10.5f}\t{2}\n".format(source[0],source[1],source[3])
			banzai_cat_file.write(line)

	print "Saving the {0} best sources to banzai.cat".format(len(sources))

	hdul.close()
	return 'banzai.cat'
