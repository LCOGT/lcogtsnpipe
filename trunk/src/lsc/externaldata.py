from astropy.io import fits
import os
import numpy as np
import lsc
import datetime

def jd2date(inputjd):
 jd0 = 2451544.5 # On Jan 1, 2000 00:00:00
 return datetime.datetime(2000,01,01,00,00,00)+datetime.timedelta(days=inputjd-jd0)


def MJDnow(datenow='',verbose=False):
   import time
   _JD0=55927.
   if not datenow:
      datenow=datetime.datetime(time.gmtime().tm_year, time.gmtime().tm_mon, time.gmtime().tm_mday, time.gmtime().tm_hour, time.gmtime().tm_min, time.gmtime().tm_sec)
   _JDtoday=_JD0+(datenow-datetime.datetime(2012, 01, 01,00,00,00)).seconds/(3600.*24)+\
             (datenow-datetime.datetime(2012, 01, 01,00,00,00)).days
   if verbose: print 'JD= '+str(_JDtoday)
   return _JDtoday


def SDSS_gain_dark(camcol, ugriz, run):
    if camcol == 1:
        if ugriz == 'u':
            gain = 1.62
            dark = 9.61
        elif ugriz == 'g':
            gain = 3.32
            dark = 15.6025
        elif ugriz == 'r':
            gain = 4.71
            dark = 1.8225
        elif ugriz == 'i':
            gain = 5.165
            dark = 7.84
        elif ugriz == 'z':
            gain = 4.745
            dark = 0.81
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'
    elif camcol == 2:
        if ugriz == 'u':
            if run < 1100:
                gain = 1.595
            elif run > 1100:
                gain = 1.825
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'
            dark = 12.6025
        elif ugriz == 'g':
                gain = 3.855
                dark = 1.44
        elif ugriz == 'r':
                gain = 4.6
                dark = 1.00
        elif ugriz == 'i':
            gain = 6.565
            if run < 1500:
                dark = 5.76
            elif run > 1500:
                dark = 6.25
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'
        elif ugriz == 'z':
            gain = 5.155
            dark = 1.0
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'
    elif camcol == 3:
        if ugriz == 'u':
            gain = 1.59
            dark = 8.7025
        elif ugriz == 'g':
            gain = 3.845
            dark = 1.3225
        elif ugriz == 'r':
            gain =  4.72
            dark = 1.3225
        elif ugriz == 'i':
            gain = 4.86
            dark = 4.6225
        elif ugriz == 'z':
            gain = 4.885
            dark = 1.0
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'
    elif camcol == 4:
        if ugriz == 'u':
            gain = 1.6
            dark = 12.6025
        elif ugriz == 'g':
            gain = 3.995
            dark = 1.96
        elif ugriz == 'r':
            gain =  4.76
            dark = 1.3225
        elif ugriz == 'i':
            gain = 4.885
            if run < 1500:
                dark = 6.25
            elif run > 1500:
                dark = 7.5625
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'
        elif ugriz == 'z':
            gain = 4.775
            if run < 1500:
                dark = 9.61
            elif run > 1500:
                dark = 12.6025
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'
    elif camcol == 5:
        if ugriz == 'u':
            gain = 1.47
            dark = 9.3025
        elif ugriz == 'g':
            gain = 4.05
            dark = 1.1025
        elif ugriz == 'r':
            gain = 4.725
            dark = 0.81
        elif ugriz == 'i':
            gain = 4.64
            dark = 7.84
        elif ugriz == 'z':
            gain = 3.48
            if run < 1500:
                dark = 1.8225
            elif run > 1500:
                dark = 2.1025
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'
    elif camcol == 6:
        if ugriz == 'u':
            gain = 2.17
            dark = 7.0225
        elif ugriz == 'g':
            gain = 4.035
            dark = 1.8225
        elif ugriz == 'r':
            gain = 4.895
            dark = 0.9025
        elif ugriz == 'i':
            gain = 4.76
            dark = 5.0625
        elif ugriz == 'z':
            gain = 4.69
            dark = 1.21
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'
    else:
        print 'ERROR in SDSS_dark_gain: CAMCOL is not set!'
    return gain, dark

def downloadsdss(_ra,_dec,_band,_radius=20, force=False):
    from astroquery.sdss import SDSS
    from astropy import coordinates as coords
    import astropy.units as u
    from astropy.io import fits
    import os
    import sys
    import string
    import numpy as np
    from scipy import interpolate
    pos = coords.SkyCoord(ra=float(_ra)*u.deg,dec=float(_dec)*u.deg)
    print 'pos =', pos
    xid = SDSS.query_region(pos, spectro=False, radius=_radius*u.arcsec)
    print xid
    if xid:
       pointing=[]
       for i in xid:
          if i['run'] > 300:
             if (i['run'],i['camcol'],i['field']) not in pointing:
                pointing.append((i['run'],i['camcol'],i['field']))
       # if too many pointing, take only first 40
       if len(pointing) > 50:
          nn=50
       else: 
          nn=len(pointing)
       filevec=[]
       print len(pointing)
       for run, camcol, field in pointing[:nn]:
          #  naomaggie image
          output1 = _band+'_SDSS_'+str(run)+'_'+str(camcol)+'_'+str(field)+'.fits'
          #  image in count
          output2 = _band+'_SDSS_'+str(run)+'_'+str(camcol)+'_'+str(field)+'c.fits'
          #  weight image
          output3 = _band+'_SDSS_'+str(run)+'_'+str(camcol)+'_'+str(field)+'.weight.fits'
          #  sky image
          output4 = _band+'_SDSS_'+str(run)+'_'+str(camcol)+'_'+str(field)+'.sky.fits'
          if os.path.isfile(output1) and not force:
              print 'already downloaded', output1
              filevec.append(output2)
              filevec.append(output3)
              continue
          im = SDSS.get_images(run=run, camcol=camcol, field=field, band=_band, cache=True)
          if os.path.isfile(output1):
             os.system('rm '+output1)
          if os.path.isfile(output2):
             os.system('rm '+output2)
          if os.path.isfile(output3):
             os.system('rm '+output3)
          if os.path.isfile(output4):
             os.system('rm '+output4)
          im[0].writeto(output1)
#         im[0][0].writeto(output2)

          FITS_file = fits.open(output1)
          new_header = FITS_file[0].header
          camcol     = FITS_file[0].header['CAMCOL']  # camcol
          ugriz      = FITS_file[0].header['FILTER']  # ugriz filter
          run1        = FITS_file[0].header['RUN']     # run
          gain, dark_var = SDSS_gain_dark(camcol, ugriz, run1)
          new_header['gain']  = gain
          new_header['dark']  = dark_var
          new_header['BUNIT']  = 'counts'
          new_header['rdnoise']  = 2
          frame_image = FITS_file[0].data.transpose()
          allsky     = FITS_file[2].data['ALLSKY'].transpose()
          allsky     = allsky[:,:,0]
          xinterp    = FITS_file[2].data['XINTERP'].transpose()
          xinterp    = xinterp[:,0]
          yinterp    = FITS_file[2].data['YINTERP'].transpose()
          yinterp    = yinterp[:,0]
          sky_function = interpolate.interp2d(np.arange(allsky.shape[1]),\
                                              np.arange(allsky.shape[0]), allsky, kind='linear')
          sky_image    = sky_function(yinterp, xinterp) # in counts
          calib     = FITS_file[1].data #  nanomaggies per count
          calib_image = np.empty_like(frame_image)
          for i in np.arange(calib_image.shape[1]):
             calib_image[:,i] = calib
          # Calculate the error in the frame image for use fitting algorithms later.
          dn_image        = frame_image / calib_image  # + sky_image # counts
          dn_err_image    = np.sqrt((dn_image + sky_image)/ gain + dark_var)
          # frame_image_err = dn_err_image * calib_image  converts to nanomaggies
          
          frame_weight = 1 / ((dn_err_image)**2)
          new_header['SKYLEVEL']  = np.mean(sky_image)
          #  save image in count
          fits.writeto(output2, dn_image.transpose(), new_header,overwrite=True)
          #  save weight image
          fits.writeto(output3, frame_weight.transpose(), new_header,overwrite=True)
          #  save sky image 
          fits.writeto(output4, sky_image.transpose(), new_header,overwrite=True)
          filevec.append(output2)
          filevec.append(output3)
       return filevec
    else:
       return ''

def sdss_swarp(imglist,_telescope='spectral',_ra='',_dec='',output='', objname='',survey='sloan',combine_type='MEDIAN', show=False):
    import re
    import datetime
    import lsc
    import time
    import string

    if _telescope == 'spectral':
        pixelscale = 0.30104  # 2 meter
        _imagesize =  2020
    elif _telescope == 'sbig':
            pixelscale = 0.467  # 1 meter
            _imagesize =  2030
    elif _telescope == 'sinistro':
            pixelscale = 0.387  # 1 meter
            _imagesize =  4020

    if survey =='sloan':
       print('no este')
       out1 = 'SDSS'
       hdr = fits.getheader(imglist[0])
       _filter = hdr.get('filter')
       _gain   = hdr.get('gain')
       _ron   = hdr.get('rdnoise')
       if not _ra:
          _ra = hdr.get('CRVAL1')
       if not _dec:
          _dec = hdr.get('CRVAL2')
       _saturate = hdr.get('SATURATE')
       if not _saturate:
          _saturate = 61000
       if 'day-obs' in hdr:
          _dayobs = hdr.get('day-obs')
       elif 'date-obs' in hdr:
          _dayobs = re.sub('-','',hdr.get('date-obs'))

       if '/' in _dayobs:
          _dayobs = '19'+''.join(string.split(_dayobs,'/')[::-1])

       try:
          _mjd = MJDnow(datetime.datetime(int(str(_dayobs)[0:4]),int(str(_dayobs)[4:6]),int(str(_dayobs)[6:8])))
       except:
          print 'warning, no mjd'
          _mjd = 0
          _dayobs = '19991231'

    elif survey =='ps1':
       out1 = 'PS1'
       imglist2=[]
       for i,img in enumerate(imglist):
           print i,img
           hdr = fits.open(img)
           if os.path.isfile(re.sub('.fits','_1.fits',img)):
              os.system('rm '+re.sub('.fits','_1.fits',img))
           fits.writeto(re.sub('.fits','_1.fits',img),hdr[1].data,hdr[1].header)
           imglist2.append(re.sub('.fits','_1.fits',img))

       imglist = imglist2
       hdr = fits.getheader(imglist[0])
       _saturate = hdr.get('HIERARCH CELL.SATURATION')
       _filter = string.split(hdr.get('HIERARCH FPA.FILTER'),'.')[0]
       _gain   = hdr.get('HIERARCH CELL.GAIN')
       _ron   = hdr.get('HIERARCH CELL.READNOISE')
       _mjd = hdr.get('MJD-OBS')
       _dayobs = jd2date(_mjd+2400000.5).strftime('%Y%d%m')
       if not _ra:
          _ra = hdr.get('RA_DEG')
       if not _dec:
          _dec = hdr.get('DEC_DEG')

    elif survey =='decam':
       print('yes este')
       out1 = 'decam'
       imglist2=[]
       for i,img in enumerate(imglist):
           print i,img
           hdr = fits.open(img)
           if os.path.isfile(re.sub('.fits','_1.fits',img)):
              os.system('rm '+re.sub('.fits','_1.fits',img))
           fits.writeto(re.sub('.fits','_1.fits',img),hdr[0].data,hdr[0].header)
           imglist2.append(re.sub('.fits','_1.fits',img))

       imglist = imglist2
       hdr = fits.getheader(imglist[0])
       #_saturate = hdr.get('HIERARCH CELL.SATURATION')
       #_gain   = hdr.get('HIERARCH CELL.GAIN')
       #_ron   = hdr.get('HIERARCH CELL.READNOISE')
       #_mjd = hdr.get('MJD-OBS')
       #_gain = 3.855 # set a random number for now to see if it works
       _filter = string.split(hdr.get('FILTER'),' ')[0]
       _dayobs = string.split(hdr.get('DATEOBS2'),'T')[0]# DATE-OBS for the last  image in the stack
       _dayobs = _dayobs.replace('-', '')
       #_telescope = string.split(hdr.get('TELESCOP'),' ')[0]
       _saturate = 61000 # used for sloan
       #_ron = 2 # from the code downloadsdss not for decam
       _gain = hdr.get('ARAWGAIN')
       #_saturate = hdr.get('SATURATE')
       _ron = hdr.get('RDNOISEA')
       _exp = 1800
       if not _ra:
          _ra = hdr.get('CRVAL1')
       if not _dec:
          _dec = hdr.get('CRVAL2')


    #  airmass
    if 'airmass' in hdr:
        _airmass = hdr.get('airmass')
    else:
        _airmass = 1

    #  filters
    filt={'U':'U','B':'B','V':'V','R':'R','I':'I','u':'up','g':'gp','r':'rp','i':'ip','z':'zs'}
    if _filter in filt.keys():
        _filter = filt[_filter]

    #if 'date-obs' in hdr:
    #   _dateobs = hdr.get('date-obs')
    #else:
    #   _dateobs = jd2date(_mjd+2400000.5).strftime('%Y-%d-%m')
    
    if not(_dayobs):
       _dateobs = jd2date(_mjd+2400000.5).strftime('%Y-%d-%m')

    if not output:
       output = _telescope+'_'+str(out1)+'_'+str(_dayobs)+'_'+str(_filter)+'_'+objname +'.fits'

    imgmask = [] 
    skylevel = []
    # if the template is from sloan, use mask and weight image
    if survey =='sloan':
       welist = [i for i in imglist if '.weight.' in i]
       if len(welist):
          # if weight images provided, use them in swarp
          # and take out them from input list image
          imglist = [j for j in imglist if j not in welist]
          imgmask = welist
          for jj in imglist:
             img_data,img_header=fits.getdata(jj, header=True)
             if 'SKYLEVEL' in img_header:
                skylevel.append(img_header['SKYLEVEL'])
    elif survey=='ps1':
       wtlist = [i for i in imglist if '.wt_' in i]
       mklist = [i for i in imglist if '.mk_' in i]
       if len(wtlist):
          for i,name in enumerate(wtlist):
             weightimg = re.sub('.wt_','.weight.',name)
             imgmask.append(weightimg)
             if re.sub('.wt_','.mk_',name) in mklist:
                weight_data, weight_header = fits.getdata(name, header=True)
                mask_data, mask_header = fits.getdata(re.sub('.wt_','.mk_',name), header=True)
                weight_data = 1/weight_data 
                weight_fits = fits.PrimaryHDU(header=weight_header, data=weight_data)
                weight_fits.writeto(weightimg, output_verify='fix', overwrite=True)
             else:
                os.system('cp '+name+' '+weightimg)
       imglist = [j for j in imglist if (j not in wtlist) and (j not in mklist)]
       # measure skylevel in all ps1 images
       for jj in imglist:
          img_data,img_header=fits.getdata(jj, header=True)
          skylevel.append(np.mean(img_data))

    elif survey=='decam':
       wtlist = [i for i in imglist if '-invvar' in i]
       mklist = [i for i in imglist if '-maskbits' in i]
       if len(wtlist):
          print('wtlist', wtlist)
          for i,name in enumerate(wtlist):
             weightimg = re.sub('-invvar','.weight.',name)
             imgmask.append(weightimg)
             if re.sub('-invvar','-maskbits',name) in mklist:
                weight_data, weight_header = fits.getdata(name, header=True)
                mask_data, mask_header = fits.getdata(re.sub('-invvar','-maskbits',name), header=True)
                weight_data = 1/weight_data
                weight_fits = fits.PrimaryHDU(header=weight_header, data=weight_data)
                weight_fits.writeto(weightimg, output_verify='fix', overwrite=True)
             else:
                os.system('cp '+name+' '+weightimg)
       imglist = [j for j in imglist if (j not in wtlist) and (j not in mklist)]
       # measure skylevel in all ps1 images
       for jj in imglist: 
          img_data,img_header=fits.getdata(jj, header=True)
          #skylevel.append(np.mean(img_data))
          if 'SKYLEVEL' in img_header:
            skylevel.append(img_header['SKYLEVEL'])
          else:
            skylevel.append(np.mean(img_data))


#    elif len(welist):
#       imglist = [j for j in imglist if j not in welist]
#       imgmask = welist


    print(imgmask, 'imgmask')
    print skylevel
    if survey == 'sloan' or survey == 'ps1':
        print('swarping')
        sampling = 'BILINEAR' # LANCZOS3
        line = 'swarp ' + ','.join(imglist) + ' -IMAGEOUT_NAME ' + str(output) + \
                ' -WEIGHTOUT_NAME ' + re.sub('.fits', '', output) + '.weight.fits' + \
                ' -RESAMPLE_DIR ./ -RESAMPLE_SUFFIX .swarptemp.fits -RESAMPLING_TYPE ' +sampling+ \
                ' -COMBINE Y  -COMBINE_TYPE ' + str(combine_type)+ \
                ' -VERBOSE_TYPE NORMAL -SUBTRACT_BACK Y -INTERPOLATE Y' + \
                ' -PIXELSCALE_TYPE MANUAL,MANUAL -PIXEL_SCALE ' + str(pixelscale) + ',' + str(pixelscale) + \
                ' -IMAGE_SIZE ' + str(_imagesize) + ',' + str(_imagesize) + \
                ' -CENTER_TYPE MANUAL,MANUAL -CENTER ' + str(_ra) + ',' + str(_dec) + \
                ' -GAIN_DEFAULT ' + str(_gain)

        if imgmask:
           line += ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ' +  ','.join(imgmask)
        print line
        os.system(line)

    # the output of swarp give the normalized weight 
    # we want to store the variance image
    # we invert the normalized weight and we "de-normalized" this image
    if survey == 'sloan' or survey =='ps1':
        hd = fits.getheader(output)
        ar = fits.getdata(output)
        hd2 = fits.getheader(re.sub('.fits', '', output) + '.weight.fits')
        ar2 = fits.getdata(re.sub('.fits', '', output) + '.weight.fits')
    else:
        print(imglist, imgmask, 'imagelist abd imgmask')
        hd = fits.getheader(imglist[3])
        ar = fits.getdata(imglist[3])
        hd2 = fits.getheader(imgmask[0])
        ar2 = fits.getdata(imgmask[0])
    print(ar, 'odaayy')
    #hd2 = fits.getheader(re.sub('.fits', '', output) + '.weight.fits')
    #ar2 = fits.getdata(re.sub('.fits', '', output) + '.weight.fits')
    variance = 1 / ar2
    #  this is to take in account that the weight is normalized
    variance *= (np.median(np.abs(ar - np.median(ar)))*1.48)**2/np.median(variance)
    varimg = re.sub('.fits', '', output) + '.var.fits'
    fits.writeto(varimg, variance, hd2, overwrite=True)

    # put the saturation all values where the weight is zero 
    ar = np.where(ar2 == 0, _saturate, ar)
    print(ar2, 'ar2')

    if len(skylevel):
        hd['SKYLEVEL'] = (np.mean(skylevel), 'average sky level')
    hd['L1FWHM']   = (9999,       'FHWM (arcsec) - computed with sextractor')
    hd['WCSERR']   = (0,          'error status of WCS fit. 0 for no error')
    #hd['MJD-OBS']  = (_mjd,       'MJD')
    hd['RA']       = (_ra,        'RA')
    hd['DEC']      = (_dec,       'DEC')
    hd['RDNOISE']  = (_ron,       'read out noise')
    hd['PIXSCALE'] = (pixelscale, '[arcsec/pixel] Nominal pixel scale on sky')
    hd['FILTER']   = (_filter,    'filter used')
    hd['DAY-OBS']   = (_dayobs,    'day of observation')
    #hd['AIRMASS']  = (_airmass,   'airmass')
    #hd['DATE-OBS'] = (_dateobs,   'date of observation')
    hd['GAIN']     = (_gain,      'gain')
    hd['SATURATE'] = (_saturate,  'saturation level ')
    if objname:
        hd['OBJECT'] = (objname, 'title of the data set')
    if survey == 'sloan':
        hd['TELESCOP'] = ('SDSS', 'name of the telescope')
        hd['INSTRUME'] = ('SDSS', 'instrument used')
        hd['SITEID'] = ('SDSS', 'ID code of the Observatory site')
    elif survey == 'ps1':
        hd['TELESCOP'] = ('PS1', 'name of the telescope')
        hd['INSTRUME'] = ('PS1', 'instrument used')
        hd['SITEID'] = ('PS1', 'ID code of the Observatory site')
    elif survey == 'decam':
        hd['TELESCOP'] = ('decam', 'name of the telescope')
        hd['INSTRUME'] = ('decam', 'instrument used')
        hd['SITEID'] = ('decam', 'ID code of the Observatory site')
        hd['EXPTIME'] = (_exp, 'exposure time')

    ar, hdr = northupeastleft(data=ar, header=hd)
    out_fits = fits.PrimaryHDU(header=hd, data=ar)
    out_fits.writeto(output, overwrite=True, output_verify='fix')
    northupeastleft(filename=varimg)
    if show:
       lsc.display_image(output,2,True,'','')
    print('outout and varimg', output,varimg)
    return output, varimg

def northupeastleft(filename='', data=None, header=None):
    if filename:
        data, header = fits.getdata(filename, header=True)
    else:
        header = header.copy()
    if abs(header['cd1_2']) > abs(header['cd1_1']):
        data = data.T
        header['cd1_1'], header['cd1_2'] = header['cd1_2'], header['cd1_1']
        header['cd2_2'], header['cd2_1'] = header['cd2_1'], header['cd2_2']
        header['crpix1'], header['crpix2'] = header['crpix2'], header['crpix1']
        header['naxis1'], header['naxis2'] = header['naxis2'], header['naxis1']
        header['datasec'] = '[' + ','.join(header['datasec'].strip('[]').split(',')[::-1]) + ']'
        print 'swapping x and y axes'
    if header['cd1_1'] > 0:
        data = data[:,::-1]
        header['cd1_1'] *= -1
        header['cd2_1'] *= -1
        print 'flipping around x'
    if header['cd2_2'] < 0:
        data = data[::-1]
        header['cd2_2'] *= -1
        header['cd1_2'] *= -1
        print 'flipping around y'
    if filename:
        fits.writeto(filename, data, header, overwrite=True, output_verify='fix')
    else:
        return data, header

##################################################################################
def sloanimage(img,survey='sloan',frames=[], show=False, force=False):
   import sys
   from lsc import readhdr, readkey3,deg2HMS,display_image
   if show:
      display_image(img,1,True,'','')
   hdr = readhdr(img)
   _ra = readkey3(hdr,'RA')
   _dec = readkey3(hdr,'DEC')
   _object = readkey3(hdr,'object')
   _instrume = readkey3(hdr,'instrume')
   _filter = readkey3(hdr,'filter')
   _radius = 1100
   filt={'up':'u','gp':'g','rp':'r','ip':'i','zs':'z'}

   if _filter in filt.keys():
      _band = filt[_filter]
   else:
      _band = _filter

   if 'fs' in _instrume:
      _telescope = 'spectral'
   elif 'fl' in _instrume:
      _telescope = 'sinistro'
   elif 'fa' in _instrume:
      _telescope = 'sinistro'
   elif 'kb' in _instrume:
      _telescope = 'sbig'
   elif 'ep' in _instrume:
       _telescope = 'muscat'

   print _ra, _dec, _band, _radius
   if survey == 'sloan':
      frames =  downloadsdss(_ra, _dec, _band, _radius, force)
   elif survey =='decam':
      print('recognizes survey decam')
      _fov = 20 # setting the field of view to 20 for now
      frames = download_decam(_ra, _dec, _fov, _band)
   elif survey == 'ps1':
      if len(frames) == 0:
         delta= 0.22
         DR = delta/np.cos(_dec*np.pi/180)
         DD = delta
         f=open(_object+'_'+_band+'_ps1request.txt','w')
         f.write('%s   %s   %s\n' %(str(_ra),str(_dec),str(_band)))
         f.write('%s   %s   %s\n' %(str(_ra+DR),str(_dec+DD),str(_band)))
         f.write('%s   %s   %s\n' %(str(_ra-DR),str(_dec-DD),str(_band)))
         f.write('%s   %s   %s\n' %(str(_ra+DR),str(_dec-DD),str(_band)))
         f.write('%s   %s   %s\n' %(str(_ra-DR),str(_dec+DD),str(_band)))
         f.close()
         print '#'*20
         print "Please submit file:"+_object+'_'+_band+'_ps1request  at this link (you need an account)\n'
         print "   http://psps.ifa.hawaii.edu/PSI/postage_stamp.php    "
         print "   select:\n       Survey ID:    3PI          3PI.PV3  \n  "
         print "   Image Type:       Total Stacked Image, 1 pixel 0.250 "
         print "   Method of Image Selection by:   Coordinate- Images are selected based on supplied Ra and DEC   "
         print "   Size of the Postage Stamp:   (x) Arc-Second             Width:  2000         Height: 2000    pixel"
         print "   Center Coordinate of Image  - Input coordinate From: "
         print "                  (x)  Upload File           "
         print "        Choose File:    "+_object+'_'+_band+'_ps1request.txt'
         print '#'*20
         answ = raw_input('Do you want to wait that the frames are downloded [y/n]  [y] ? ')
         if not answ: 
             answ= 'y'
         if answ in ['Yes','yes','Y','y']:
              print 'Which is the last req_name at this page: '
              answ1 = raw_input('http://psps.ifa.hawaii.edu/PSI/postage_stamp_results.php ? ' )
              if not answ1:  
                  sys.exit('no name provided, no PS1 images have been downloaded')
              else:
                  print 'try download........'
                  frames = downloadPS1('./',answ1)
      else:
          frames2 = []
          for img in frames:
               if '_'+_band+'_' in img:
                    frames2.append(img)

          frames = frames2
   if len(frames):
        if survey == 'sloan' or survey == 'ps1' :
             print('input for sdss_swarp',frames,_telescope,_ra,_dec,'',_object, survey)
             out, varimg = sdss_swarp(frames,_telescope,_ra,_dec,'',_object, survey, show=show)
        else:
             out, varimg = name_change(frames, _object, _telescope)
   else:
       sys.exit('exit, no PS1 images have been downloaded')
   print(out, varimg, 'out and varimg odaayyy')
   return out, varimg
####################################################


def downloadPS1(homedir,filename):
    import sys
    import string
    import os
    import re
    import urllib
    import urllib2
    parent_dir = os.getcwd()+'/'
    directory=filename+'/'
    os.system('rm '+str(homedir)+'pssindex.txt')
    dataStoreBase = 'http://datastore.ipp.ifa.hawaii.edu/pstampresults/'
    dataStoreProduct = directory
    urlOfFileToDownload = dataStoreBase + dataStoreProduct + 'index.txt'
    print urlOfFileToDownload
    localFilename = homedir+'pss' + os.path.basename(urlOfFileToDownload)
    try:
       urllib.urlretrieve(urlOfFileToDownload, localFilename)
    except:
        print 'req_name found on poststamp datastore'
        sys.exit()
    frames = []
    try:
        f=file(homedir+'pssindex.txt','r')
        _index=f.readlines()
        f.close()
        for i in range(0,len(_index)):
            print string.split(_index[i],'|')[0],string.count(string.split(_index[i],'|')[0],'results.fits')
            if  string.count(string.split(_index[i],'|')[0],'.fits') and string.count(string.split(_index[i],'|')[0],'results.fits')==0:

                urlOfFileToDownload3 = dataStoreBase + dataStoreProduct +string.split(_index[i],'|')[0]
                urlOfFileToDownload3=re.sub('.txt','',urlOfFileToDownload3)
                localFilename3 = parent_dir+ os.path.basename(urlOfFileToDownload3)
                if not os.path.isfile(localFilename3):
                    print 'downloading file: '+string.split(_index[i],'|')[0]
                    try:
                        urllib.urlretrieve(urlOfFileToDownload3, localFilename3)
                    except:
                        print 'problem '+str(urlOfFileToDownload3)
                        pass
                else:
                    print 'file already downloaded '

                if not os.path.isdir(homedir+filename):
                    os.mkdir(homedir+filename)
                os.system('mv '+string.split(localFilename3,'/')[-1]+' '+homedir+filename+'/')
                if 'unconv.fits' in string.split(localFilename3,'/')[-1]:
                           frames.append(homedir+filename++'/'+string.split(localFilename3,'/')[-1])
    except:
        print 'stamp_directory not found '
        sys.exit()
    os.system('rm '+str(homedir)+'pssindex.txt')
    print 'download COMPLETE '
    return frames

################################################################################################################

def download_decam(_ra, _dec, _fov, _band):

        #from astroquery.noirlab import Noirlab
        import astropy.coordinates as coord
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        from astropy.visualization import ZScaleInterval

        # std lib
        from getpass import getpass
        import warnings
        from astropy.utils.exceptions import AstropyWarning
        warnings.simplefilter('ignore', category=AstropyWarning) # to quiet Astropy warnings

        # 3rd party
        import numpy as np
        from numpy.core.defchararray import startswith
        import pylab as plt
        import matplotlib

        from pyvo.dal import sia
        from astropy.utils.data import download_file
        from astropy.io import fits
        from astropy.wcs import WCS
        from astropy.visualization import make_lupton_rgb
        import urllib
        from urllib import urlretrieve
        
        frames=[]
        #DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/ls_dr7"
        #svc_ls_dr7 = sia.SIAService(DEF_ACCESS_URL)
        
        DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/ls_dr9"
        svc_ls_dr9 = sia.SIAService(DEF_ACCESS_URL)


        imgTable = svc_ls_dr9.search((_ra,_dec), (_fov/np.cos(_dec*np.pi/180), _fov), verbosity=2).to_table()
        if len(imgTable)==0:
                print('No data with these coordinates available')

        base_path= os.getcwd()
                            
        sel = (imgTable['proctype'] == 'Stack') & (imgTable['instrument_name'] == 'DECam') & \
                            (startswith(imgTable['obs_bandpass'].astype(str), _band))

        row = imgTable[sel]
        for i in range(len(row)):
                row = imgTable[sel][i]
                name_image= row['obs_publisher_did'].split('/')[-1][:-3]
                output = base_path+'/'+name_image
                url = row['access_url'] # get the download URL
                #output = urllib.request.urlretrieve(url, output)
                #output =urllib.urlopen(url, output)
                try:
                    requests = urllib.urlretrieve(url, output)
                except:
                    print 'problem '+str(url)
                    continue
                
                frames.append(output)
        print(frames, 'este')
        
        '''
        frames = ['/Users/estefaniapadilla/pipeline_data/data_reduction/22joj/legacysurvey-2203p030-image-g.fits',
                '/Users/estefaniapadilla/pipeline_data/data_reduction/22joj/legacysurvey-2203p030-invvar-g.fits',
                 '/Users/estefaniapadilla/pipeline_data/data_reduction/22joj/legacysurvey-2203p030-model-g.fits' ,
                 '/Users/estefaniapadilla/pipeline_data/data_reduction/22joj/legacysurvey-2203p030-galdepth-g.fits',
                  '/Users/estefaniapadilla/pipeline_data/data_reduction/22joj/legacysurvey-2203p030-nexp-g.fits',
                '/Users/estefaniapadilla/pipeline_data/data_reduction/22joj/legacysurvey-2203p030-depth-g.fits',
                '/Users/estefaniapadilla/pipeline_data/data_reduction/22joj/legacysurvey-2203p030-chi2-g.fits']
        '''
        return frames


######################################################################################################################
# we are going to change the name of the downloaded files and output a var image and an image with the desired keywords


def name_change(imglist, objname ='', _telescope=''):
    import string
    import re
    import os
    
    ##TO DO, instead of dropping blobmodel, check header instead and do try except
    #hdr = fits.getheader(imglist[0])
    imglist = [i for i in imglist if not 'blobmodel' in i]
    


    output=''
    out1 = 'decam'
    if _telescope == 'spectral':
        pixelscale = 0.30104  # 2 meter
        _imagesize =  2020
    elif _telescope == 'sbig':
            pixelscale = 0.467  # 1 meter
            _imagesize =  2030
    elif _telescope == 'sinistro':
            pixelscale = 0.387  # 1 meter
            _imagesize =  4020

    
    #create a copy with _1
    imglist2=[]
    for i,img in enumerate(imglist):
        
        hdr = fits.open(img)
        if os.path.isfile(re.sub('.fits','_1.fits',img)):
            os.system('rm '+re.sub('.fits','_1.fits',img))
        fits.writeto(re.sub('.fits','_1.fits',img),hdr[0].data,hdr[0].header)
        imglist2.append(re.sub('.fits','_1.fits',img))
        
    imglist = imglist2
    
    #get the wtmap
    imgmask = []
    skylevel = []
    wtlist = [i for i in imglist if '-invvar' in i]
    mklist = [i for i in imglist if '-maskbits' in i]
    explist = [i for i in imglist if '-nexp' in i]

    if len(wtlist):
        for i,name in enumerate(wtlist):
            #print(wtlist, 'wtlist')
            weightimg = re.sub('-invvar','.weight.',name)
            #print(weightimg, 'wetting')
            #hdr.writeto(weightimg, overwrite=True)
            imgmask.append(weightimg)
            if re.sub('-invvar','-maskbits',name) in mklist:
                print('in mklist', mklist)
                weight_data, weight_header = fits.getdata(name, header=True)
                mask_data, mask_header = fits.getdata(re.sub('-invvar','-maskbits',name), header=True)
                weight_data = 1/weight_data
                weight_fits = fits.PrimaryHDU(header=weight_header, data=weight_data)
                weight_fits.writeto(weightimg, output_verify='fix', overwrite=True)
            else:
                os.system('cp '+name+' '+weightimg)
    
    print('made it through weights, new update')
    #get the image list
    imglist = [i for i in imglist if 'image' in i]
    

    # Open the FITS file of image
    print(imglist, imgmask)
    image_file = fits.open(imglist[0])
    wtmap_file = fits.open(imgmask[0])
    exp_file = fits.open(explist[0])
    
    # Read the data
    data_image = image_file[0].data#*exp_file[0].data
    data_weights = wtmap_file[0].data
    
    
    
    #get header key word of imglist
    hdr = fits.getheader(imglist[0])
    _filter = hdr.get('FILTER').split(' ')[0]
    #_filter = string.split(hdr.get('FILTER'),' ')[0]
    #_dayobs = string.split(hdr.get('DATEOBS2'),'T')[0]# DATE-OBS for the last  image in the stack
    _dayobs = hdr.get('DATEOBS2').split('T')[0]
    _dayobs = _dayobs.replace('-', '')
    _saturate = 61000 # used for sloan
    _gain = hdr.get('ARAWGAIN')
    _ron = hdr.get('RDNOISEA')
    _ra = hdr.get('RA')
    _dec = hdr.get('DEC')
    
    if not _ra:
        _ra = hdr.get('CRVAL1')
    if not _dec:
        _dec = hdr.get('CRVAL2')
    
    if 'airmass' in hdr:
        _airmass = hdr.get('airmass')
    else:
        _airmass = 1

    #  filters
    filt={'U':'U','B':'B','V':'V','R':'R','I':'I','u':'up','g':'gp','r':'rp','i':'ip','z':'zs'}
    if _filter in filt.keys():
        _filter = filt[_filter]
        
    #update the header key words
    if len(skylevel):
        hdr['SKYLEVEL'] = (np.mean(skylevel), 'average sky level')
    hdr['L1FWHM']   = (9999,       'FHWM (arcsec) - computed with sextractor')
    hdr['WCSERR']   = (0,          'error status of WCS fit. 0 for no error')
    #hd['MJD-OBS']  = (_mjd,       'MJD')
    hdr['RA']       = (_ra,        'RA')
    hdr['DEC']      = (_dec,       'DEC')
    hdr['RDNOISE']  = (_ron,       'read out noise')
    hdr['PIXSCALE'] = (pixelscale, '[arcsec/pixel] Nominal pixel scale on sky')
    hdr['FILTER']   = (_filter,    'filter used')
    hdr['DAY-OBS']   = (_dayobs,    'day of observation')
    #hd['AIRMASS']  = (_airmass,   'airmass')
    #hd['DATE-OBS'] = (_dateobs,   'date of observation')
    hdr['GAIN']     = (_gain,      'gain')
    hdr['SATURATE'] = (_saturate,  'saturation level ')
    if objname:
        hdr['OBJECT'] = (objname, 'title of the data set')
    hdr['TELESCOP'] = ('decam', 'name of the telescope')
    hdr['INSTRUME'] = ('decam', 'instrument used')
    hdr['SITEID'] = ('decam', 'ID code of the Observatory site')

    # Close the FITS file
    image_file.close()

    # Rename the file
    if not output:
        output = _telescope+'_'+str(out1)+'_'+str(_dayobs)+'_'+str(_filter)+'_'+objname
        
    #create new file
    #os.rename(image_file, output+'.fits')
    
    #create variance file
    output_img = output+'.fits'
    varimg = output+'var.fits'
    image_fits = fits.PrimaryHDU(header=hdr, data = data_image)
    image_fits.writeto(output_img, overwrite=True, output_verify='fix')
    var_fits = fits.PrimaryHDU(header=hdr, data= 1/data_weights)
    var_fits.writeto(varimg, overwrite=True, output_verify='fix')
    print(image_fits, var_fits, 'image_fits and var_fits')
    
    return output_img, varimg













