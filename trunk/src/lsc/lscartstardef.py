import re
import sys
import string
import os
import numpy as np
import lsc
import tempfile



    
def artstar(img,_ra,_dec,_num=1, _app = 17 , _random=False,_verbose=True):
    from astropy.io import fits as pyfits
    from astropy import wcs as pywcs
    from pyraf import iraf
    from iraf import digiphot
    from iraf import daophot

#    try:
    hostname, username, passwd, database=lsc.mysqldef.getconnection('lcogt2')
    conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
#    except:
#        print('problem with the connection to the database')
    print(conn)

    hdr = pyfits.open(img)[0].header
    ggg = lsc.mysqldef.getfromdataraw(conn, 'photlco', 'filename', os.path.basename(img), '*')
    
    if 'exptime' in hdr:
        _exptime = hdr['exptime']
    elif ggg:
        if 'exptime' in ggg[0]:
            _exptime = ggg[0]['exptime']
        else:
            _exptime = None
    else:
            _exptime = None
            
    if  'ZP' in hdr:
        _zeropoint = hdr['ZP']
    elif ggg:
        print(ggg)
        if 'zn' in ggg[0]:
            _zeropoint = ggg[0]['zn']
        else:
            _zeropoint = None
    else:
            _zeropoint = None
        
    if _zeropoint is None:
        sys.exit('zeropoint not found')
    if _exptime is None:
        sys.exit('exptime not found')
        
    if _zeropoint is not None:
        _mag =  _app  - _zeropoint
        print _app, _mag, _zeropoint
    else:
        sys.exit('zeropoint or exptime header not found')

    psffile = re.sub('.fits','.psf.fits',img)
    if os.path.isfile(psffile):
        print('psf file found')
        _ra0 = lsc.util.readkey3(hdr,'RA')
        _dec0 = lsc.util.readkey3(hdr,'DEC')
        xx = ''
        yy = ''

        if _ra and _dec:
            print _ra0, _dec0, _ra,_dec
            
            w = pywcs.WCS(hdr)
            apass = zip([_ra],[_dec])
            _ra = [_ra]
            _dec = [_dec]
            pixel_apass = w.wcs_world2pix(apass, 1)
            xx = [pixel_apass[0][0]]
            yy = [pixel_apass[0][1]]
            #xx, yy = zip(*pixel_apass)
            print(xx,yy)
            if xx[0] < 0 or xx[0] > hdr['NAXIS1']:
                print('WRONG coordinate')
                sys.exit()
            if yy[0] < 0 or yy[0] > hdr['NAXIS2']:
                print('WRONG coordinate')
                sys.exit()

        else:
#            print '_ra and _dec not defined' 
#            print 'chose a random position'
#            xx = 50+(hdr['NAXIS1']-100)*np.random.random(int(_num))
#            yy = 50+(hdr['NAXIS2']-100)*np.random.random(int(_num))
#############################################
            import math
            print '_ra and _dec not defined.' 
            print 'Choosing random positions.'
            print 'EDITED VERSION 2.3'
            #Gathering all of the randomized valid positions
            ###################################################
            #Masked image area has values of 1
            #Setup art star lists
            coord_list = []
            xx = np.array([50 + (hdr['NAXIS1']-100) * np.random.random(int(1))[0]])
            yy = np.array([50 + (hdr['NAXIS2']-100) * np.random.random(int(1))[0]])

            #Find the mask data if file exists
            imgmask = re.sub('.fits', '.mask.fits', img)
            masked = False
            if os.path.isfile(imgmask):
                masked = True
                hdu = pyfits.open(imgmask)
                data = hdu[0].data
            print("Masked file:", masked)
            #Distance threshold: How far radially all artstars must be away from all other artstars at minimum
            threshold = 20
            loopcount = 0
            while len(xx) < _num:
                loopcount += 1
                #Grab new position
                tempx = 50+(hdr['NAXIS1']-100)*np.random.random(int(1))[0]
                tempy = 50+(hdr['NAXIS2']-100)*np.random.random(int(1))[0]
                #Find if distance to all other artstars and check versus threshold distance
                tempxdif = xx - tempx
                tempydif = yy - tempy
                #Distance is list of T/F. False means it's safe distance
                distance = np.sqrt(tempxdif**2+tempydif**2)
                dlt_thresh = distance  <= threshold
                
                #Failsafe against loop going forever with too many artstars
                #If it loops more than threshold*num total stars then just add the star
                #10+(logx)^4 gives values in low 10s for low star count
                #up to around 100-150 for num stars around 1000-2000
                if (loopcount <= (10+int(pow(math.log(_num,10),4)))):
                    #Check if on the mask
                    if masked:
                        if data[int(tempx),int(tempy)]:
                            #print("Position on the mask. Repositioning.")
                            continue
                        #Check if too close
                    if np.any(dlt_thresh):
                        #print('Too close to another star. Repositioning')
                        continue
                else:
                    print("Loop approaching infinity. Adding random star.")
                #Star should be good, add to final list
                print(tempx,tempy)
                print(data[int(tempx),int(tempy)])
                xx = np.append(xx,tempx)
                yy = np.append(yy,tempy)
                loopcount = 0
                ############################################






############################################            
            w = pywcs.WCS(hdr)
            apass = zip(xx,yy)
            coor_apass = w.wcs_pix2world(apass, 1)
            _ra,_dec = zip(*coor_apass)
            

        if len(xx) and len(yy):
          print(xx)
          print(yy)
          if float(_mag) != 0.0:
              if _random:
                  #
                  # take a random distribution of mag arounf _mag 
                  _mag = ((3 * np.random.random(_num))-3./2 )+_mag
              else:
                  _mag = np.zeros(len(xx)) + _mag
                  
              ddd = ''
              for jj,kk in enumerate(xx): 
                  print jj,kk
                  ddd = ddd + '%s  %s  %s\n' % (xx[jj],yy[jj],_mag[jj])
                  #ddd = ddd + '%s  %s  %s\n' % (xx[jj],yy[jj],_mag)

              addobj = next(tempfile._get_candidate_names())
              tmp1 = next(tempfile._get_candidate_names())+'.fits'
              tmp2 = next(tempfile._get_candidate_names())+'.fits'

              ff=open(addobj,'w')
              ff.write(ddd)
              ff.close()

              print 'add star'
              imgout=re.sub('.fits','.art.fits',img)
              print imgout
              os.system('rm -rf '+ tmp1)
              os.system('rm -rf '+ tmp2)
              os.system('rm -rf '+ imgout)
              iraf.imarith(img, '-', img, tmp1, verbose='yes')
              iraf.daophot.addstar(tmp1, addobj, psffile, tmp2, nstar=1, veri='no', simple='yes',
                                   verb='yes')
              iraf.imarith(img, '+', tmp2, imgout, verbose='yes')

              lsc.util.delete(addobj)                        
              lsc.util.delete(tmp1) 
              lsc.util.delete(tmp2)
              if os.path.isfile(tmp1+'.art'):
                  lsc.util.delete(tmp1+'.art')
              if os.path.isfile(tmp2+'.art'):
                  lsc.util.delete(tmp2+'.art')

              ################# insert new image in the archive
              _object = hdr.get('object')
              _gain = hdr.get('gain')
              _rdnoise = hdr.get('rdnoise')
              _saturate = hdr.get('saturate')
              _wcs = hdr.get('wcs')
              _targetid = lsc.mysqldef.targimg(imgout)

              if ggg:
                  _dayobs = ggg[0]['dayobs']
                  telescopeid = ggg[0]['telescopeid']
                  instrumentid = ggg[0]['instrumentid']
                  _wcs = ggg[0]['wcs']
              else:
                  _dayobs = None
                  telescopeid = None
                  instrumentid = None 
                  _wcs = None
                  
              
              dictionary={'telescope':lsc.util.readkey3(hdr,'telescop'),'instrument':lsc.util.readkey3(hdr,'instrume'),'dec0':lsc.util.readkey3(hdr,'DEC'),\
                            'ra0': lsc.util.readkey3(hdr,'RA'),'ut': lsc.util.readkey3(hdr,'ut'), 'dateobs': lsc.util.readkey3(hdr,'date-obs'),\
                            'exptime': lsc.util.readkey3(hdr,'exptime'), 'filter': lsc.util.readkey3(hdr,'filter'),'mjd': lsc.util.readkey3(hdr,'MJD-OBS'),\
                            'airmass': lsc.util.readkey3(hdr,'airmass'), 'objname':_object,'targetid': _targetid,'filetype' : 5,\
                            'dayobs': _dayobs,'wcs':_wcs, 'telescopeid':telescopeid, 'instrumentid':instrumentid,\
                            #'obsid': lsc.util.readkey3(hdr,'OBSID'),\
                            'wcs' : _wcs}
              dictionary['filename']=string.split(imgout,'/')[-1]
              dictionary['filepath']=re.sub(dictionary['filename'],'',imgout)
              
              if not lsc.mysqldef.getfromdataraw(conn,'photlco','filename', string.split(imgout,'/')[-1],column2='filename'):
                  lsc.mysqldef.insert_values(conn,'photlco',dictionary)
              else:
                  for voce in dictionary:
                      if voce!='id' and voce!='filename':
                          lsc.mysqldef.updatevalue('photlco', voce,dictionary[voce],string.split(imgout,'/')[-1])
#########################                          
#              #################    insert fake detection  ################################
              _file = os.path.basename(imgout)
              _filepath = re.sub(_file,'',imgout)
              ggg1 = lsc.mysqldef.getfromdataraw(conn, 'fakecandidates', 'filename', os.path.basename(imgout), 'filename')        
              if ggg1:
                  lsc.mysqldef.deleteredufromarchive(os.path.basename(imgout),'fakecandidates','filename')
                  print('cancel old fake candidates from this frames')

              if _verbose:
                  X, hdr = pyfits.getdata(imgout, header=True)
                  _z1,_z2 = lsc.zscale(X)
                  from matplotlib import pylab as plt
                  plt.ion()
                  fig = plt.figure()
                  ax = fig.add_axes([0.1,0.1,0.85,0.85])
                  image = ax.imshow(X,interpolation='nearest', origin='upper', cmap='gray_r', vmin=_z1, vmax=_z2)
                  
              for jj,kk in enumerate(xx):
                  if _verbose:
                      ax.plot(xx[jj],yy[jj],'^', markeredgecolor = 'r', mfc='r',markersize=10)
                      fig.canvas.draw()
                      fig.canvas.flush_events()
                  
                  dictionary={'ra0' : _ra[jj], 'dec0': _dec[jj], 'xpos': xx[jj] ,'ypos': yy[jj], 'filename': _file,
                              'filepath': _filepath, 'magauto': _mag[jj], 'mjd': hdr['MJD-OBS']}
                  lsc.mysqldef.insert_values(conn,'fakecandidates', dictionary)
##########################
               #
              #_app = _mag  (np.log10(_exptime) * 2.5 + _zeropoint) 
              return _ra,_dec,xx,yy,_mag,_exptime,_zeropoint
          else:
              return '','','','','','',''
        else:
            return '','','','','','',''
    else:
        print ('psf file not found')
        return '','','','','','',''
