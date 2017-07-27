import sys
import os
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import lsc

workdirectory = os.getenv('LCOSNDIR')
if workdirectory is None:
    workdirectory = '/supernova/'
configfile = os.path.join(workdirectory, 'configure')

def readpasswd(configfile):
   """ read all information to connect to database from configuration file  
   """
   from numpy import genfromtxt
   data = genfromtxt(configfile, str)
   gg= {}
   for i in data:
      try:
         gg[i[0]] = eval(i[1])
      except:
         gg[i[0]] = i[1]
   return gg

#############################################################################################
#   read information from configuration file
#######################################################################################
readpass = readpasswd(configfile)
proposal = readpass['proposal']
users = readpass['users']
triggerpass = readpass['triggerpass']
extraobject = readpass['extraobject']
skipobjects = readpass['skipobjects']
proposalingestion = readpass['proposalingestion']
#############################################################################################

def ReadAscii2(ascifile):
   import string
   f=open(ascifile,'r')
   ss=f.readlines()
   f.close()
   vec1,vec2=[],[]
   for line in ss:
      if line[0]!='#':
         vec1.append(float(string.split(line)[0]))
         vec2.append(float(string.split(line)[1]))
   return vec1,vec2
#########################################################################
def readlist(listfile):
    import string,os,sys,re,glob
    if '*' in listfile:
        imglist=glob.glob(listfile)
    elif ',' in listfile: imglist = string.split(listfile,sep=',')
    else:
        try:            hdulist= fits.open(listfile)
        except:           hdulist=[]
        if hdulist:            imglist = [listfile]
        else:
           try:
              ff = open(listfile,'r')
              files = ff.readlines()
              ff.close()
              imglist = []
              for ff in files: 
                 ff=re.sub(' ','',ff)
                 if not ff=='\n' and ff[0]!='#':
                    ff=re.sub('\n','',ff)
                    try:
                       hdulist= fits.open(ff)
                       imglist.append(ff)
                    except Exception as e:
                        print 'problem reading header of', ff
                        print e
           except:              sys.exit('\n##### Error ###\n file '+str(listfile)+' do not  exist\n')
    if len(imglist)==0:
           sys.exit('\n##### Error ###\nIf "'+str(listfile)\
                                +'" is an image, it is corrupted \n or is not a list of image\n')
    return imglist
##############################################################################
def delete(listfile):
    import os,string,re,glob
    if listfile[0]=='@':   
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files: 
            ff=re.sub(' ','',ff)
            if not ff=='\n' and ff[0]!='#':
                ff=re.sub('\n','',ff)
                imglist.append(ff)
    elif ',' in listfile: imglist = string.split(listfile,sep=',')
    else:       imglist=[listfile]    
    lista=[]
    for _file in imglist:   lista=lista+glob.glob(_file)
    if lista:
        for _file in lista:
            try:          os.system('rm '+_file)
            except:       pass

def imcopy(imgin, imgout, center=None, cutout_size=None, ext=0):
    '''iraf.imcopy equivalent
       center = (x, y) of center of cutout
       cutout_size = (dy, dx) or single value for square
       By default, the entire 0th extension is copied.'''
    hdulist = fits.open(imgin)
    data = hdulist[ext].data
    hdr = hdulist[ext].header
    if center is not None and cutout_size is not None:
        cutout = Cutout2D(data, center, cutout_size, WCS(hdr))
        data = cutout.data
        hdr['CRPIX1'] = cutout.wcs.wcs.crpix[0]
        hdr['CRPIX2'] = cutout.wcs.wcs.crpix[1]
    fits.writeto(imgout, data, hdr)
    hdulist.close()

###############################################################
def readhdr(img):
    try:
        hdr = fits.getheader(img)
    except Exception as e:
        print "Couldn't read header of {}. Try deleting it and starting over.".format(img)
        raise e
    return hdr
    
missingvalues = ['NaN', 'UNKNOWN', None, '', 'UNSPECIFIED', 'N/A']

def readkey3(hdr,keyword):
    from astropy.coordinates import Angle
    from astropy import units as u

    try:    
       _instrume=hdr.get('INSTRUME').lower()
    except: 
       _instrume='none'
    if 'kb' in _instrume: # SBIG
        useful_keys = {'object'    : 'OBJECT',\
                           'date-obs'  : 'DATE-OBS',\
                           'ut'        : 'DATE-OBS',\
                           'date-night': 'DAY-OBS',\
                           'RA'        : 'RA',\
                           'DEC'       : 'DEC',\
                           'CAT-RA'    : 'CAT-RA',\
                           'CAT-DEC'   : 'CAT-DEC',\
                           'datamin'   :  -100.0,\
                           'datamax'   : 'SATURATE',\
                           'observer'  : 'OBSERVER',\
                           'exptime'   : 'EXPTIME',\
                           'wcserr'    : 'WCSERR',\
                           'instrume'  : 'INSTRUME',\
                           'JD'        : 'MJD-OBS',\
                           'mjd'        : 'MJD-OBS',\
                           'filter'    : 'FILTER',\
                           'gain'      : 'GAIN',\
                           'ron'       : 'RDNOISE',\
                           'airmass'   : 'AIRMASS',\
                           'type'      : 'OBSTYPE',\
                           'propid'      : 'PROPID',\
                           'userid'      : 'USERID',\
                           'telescop'  : 'TELESCOP'} 
    elif 'fl' in _instrume: # sinistro
        useful_keys = {'object'    : 'OBJECT',\
                           'date-obs'  : 'DATE-OBS',\
                           'ut'        : 'DATE-OBS',\
                           'date-night': 'DAY-OBS',\
                           'RA'        : 'RA',\
                           'DEC'       : 'DEC',\
                           'CAT-RA'    : 'CAT-RA',\
                           'CAT-DEC'   : 'CAT-DEC',\
                           'datamin'   :  -100.0,\
                           'datamax'   : 'SATURATE',\
                           'observer'  : 'OBSERVER',\
                           'exptime'   : 'EXPTIME',\
                           'wcserr'    : 'WCSERR',\
                           'instrume'  : 'INSTRUME',\
                           'JD'        : 'MJD-OBS',\
                           'mjd'       : 'MJD-OBS',\
                           'filter'    : 'FILTER',\
                           'gain'      : 'GAIN',\
                           'ron'       : 'RDNOISE',\
                           'airmass'   : 'AIRMASS',\
                           'type'      : 'OBSTYPE',\
                           'propid'    : 'PROPID',\
                           'userid'    : 'USERID',\
                           'telescop'  : 'TELESCOP'}
    elif 'fs' in _instrume or 'em' in _instrume:
       if hdr.get('DATE-OBS') < '2014-04-01':
          if 'RDNOISE' in hdr: ron='RDNOISE'
          elif 'READNOIS' in hdr: ron='READNOIS'
          else: ron='ron'
          useful_keys = {'object'    : 'OBJECT',\
                         'date-obs'  : 'DATE-OBS',\
                         'ut'        : 'DATE-OBS',\
                         'date-night': 'DAY-OBS',\
                         'RA'        : 'RA',\
                         'DEC'       : 'DEC',\
                         'CAT-RA'    : 'CAT-RA',\
                         'CAT-DEC'   : 'CAT-DEC',\
                         'datamin'   : -100.0,\
                         'datamax'   : 60000.0,\
                         'wcserr'    : 'WCS_ERR',\
                         'observer'  : 'OBSERVER',\
                         'exptime'   : 'EXPTIME',\
                         'instrume'  : 'INSTRUME',\
                         'JD'        : 'MJD',\
                         'mjd'       : 'MJD',\
                         'filter'    : 'FILTER',\
                         'gain'      : 'GAIN',\
                         'ron'       : ron,\
                         'airmass'   : 'AIRMASS',\
                         'userid'    : 'USERID',\
                         'propid'    : 'PROPID',\
                         'type'      : 'OBSTYPE',\
                         'telescop'  : 'TELID'} 
          if _instrume == 'fs02': # OGG
              useful_keys['pixscale'] = 0.30104
          elif _instrume in ['fs01', 'fs03']: # COJ
              useful_keys['pixscale'] = 0.304
       else:
          useful_keys = {'object'    : 'OBJECT',\
                            'date-obs'  : 'DATE-OBS',\
                            'ut'        : 'DATE-OBS',\
                            'date-night': 'DAY-OBS',\
                            'RA'        : 'RA',\
                            'DEC'       : 'DEC',\
                            'CAT-RA'    : 'CAT-RA',\
                            'CAT-DEC'   : 'CAT-DEC',\
                            'datamin'   :  -100.0,\
                            'datamax'   : 'SATURATE',\
                            'observer'  : 'OBSERVER',\
                            'exptime'   : 'EXPTIME',\
                            'wcserr'    : 'WCSERR',\
                            'instrume'  : 'INSTRUME',\
                            'JD'        : 'MJD-OBS',\
                            'mjd'        : 'MJD-OBS',\
                            'filter'    : 'FILTER',\
                            'gain'      : 'GAIN',\
                            'ron'       : 'RDNOISE',\
                            'airmass'   : 'AIRMASS',\
                           'type'      : 'OBSTYPE',\
                           'propid'      : 'PROPID',\
                           'userid'      : 'USERID',\
                           'telescop'  : 'TELESCOP'}
    else: 
       useful_keys = {'object'    : 'OBJECT',\
                      'RA'        : 'RA',\
                      'DEC'       : 'DEC',\
                      'CAT-RA'    : 'RA',\
                      'CAT-DEC'   : 'DEC',\
                      'ron'       : 'RDNOISE',\
                      'date-obs'  : 'DATE-OBS',\
                      'date-night': 'DAY-OBS',\
                      'datamax'   : 'SATURATE'}
    if keyword in useful_keys:
       if type(useful_keys[keyword])==float:
          value=useful_keys[keyword]
       else:
          value=hdr.get(useful_keys[keyword])
          if keyword=='date-obs':
             try:
                value = value.split('T')[0].replace('-', '')
             except:
                pass
          elif keyword=='ut':
             try:
                value = value.split('T')[1]
             except:
                pass
          elif keyword=='object':
             value = value.translate(None, ' }{][)(')
          elif keyword=='JD':       
             value=value+0.5
          elif keyword=='instrume':      value=value.lower()
          elif keyword=='filter' and value in [None, 'air']:
             for key in ['FILTER2', 'FILTER1', 'FILTER3']:
                if hdr.get(key) not in [None, 'air']:
                   value = hdr[key]
                   break
          elif keyword in ['RA', 'CAT-RA'] and type(value) == str and ':' in value:
             value = Angle(value, u.hourangle).deg
          elif keyword in ['RA', 'CAT-RA', 'DEC', 'CAT-DEC']:
             value = Angle(value, u.deg).deg
    elif keyword in hdr:
       value=hdr.get(keyword)
    else:
       value=''
    if type(value) == str:
       value = value.replace('\#', '')
    if value == 'ftn':
      value = '2m0-01'
    elif value == 'fts':
      value = '2m0-02'
    return value

#######################################################
def writeinthelog(text,logfile):
    f=open(logfile,'a')
    f.write(text)
    f.close()

def updateheader(filename, dimension, headerdict):
    tupledict = {key: tuple(value) for key, value in headerdict.items()}
    try:
        hdulist = fits.open(filename, mode='update')
        header = hdulist[dimension].header
        header.update(tupledict)
        hdulist.close()
    except Exception as e:
        print 'header of', filename, 'not updated:'
        print e
#################################################################################################
def display_image(img,frame,_z1,_z2,scale,_xcen=0.5,_ycen=0.5,_xsize=1,_ysize=1,_erase='yes'):
    goon='True'
    import glob, subprocess, os, time
    ds9 = subprocess.Popen("ps -U {:d} u | grep -v grep | grep ds9".format(os.getuid()),shell=True,stdout=subprocess.PIPE).stdout.readlines()
    if len(ds9)== 0 :   
       subproc = subprocess.Popen('ds9',shell=True)
       time.sleep(3)

    if glob.glob(img):
       from pyraf import iraf
       iraf.images(_doprint=0)
       iraf.tv(_doprint=0)
       import string,os
       if _z2: 
          try:
              sss=iraf.display(img + '[0]', frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,\
                                   fill='yes', zscale='no', zrange='no', z1=_z1, z2=_z2,Stdout=1)
          except:
              print ''
              print '### ERROR: PROBLEM OPENING DS9'
              print ''
              goon='False'                 
       else:
        try:  
            sss=iraf.display(img + '[0]', frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase, fill='yes', Stdout=1)
        except:
            print ''
            print '### ERROR: PROBLEM OPENING DS9'
            print ''
            goon=False
 
       if scale and goon:
          answ0 = raw_input('>>> Cuts OK ? [y/n] ? [y] ')
          if not answ0: answ0='y'
          elif answ0=='no' or answ0=='NO': answ0='n' 

          while answ0=='n':
              _z11=float(string.split(string.split(sss[0])[0],'=')[1])
              _z22=float(string.split(string.split(sss[0])[1],'=')[1])
              z11 = raw_input('>>> z1 = ? ['+str(_z11)+'] ? ')
              z22 = raw_input('>>> z2 = ? ['+str(_z22)+'] ? ')
              if not z11: z11=_z11
              else: z11=float(z11)
              if not z22: z22=_z22
              else: z22=float(z22)
              print z11,z22
              sss=iraf.display(img + '[0]',frame,fill='yes', xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,\
                                   zrange='no', zscale='no', z1=z11, z2=z22, Stdout=1)
              answ0 = raw_input('>>> Cuts OK ? [y/n] ? [y] ')
              if not answ0: answ0='y'
              elif answ0=='no' or answ0=='NO': answ0='n'
       if goon:
          _z1,_z2=string.split(string.split(sss[0])[0],'=')[1],string.split(string.split(sss[0])[1],'=')[1]
    else:
        print 'Warning: image '+str(img)+' not found in the directory '
    return _z1,_z2,goon

###########################################################################
def readstandard(standardfile):
    import lsc
    from numpy import array,abs
    import string,os
    if os.path.isfile(standardfile):
       listastandard=standardfile
    elif standardfile[0]=='/':
       listastandard=standardfile
    else:
       listastandard=lsc.__path__[0]+'/standard/stdlist/'+standardfile
    f=open(listastandard,'r')
    liststd=f.readlines()
    f.close()
    star,ra,dec=[],[],[]
    magnitude=[]
    for i in liststd:
       if i[0]!='#':
          star.append(string.split(i)[0])
          _ra=string.split(string.split(i)[1],':')
          _dec=string.split(string.split(i)[2],':')
          ra.append((float(_ra[0])+((float(_ra[1])+(float(_ra[2])/60.))/60.))*15)
          if '-' in str(_dec[0]):
             dec.append((-1)*(abs(float(_dec[0]))+((float(_dec[1])+(float(_dec[2])/60.))/60.)))
          else:
             dec.append(float(_dec[0])+((float(_dec[1])+(float(_dec[2])/60.))/60.))
          try:   magnitude.append(string.split(i)[3])
          except:  magnitude.append(999)
    return array(star),array(ra),array(dec),array(magnitude)

#################################################################################################
def readspectrum(img):
    from numpy import array
    import string
    fl=''
    lam=''
    graf=1
    spec = fits.open(img)
    head = spec[0].header
    try:
        if spec[0].data.ndim == 1: fl = spec[0].data
        elif spec[0].data.ndim == 2: fl = spec[0].data[:,0]
        elif spec[0].data.ndim == 3: fl = spec[0].data[0,0,:]
    except:
        if spec[0].data.rank == 1: fl = spec[0].data
        elif spec[0].data.rank == 2: fl = spec[0].data[:,0]
        elif spec[0].data.rank == 3: fl = spec[0].data[0,0,:]
    naxis1 = head['naxis1']
    try:
        crpix1 = head['crpix1']
        crval1 = head['crval1']
        try: cdelt1 = head['cdelt1']
        except: cdelt1 = head['cd1_1']
        pix = array(range(1,naxis1+1,1))
        pix = array(range(1,len(fl)+1,1))
        lam = (pix-crpix1)*cdelt1+crval1
    except:
        try:
           WAT= head['WAT2_001']
           pix = array(range(1,naxis1+1,1))
           crpix1=string.split(string.split(WAT,'"')[1])[0]
           crval1=string.split(string.split(WAT,'"')[1])[3]
           cdelt1=string.split(string.split(WAT,'"')[1])[4]
           lam = (pix-float(crpix1))*float(cdelt1)+float(crval1)
        except:
           graf=0
    return lam,fl
###########################################################################
def pval(_xx, p):
    _y=+p[0]+p[1]*_xx
    return _y

def residual(p,y,x):
    for i in range(len(p)):
        err = (y-p[i]*x**i)
    return err
#########################################################################
def defsex(filename):
    import lsc
    import string,re,os
    sexfile=lsc.__path__[0]+'/standard/sex/default.sex'
    f=open(sexfile,'r')
    ss=f.readlines()
    f.close()   
    ff=open(filename,'w')
    for i in ss:
        if string.count(i,'PARAMETERS_NAME')==1:
            ff.write('PARAMETERS_NAME  "'+lsc.__path__[0]+'/standard/sex/default.param"\n')
        elif string.count(i,'FILTER_NAME')==1:
            ff.write('FILTER_NAME  "'+lsc.__path__[0]+'/standard/sex/default.conv"\n')
        elif string.count(i,'STARNNW_NAME')==1:
            ff.write('STARNNW_NAME "'+lsc.__path__[0]+'/standard/sex/default.nnw"\n')
        else:
            ff.write(i)
    ff.close()
    return filename

############################################################

def defswarp(filename,imgname,_combine,gain='',ron='',pixelscale=0.4699,_ra='',_dec=''):
    import lsc
    import string,re,os
    if _combine.lower() in ['median']: _combine='MEDIAN'
    elif _combine.lower() in ['average']: _combine='AVERAGE'
    elif _combine.lower() in ['sum']: _combine='SUM'
    swarpfile=lsc.__path__[0]+'/standard/sex/default.swarp'
    f=open(swarpfile,'r')
    ss=f.readlines()
    f.close()   
    ff=open(filename,'w')
    for i in ss:
        if string.count(i,'IMAGEOUT_NAME')==1:
            ff.write('IMAGEOUT_NAME    '+str(imgname)+'  # Output filename \n')
        elif string.count(i,'WEIGHTOUT_NAME')==1:
            ff.write('WEIGHTOUT_NAME   '+str(re.sub('.fits','.weight.fits',imgname))+'  # Output weight-map filename  \n')
        elif string.count(i,'COMBINE_TYPE')==1:
            ff.write('COMBINE_TYPE    '+str(_combine)+'  # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CHI2 \n')
        elif string.count(i,'GAIN_DEFAULT')==1:
           if gain:
              ff.write('GAIN_DEFAULT    '+str(gain)+'  # Default gain if no FITS keyword found \n')
           else:     ff.write(i)
        elif string.count(i,'RDNOISE_DEFAULT')==1:
           if ron:
              ff.write('RDNOISE_DEFAULT    '+str(ron)+'  # Default ron if no FITS keyword found \n')
           else:     ff.write(i)
        elif string.count(i,'PIXEL_SCALE')==1:
              ff.write('PIXEL_SCALE    '+str(pixelscale)+','+str(pixelscale)+'  #  \n')
        elif string.count(i,'PIXELSCALE_TYPE')==1:
              ff.write('PIXELSCALE_TYPE MANUAL,MANUAL  #  \n')
        elif string.count(i,'CENTER_TYPE')==1:
              ff.write('CENTER_TYPE MANUAL,MANUAL  #  \n')
        elif string.count(i,'Coordinates of the image center')==1:
              ff.write('CENTER '+str(_ra)+','+str(_dec)+'  #  \n')
        else:
            ff.write(i)
    ff.close()
    return filename

#################################################################################
def airmass(img,overwrite=True,_observatory='lasilla'):
   import lsc
   from lsc.util import readhdr,readkey3, delete, updateheader
   from pyraf import iraf
   iraf.astutil(_doprint=0)
   hdr=readhdr(img)
   if readkey3(hdr,'UTC'):    
      _UT=(readkey3(hdr,'UTC')+(readkey3(hdr,'exptime')/2))/3600
      _date=readkey3(hdr,'date-obs')
      _date=_date[0:4]+'-'+_date[4:6]+'-'+_date[6:8]
      _RA=readkey3(hdr,'RA')/15
      _DEC=readkey3(hdr,'DEC')
      f = file('airmass.txt','w')
      f.write('mst = mst ("'+str(_date)+'",'+str(_UT)+', obsdb ("'+str(_observatory)+'", "longitude"))\n')
      f.write('air = airmass ('+str(_RA)+','+str(_DEC)+',mst, obsdb ("'+str(_observatory)+'", "latitude"))\n')
      f.write('print(air)\n')
      f.close()
      _air=iraf.astcalc(image=img, command="airmass.txt",Stdout=1)[0]
      try: _air=float(_air)
      except: _air=999
      delete('airmass.txt')
      if overwrite and _air<99.:
         updateheader(img,0,{'AIRMASS':[_air,'mean airmass computed with astcalc']})
   else:   _air=''
   return _air
####################################################################################

def  name_duplicate(img,nome,ext):  ###########################
   import re,string,os,glob
   import lsc
   from lsc.util import readhdr,readkey3, delete
   dimg=readkey3(readhdr(img),'DATE-OBS')
   listafile=glob.glob(nome+'_?'+ext+'.fits')+glob.glob(nome+'_??'+ext+'.fits')
   if len(listafile) == 0: nome = nome+"_1"+ext+'.fits'
   else:
      date=[]
      for l in listafile:
         date.append(readkey3(readhdr(l),'DATE-OBS'))
      if dimg in date:
         nome=listafile[date.index(dimg)]
#         if overwrite:
#            delete(nome)
      else:
         n=1
         while nome+'_'+str(n)+str(ext)+'.fits' in listafile:
            n=n+1
         nome=nome+'_'+str(n)+str(ext)+'.fits'
   return nome
###############################################################################
def correctobject(img,coordinatefile):
    import os,re,sys,string
    from numpy import arccos, sin,cos,pi,argmin
    import lsc
    from lsc.util import readstandard, readhdr,readkey3, updateheader
    scal=pi/180.
    std,rastd,decstd,magstd=readstandard(coordinatefile)
    img=re.sub('\n','',img)
    hdr=readhdr(img)
    _ra=readkey3(hdr,'RA')
    _dec=readkey3(hdr,'DEC')
    dd=arccos(sin(_dec*scal)*sin(decstd*scal)+cos(_dec*scal)*cos(decstd*scal)*cos((_ra-rastd)*scal))*((180/pi)*3600)
    if min(dd)<200:
       updateheader(img,0,{'OBJECT':[std[argmin(dd)],'Original target.']})
       aa,bb,cc=rastd[argmin(dd)],decstd[argmin(dd)],std[argmin(dd)]
    else: aa,bb,cc='','',''
    return aa,bb,cc

##################################################################################################
def repstringinfile(filein,fileout,string1,string2):
    import re
    f=open(filein,'r')
    ss=f.readlines()
    f.close()
    f=open(fileout,'w')
    for n in range(len(ss)):
        if string1 in ss[n]:   f.write(re.sub(string1,string2,ss[n]))
        else:                 f.write(ss[n])
    f.close()
###################################################

def limmag(img):
   import lsc
   from lsc.util import readhdr
   hdr=readhdr(img)
   _ZP=readkey3(hdr,'PHOTZP')
   _gain=readkey3(hdr,'gain')
   _exptime=readkey3(hdr,'exptime')
   _fwhm=readkey3(hdr,'PSF_FWHM')  
   _mbkg=readkey3(hdr,'MBKG')   # background from sextractor
   _instrume=readkey3(hdr,'instrume')  
   check=1
   if not _ZP:     check=0
   if not _gain:   check=0
   if not _fwhm:   check=0
   if not _mbkg:   check=0
   else:  
      if _mbkg<=0: _mbkg=0
   if check==1:
      # formula from McLean 1997)
      from numpy import pi,log10
      if _instrume=='efosc':   ps=readkey3(hdr,'binx')*.12
      else:   ps=0.288
      n=pi*((_fwhm/ps)**2) 
      sn=5       # signal to noise
      maglim=_ZP -2.5 * log10(sn * (1/_gain) * ((n*_mbkg/_exptime)**(.5)) )
      return maglim
   else: return ''
##########################################################################
def marksn2(img,fitstab,frame=1,fitstab2='',verbose=False):
   from pyraf import iraf
   from numpy import array   #,log10
   import lsc
   iraf.noao(_doprint=0)
   iraf.digiphot(_doprint=0)
   iraf.daophot(_doprint=0)
   iraf.images(_doprint=0)
   iraf.imcoords(_doprint=0)
   iraf.proto(_doprint=0)
   iraf.set(stdimage='imt1024')
   hdr=lsc.util.readhdr(fitstab)
   _filter=lsc.util.readkey3(hdr,'filter')
   column=lsc.lscabsphotdef.makecatalogue([fitstab])[_filter][fitstab]

   rasex=array(column['ra0'],float)
   decsex=array(column['dec0'],float)


   if fitstab2:
      hdr=lsc.util.readhdr(fitstab2)
      _filter=lsc.util.readkey3(hdr,'filter')
      _exptime=lsc.util.readkey3(hdr,'exptime')
      column=lsc.lscabsphotdef.makecatalogue([fitstab2])[_filter][fitstab2]
      rasex2=array(column['ra0'],float)
      decsex2=array(column['dec0'],float)

   iraf.set(stdimage='imt1024')
   iraf.display(img + '[0]',frame,fill=True,Stdout=1)
   vector=[]
   for i in range(0,len(rasex)):
      vector.append(str(rasex[i])+' '+str(decsex[i]))

   xy = iraf.wcsctran('STDIN',output="STDOUT",Stdin=vector,Stdout=1,image=img+'[0]',inwcs='world',units='degrees degrees',outwcs='logical',\
                         formats='%10.1f %10.1f',verbose='yes')[3:]
   iraf.tvmark(frame,'STDIN',Stdin=list(xy),mark="circle",number='yes',label='no',radii=10,nxoffse=5,nyoffse=5,color=207,txsize=2)

   if verbose:
 #     print 2.5*log10(_exptime)
      for i in range(0,len(column['ra0'])):
         print xy[i],column['ra0'][i],column['dec0'][i],column['magp3'][i],column['magp4'][i],column['smagf'][i],column['magp2'][i]

   if fitstab2:
      vector2=[]
      for i in range(0,len(rasex2)):
         vector2.append(str(rasex2[i])+' '+str(decsex2[i]))
      xy1 = iraf.wcsctran('STDIN',output="STDOUT",Stdin=vector2,Stdout=1,image=img+'[0]',inwcs='world',units='degrees degrees',outwcs='logical',\
                            formats='%10.1f %10.1f',verbose='yes')[3:]
      iraf.tvmark(frame,'STDIN',Stdin=list(xy1),mark="cross",number='yes',label='no',radii=10,nxoffse=5,nyoffse=5,color=205,txsize=2)


###############################

def Docosmic(img,_sigclip=5.5,_sigfrac=0.2,_objlim=4.5):
   import time
   start=time.time()
   import lsc
   import re,os,string
   import numpy as np
   import tempfile

   ar, hd = fits.getdata(img, header=True)

   if 'TELID' in hd:
      _tel=hd['TELID']
   elif 'telescop' in hd:
      _tel = hd['telescop']
   else:
      _tel='extdata'

   if _tel in ['fts','ftn']:
      temp_file0 = next(tempfile._get_candidate_names())
      lsc.delete(temp_file0)
      out_fits = fits.PrimaryHDU(header=hd,data=ar)
      try:
         out_fits.scale('float32',bzero=0,bscale=1)
      except TypeError as e:
         print "FITS rescaling failed (but it probably doesn't matter). See Astropy Issue #5955."
         print e
      out_fits.writeto(temp_file0, overwrite=True, output_verify='fix')
      ar = fits.getdata(temp_file0)
      lsc.delete(temp_file0)
      gain    = hd['GAIN']
      sat     = 35000
      rdnoise = hd['RDNOISE']
   else:
      if 'gain' in hd:
         gain    = hd['GAIN']
      else:
         print 'warning GAIN not found'
         gain = 1
      if 'saturate' in hd:
         sat     = hd['SATURATE']
      else:
         print 'warning SATURATE not found'
         sat = 60000
      if 'RDNOISE' in hd:
         rdnoise = hd['RDNOISE']
      else:
         print 'warning RDNOISE not found'
         rdnoise = 1
   # need to trick LACosmic into using the right sigma for a sky-subtracted image
   med = np.median(ar)                           # median pixel of image (in ADU)
   noise = 1.4826*np.median(np.abs(ar - med))    # using median absolute deviation instead of sigma
   _pssl = gain*noise**2 - rdnoise**2/gain - med # previously subtracted sky level
   ar[ar < -_pssl] = sat                         # change (what will be) negative values to saturated

#   # same as above but using stats from ORAC pipeline
#   if 'L1SIGMA' in hd and 'L1MEAN' in hd:
#      _pssl=((gain*hd['L1SIGMA'])**2-rdnoise**2)/gain-hd['L1MEAN']
#   else:
#      _pssl=0.0

   print 'gain    sat     noise   sigclip objlim  sigfrac pssl'
   print '{:<7.1f} {:<7.0f} {:<7.1f} {:<7.1f} {:<7.0f} {:<7.1f} {:<7.2f}'.format(gain, sat, rdnoise, _sigclip, _objlim, _sigfrac, _pssl)

   niter = 1
   c = lsc.cosmics.cosmicsimage(ar, pssl=_pssl, gain=gain, readnoise=rdnoise, sigclip=5, sigfrac=0.3 , objlim=5, satlevel=sat)
   c.run(maxiter = niter)

   out=re.sub('.fits','.clean.fits',string.split(img,'/')[-1])
   outmask=re.sub('.fits','.mask.fits',string.split(img,'/')[-1])
   outsat=re.sub('.fits','.sat.fits',string.split(img,'/')[-1])

   out1=c.cleanarray
   out2=c.cleanarray-c.rawarray
   out3=c.getsatstars()

   out_fits = fits.PrimaryHDU(header=hd,data=out1)
   out_fits.writeto(out, overwrite=True, output_verify='fix')

   # we are going to register the mask for the template image,
   # so it makes sense to save it as a float instead of an int
   if 'temp' in img: pixtype = 'float32'
   else:             pixtype = 'uint8'
   out_fits = fits.PrimaryHDU(header=hd,data=(out2!=0).astype(pixtype))
   out_fits.writeto(outmask, overwrite=True, output_verify='fix')

   out_fits = fits.PrimaryHDU(header=hd,data=(out3!=0).astype('uint8'))
   out_fits.writeto(outsat, overwrite=True, output_verify='fix')

   print 'time to do cosmic ray rejection:', time.time()-start
   return out,outmask,outsat

##############################################

def checksnlist(img,listfile):
    import lsc
    import string
    from lsc.util import readkey3,readhdr
    from numpy import cos,sin,arccos,pi, argmin
    scal=pi/180.    
    std,rastd,decstd,magstd=lsc.util.readstandard(listfile)
    hdrt=readhdr(img)
    _ra=readkey3(hdrt,'RA')
    _dec=readkey3(hdrt,'DEC')
    _object=readkey3(hdrt,'object')
    _xdimen=readkey3(hdrt,'XDIM')
    _ydimen=readkey3(hdrt,'YDIM')
    if not _xdimen: _xdimen=readkey3(hdrt,'NAXIS1')
    if not _ydimen: _ydimen=readkey3(hdrt,'NAXIS2')
    dd=arccos(sin(_dec*scal)*sin(decstd*scal)+cos(_dec*scal)*cos(decstd*scal)*cos((_ra-rastd)*scal))*((180/pi)*3600)
    lll=[str(rastd[argmin(dd)])+' '+str(decstd[argmin(dd)])]
    from pyraf import iraf
    bbb=iraf.wcsctran('STDIN','STDOUT',img+'[0]',Stdin=lll,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.5f %10.5f',Stdout=1)[3]
    if 'INDEF' not in bbb and float(string.split(bbb)[0])<=_xdimen and float(string.split(bbb)[1])<=_ydimen and float(string.split(bbb)[0])>=0 and float(string.split(bbb)[1])>=0:
        #print str(std[argmin(dd)])+' in the field '+str(bbb)
        _RA=rastd[argmin(dd)]
        _DEC=decstd[argmin(dd)]
        _SN=std[argmin(dd)]
    else: 
        #print 'out '+str(bbb)
        _RA,_DEC,_SN='','',''
    return _RA,_DEC,_SN

##########################################################################################################
def checksndb(img):
    filename = os.path.basename(img)
    command = '''select targets.ra0, targets.dec0, classificationid
                 from targets, photlcoraw
                 where targetid=targets.id and filename="{}"'''.format(filename)
    target = lsc.mysqldef.query([command], lsc.conn)
    if target:
        return target[0]['ra0'], target[0]['dec0'], target[0]['classificationid']
    else:
        return '', '', ''
##################################################################3

def getcatalog(name_or_filename, field):
    catalog = ''
    catalog_path = lsc.__path__[0] + '/standard/cat/' + field + '/'
    # get the catalog from the database
    entries = lsc.mysqldef.query(['''select name, sloan_cat, landolt_cat, apass_cat
                                     from targets, targetnames, photlco
                                     where targets.id=targetnames.targetid and targets.id=photlco.targetid
                                     and (name="{0}" or filename="{0}")'''.format(os.path.basename(name_or_filename))], lsc.conn)
    if entries:
        if field + '_cat' in entries[0] and entries[0][field + '_cat']:
            catalog = catalog_path + entries[0][field + '_cat']
    # if not in database, search for the catalog in the directory
    if not catalog:
        for entry in entries:
            catlist = os.listdir(catalog_path)
            targetnames = [os.path.basename(cat).split('_' + field)[0].lower() for cat in catlist]
            targetname = entry['name'].replace(' ', '').lower()
            if targetname in targetnames:
               catalog = catalog_path + catlist[targetnames.index(targetname)]
    return catalog

######################################################################################################
