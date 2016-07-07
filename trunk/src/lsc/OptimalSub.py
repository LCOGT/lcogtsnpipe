#!/home/dguevel/conda-envs/david/bin python

import numpy as np
from scipy.optimize import minimize
from astropy.io import fits
from pyraf import iraf
import time

def Subtract(Nf,Rf,Pnf,Prf,Df):

    #open fits files
    N_hdulist=fits.open(Nf)
    R_hdulist=fits.open(Rf)

    N=N_hdulist[0].data
    R=R_hdulist[0].data

    #fix and roll PSFs, crop the images to the same shape
    Pn=ExtractPSF(Pnf)
    Pr=ExtractPSF(Prf)

    N,R,Pn,Pr=CropImages(N,R,Pn,Pr)


    print 'Doing Fourier transforms...'

    sn=np.std(N)
    sr=np.std(R)

    print '1/4'
    N_hat=np.fft.fft2(N)    
    print '2/4'
    R_hat=np.fft.fft2(R)
    print '3/4'
    Pr_hat=np.fft.fft2(Pr)
    print '4/4'
    Pn_hat=np.fft.fft2(Pn)

    print 'Fitting B, gamma...'
    #fit B using robust fitting
    t0=time.time()

    #N=np.subtract(N,Background(N))
    #R=np.subtract(N,Background(R))

    B,gamma=FitB(N_hat,R_hat,Pn_hat,Pr_hat,sn,sr)

    print 'B = '+str(B)
    print 'Gamma = '+str(gamma)
    print 'Fit time: '+str(time.time()-t0)
    print 'Subtracting images...'

    N=np.subtract(N,gamma)

    #subtract images
    D_hat=(Pr_hat*N_hat-B*Pn_hat*R_hat)/np.sqrt(sn**2*abs(Pr_hat**2)+sr**2*B**2*abs(Pn_hat**2))
    D=np.real(np.fft.ifft2(D_hat))

    hdu = fits.PrimaryHDU(D)
    hdu.header=N_hdulist[0].header
    hdu.header['CONVOL00']='TEMPLATE'
    hdu.writeto(Df)

    print 'Done!'
    return D

def Background(im,num=100):
    back=[]
    s=np.std(im)
    ind=[int(im.shape[0]*np.random.rand(1)),int(im.shape[1]*np.random.rand(1))]
    for i in range(num):
        ind=[int(im.shape[0]*np.random.rand(1)),int(im.shape[1]*np.random.rand(1))]
        n=im[ind[0],ind[1]]
        if abs(n)<=np.mean(im)+s:
            back.append(n)
    bg=np.mean(back)
    return bg

def Diff(Bg,N_hat,R_hat,Pn_hat,Pr_hat,sn,sr):
    B,gamma=Bg
    Dn_hat=Pr_hat*N_hat/np.sqrt(sn**2*abs(Pr_hat**2)+B**2*sr**2*abs(Pn_hat)**2)
    Dr_hat=Pn_hat*R_hat/np.sqrt(sn**2*abs(Pr_hat**2)+B**2*sr**2*abs(Pn_hat)**2)
    Dn=np.real(np.fft.ifft2(Dn_hat))
    Dr=np.real(np.fft.ifft2(Dr_hat))
    gp=gamma/np.sqrt(sn**2+B**2*sr**2)
    eqn=(B*Dr+gp-Dn)

    return 0.5*np.sum(eqn**2)

def FitB(N_hat,R_hat,Pn_hat,Pr_hat,sn,sr):
    #nonlinear robust least squares
    bound=[(.1,5.),(-2000.,2000.)]
    status=False
    mi=1000
    n=N_hat.shape
    a,b=[n[0]-50,n[1]+50]
    N_hat,R_hat,Pn_hat,Pr_hat,sn,sr=N_hat[a:b,a:b],R_hat[a:b,a:b],Pn_hat[a:b,a:b],Pr_hat[a:b,a:b],sn,sr
    m=minimize(Diff,[1.,0.],args=(N_hat,R_hat,Pn_hat,Pr_hat,sn,sr),bounds=bound,options={'maxiter':mi})
    status=m.success
    print m
    B,gamma=m.x
    return B,gamma


def CropImages(N,R,Pn,Pr):
    if N.shape[0]<=R.shape[0] & N.shape[1]<=R.shape[1]:
        s=N.shape
    elif N.shape[0]>R.shape[0] & N.shape[1]>R.shape[1]:
        s=R.shape
    else:
        s=[np.min([N.shape[0],R.shape[0]]),np.min([N.shape[1],R.shape[1]])]

    N=N[:s[0],:s[1]]
    R=R[:s[0],:s[1]]
    p=Pn.shape

    Pn_ext=np.zeros(s)
    Pn_ext[s[0]/2-p[0]/2-1:s[0]/2+p[0]/2,s[1]/2-p[1]/2-1:s[1]/2+p[1]/2]=Pn
    Pn_ext=np.roll(Pn_ext,s[0]/2,0)
    Pn_ext=np.roll(Pn_ext,s[1]/2,1)

    p=Pr.shape

    Pr_ext=np.zeros(s)
    Pr_ext[s[0]/2-p[0]/2-1:s[0]/2+p[0]/2,s[1]/2-p[1]/2-1:s[1]/2+p[1]/2]=Pr
    Pr_ext=np.roll(Pr_ext,s[0]/2,0)
    Pr_ext=np.roll(Pr_ext,s[1]/2,1)

    return N,R,Pn_ext,Pr_ext

def ExtractPSF(Pf):

    iraf.noao()
    iraf.digiphot()
    iraf.daophot()
    iraf.seepsf(Pf,'temp.psf.fits')
    P=fits.open('temp.psf.fits')[0].data
    system('rm temp.psf.fits')

    #normalize PSF
    P=P/np.sum(P)

    return P





