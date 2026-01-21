#!/usr/bin/env python

import numpy as np
#import matplotlib
#matplotlib.rcParams["font.size"]=20
#matplotlib.rcParams["text.usetex"]=True
#import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm

ax_range = [150,-150,-150,150]

f = open("image.out")
f.readline()
npx, npy = [int(t) for t in f.readline().split()]
nlam = float(f.readline())
au = 1.49597871e13
dpx, dpy = [float(t)/au for t in f.readline().split()] # dpx, dpy in AU
f.readline()
data = np.loadtxt(f)
f.close()

I = data[:,0]
Q = data[:,1]
U = data[:,2]
V = data[:,3]

I.shape = npx, npy
Q.shape = npx, npy
U.shape = npx, npy
V.shape = npx, npy

x = (np.arange(npx) - npx//2 - ((npx+1)%2)*0.5) * dpx
y = (np.arange(npy) - npy//2 - ((npy+1)%2)*0.5) * dpy
XX,YY = np.meshgrid(x,y)

Bma, Bmn = 1.039104329215E-05*3600, 8.840234950185E-06*3600 # arcsec
distance = 140 # pc = 63 light year
ma,mn = Bma*distance, Bmn*distance # AU
#ma,mn = 20, 10 # AU
Bpa =  (2.576065635681E+01+90) / 180.*np.pi # position angle of the beam, in radian
a = np.cos(-Bpa)**2/(2*(0.5*ma)**2) + np.sin(-Bpa)**2/(2*(0.5*mn)**2);
b = -np.sin(-2*Bpa)/(4*(0.5*ma)**2) + np.sin(-2*Bpa)/(4*(0.5*mn)**2);
c = np.sin(-Bpa)**2/(2*(0.5*ma)**2) + np.cos(-Bpa)**2/(2*(0.5*mn)**2);
A = 1/(2*np.pi*0.5*ma*0.5*mn)
kernel = A*np.exp( - (a*XX**2 + 2*b*XX*YY + c*YY**2))
kernel /= kernel.sum()

Omega_beam = np.pi*(Bma/2 /3600/180*np.pi)*(Bmn/2 /3600/180*np.pi) # ster
I_to_Jy = Omega_beam * 1e23  # F_mJy = I_cgs * I_to_mJy

from astropy.io import fits

# Write to FITS
fits.writeto('stokes_I.fits', I_to_Jy*I, overwrite=True)
fits.writeto('stokes_Q.fits', I_to_Jy*Q, overwrite=True)
fits.writeto('stokes_U.fits', I_to_Jy*U, overwrite=True)
fits.writeto('stokes_V.fits', I_to_Jy*V, overwrite=True)

# Stack into cube: shape (4, ny, nx) for IQUV (or just IQU)
stokes_cube = np.stack([I_to_Jy*I, I_to_Jy*Q, I_to_Jy*U, I_to_Jy*V])  # Shape (3, ny, nx)

# Define header with appropriate Stokes axis
hdr = fits.getheader('stokes_I.fits')
hdr['NAXIS'] = 4
hdr['NAXIS4'] = 3
hdr['CTYPE4'] = 'STOKES'
hdr['CRVAL4'] = 1  # 1: I, 2: Q, 3: U, etc.
hdr['CDELT4'] = 1
hdr['CRPIX4'] = 1

fits.writeto('model_cube.fits', stokes_cube, hdr, overwrite=True)

