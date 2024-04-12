#%%
import os
import re
import sys
sys.path.append("/Users/jing/academic/module/")
print(sys.path)
print("What happend?")
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from glob import glob
from astropy.io import fits as pf
import numpy as np
from astropy.wcs import WCS
from astropy import wcs
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from reproject import reproject_interp
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse, Rectangle
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import interpolate
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator
from copy import deepcopy
from myfits import load_fits_image, crop
from astropy.nddata import Cutout2D
#%%


def crop(header0, data0, RA, DEC, const, rotate = None, writeto = False ):
    # Define the header and data
    hdu = pf.ImageHDU(data0, header=header0)
    hdr0 = deepcopy(header0)
    size = round(abs(const/hdr0['CDELT1']))
    hdr0['NAXIS1'] = size  # 右
    hdr0['NAXIS2'] = size  # 上
    hdr0['CRVAL1'] = RA
    hdr0['CRVAL2'] = DEC
    hdr0['CRPIX1'] = size/2  # 左
    hdr0['CRPIX2'] = size/2  # 下
    hdr0['CROTA1'] = 0  # 旋转前
    if rotate != None:
        hdr0['CROTA2'] = rotate  # 旋转后
    data, footprint = reproject_interp(hdu, hdr0)
    wcs = WCS(hdr0)
    if writeto != False:
        pf.writeto(writeto, data, hdr0, overwrite = True)
    return hdr0, data, wcs


hdr = pf.getheader(file)
wcs = WCS(hdr).dropaxis(2)
print(wcs)
data = pf.getdata(file)

rms = 5.11448088218458e-05

# the center of SNR
ra0 = '14:00:46' 
dec0 = '-63:25:43' 
c0 = SkyCoord(ra0, dec0, frame='icrs', unit=('hour', 'degree'))
ra0, dec0 = c0.ra.deg, c0.dec.deg
print(ra0, dec0)

hdr, data, wcs = crop(hdr, data, ra0, dec0, 500/3600)
fig = plt.figure()

wcs = wcs.dropaxis(2)
print(wcs)
ax = fig.add_subplot(1,1,1, projection = wcs)
vmin = np.nanpercentile(data, 10)
vmax = np.nanpercentile(data, 95)
im = ax.imshow(data, origin='lower', cmap = 'viridis', vmin=vmin, vmax=vmax) 
outname = file.replace('.fits','_crop-500.fits')
print(outname)

# pf.writeto(outname, data, hdr, overwrite=True)
print("The file has been restored as {outname }")

#%%
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
position = SkyCoord('13h11m29.96s -01d19m18.7s', frame='icrs')
wcs = WCS(naxis=2)
rho = np.pi / 3.
scale = 0.05 / 3600.
wcs.wcs.cd = [[scale*np.cos(rho), -scale*np.sin(rho)],
              [scale*np.sin(rho), scale*np.cos(rho)]]
wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
wcs.wcs.crval = [position.ra.to_value(u.deg),
                 position.dec.to_value(u.deg)]
wcs.wcs.crpix = [50, 100]
cutout = Cutout2D(data, position, (30, 40), wcs=wcs)
plt.imshow(cutout.data, origin='lower')
#%%
#######################
### The problem is "large_array_shape" and "small_array_shape" must have the same number of dimensions.
########################

def cutout_largefile(filename, ra0, dec0, r_deg):
    hdu = pf.open(filename)[0]
    hdr = hdu.header
    wcs = WCS(hdr)
    

    position = SkyCoord(ra0, dec0, frame='icrs', unit =('hour', 'degree'))
    print(position)
    size = round(abs(r_deg/hdr['CDELT1']))

    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data

    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    cutout_filename = 'example_cutout.fits'
    hdu.writeto(cutout_filename, overwrite=True)
if __name__ == "__main__":
    filename = "../Fits/image.i.EMU_1356-64.SB53310.beam27.spw-1.taylor.0.restored.fits"
    ra0 = '14:00:46' 
    dec0 = '-63:25:43' 
    r_deg = 500/3600
    cutout_largefile(filename, ra0, dec0, r_deg)




# %%
