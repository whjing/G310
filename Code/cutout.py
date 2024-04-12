#%%
# import

import os
import re
import sys
sys.path.append("/Users/jing/academic/module/")
print(sys.path)
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from glob import glob
from astropy.io import fits as pf
import numpy as np
from astropy.wcs import WCS
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
# SETUP




plt.rcParams['font.size']
plt.rcParams['font.size'] = 20

plt.rc('font',family='Times New Roman') 
# Function

def modify_hdr(header):
    stripHeadKeys = ['NAXIS','CRVAL','CRPIX','CDELT','CTYPE','CROTA',
                         'CD1_','CD2_','CUNIT','PC1_','PC2_','PC3_','PC4_']
    stripHeadKeys1 = ['PC1_','PC2_','PC3_','PC4_','PC01_','PC02_','PC03_','PC04_']
    for i in range(3,6):
        for key in stripHeadKeys:
            if (key+str(i)) in header: del header[key+str(i)]
    for j in range(1,6):            
        for key1 in stripHeadKeys1:
            if (key1+str(j)) in header: del header[key1+str(j)]
            if (key1+ '0' + str(j)) in header: del header[key1+ '0' + str(j)]
    header['NAXIS'] = 2
    return header



def smooth(data_in, hpbw=2.):
    g = Gaussian2DKernel(hpbw)
    data_out = convolve(data_in, g, boundary='extend')
    return data_out


def beampatch(hdr, color = 'yellow'):
    pix_size1 = abs(hdr['CDELT1'])
    pix_size2 = abs(hdr['CDELT2'])
    p_pix = [0.005 / pix_size1, 0.005 / pix_size2]
    beam_size_x = hdr_askap['BMAJ'] / pix_size1
    beam_size_y = hdr_askap['BMIN'] / pix_size2
    beam = Ellipse((p_pix[0], p_pix[1]), beam_size_x, beam_size_y, angle=hdr_askap['BPA'], facecolor='none', edgecolor=color, alpha=1, zorder=200, linewidth=2)
    res=np.sqrt(hdr_askap['BMAJ'] *hdr_askap['BMIN'])*3600
    print('resolution'+str(res))
    return beam


def show_cb(im, vmin, vmax, stp = 7, n = 1000, label = r'mJy beam$^{-1}$'):
    cb = plt.colorbar(im,  pad=0, fraction=0.0476, location='top')#orientation='horizontal')fraction=0.05,
    colors_range = np.linspace(vmin, vmax, stp)
    cb.set_ticks(colors_range)
    tick_label = []
    for l in range(stp):
        tick_label.append(round(colors_range[l] * n, 1))
    cb.set_ticklabels(tick_label)
    cb.ax.tick_params(labelsize=12)
    cb.set_label(label)
    return cb


def set_ax(ax, num_ticks=5, minpad = 1 , ylabel = 'Declination (J2000)' ):
    ra = ax.coords[0]
    dec = ax.coords[1]
    ra.set_axislabel('Right Ascension (J2000)',minpad = minpad)
    dec.set_axislabel(ylabel , minpad = minpad,) 
    plt.xticks(fontproperties = 'Times New Roman', size = 14)
    plt.yticks(fontproperties = 'Times New Roman', size = 14)

    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_ticks(number = num_ticks)
    ra.set_ticklabel_visible(False)
    dec.set_ticks(number = num_ticks)
    dec.set_ticklabel_visible(False)
    dec.set_ticks_position('right')
    # dec.set_ticklabel(rotation='vertical', va='center')

def show_ctr(ax, data_orig, hdr_ctr, hdr_img, rms, num = 6, hpbw = 2., color = 'tomato', linewidths = 2):
    data_s = smooth(data_orig, hpbw)
    levels=[]
    for i in range(num): 
        level = rms* 3 * (2** (i))
        levels.append(level)
    data_ctr, _ = reproject_interp(pf.PrimaryHDU(data_s, hdr_ctr), hdr_img)
    contour = ax.contour(data_ctr, levels=levels, colors = color, linewidths = linewidths)
    return contour

def sv_fig(name):
    sv = plt.savefig(name + '.png', bbox_inches='tight', dpi = 300)
    sv = plt.savefig(name + '.pdf', bbox_inches='tight')
    return sv

# Parameters

rms=5.11448088218458e-05

# the center of SNR
ra0 = '14:00:46' 
dec0 = '-63:25:43' 
c0 = SkyCoord(ra0, dec0, frame='icrs', unit=('hour', 'degree'))
ra0, dec0 = c0.ra.deg, c0.dec.deg
print(ra0, dec0)

print('YES! ALL import!')

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


def crop_cube(RA, Dec, r_deg):
    hdu = pf.open(fits)[0]
    hdr = pf.getheader(fits)
    data = pf.getdata(fits)
    wcs = WCS(hdr)
    hdr0 = deepcopy(hdr)
    size = round(abs(r_deg/hdr0['CDELT1']))
    hdr0['NAXIS1'] = size  # 右
    hdr0['NAXIS2'] = size  # 上
    hdr0['CRVAL1'] = RA
    hdr0['CRVAL2'] = Dec
    hdr0['CRPIX1'] = size/2  # 左
    hdr0['CRPIX2'] = size/2  # 下

    print(data.shape)
    print(repr(hdr0))

    data, footprint = reproject_interp(hdu, hdr0)
    outname = fits.replace('.fits','_crop.fits')
    # pf.writeto(outname, data, hdr0, overwrite = True)





# fitspath = '../Fits/'
# fitsname = 'image.*.cube.cutout.total_3D.fits'
# fitslist = glob(fitspath + fitsname)
# fitsname = "../Chandra_process/process/flux_1/hard_flux.img.fits"
fitsList = glob("../Chandra_process/process/flux_1/hard_flux.img.fits")

# fitspath = '../Fits_4-subband-useful'
# fitsname = 'shell-Chandra.mask.fits'
# fitsList = glob(f"{fitspath}/{fitsname}")
print(fitsList)
RA = ra0
Dec = dec0
r_deg = 250

for fits in fitsList:
    # crop_cube(RA, Dec, r_deg)
    print(fits)
    # # hdu = pf.open(fits)[1]
    # print()
    hdr = pf.getheader(fits)
    print(hdr)
    
    hdr, data, wcs = load_fits_image(fits)
    hdr, data, wcs = crop(hdr, data, ra0, dec0, r_deg/3600)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = wcs)
    vmin = np.nanpercentile(data, 10)
    vmax = np.nanpercentile(data, 95)
    im = ax.imshow(data, origin='lower', cmap = 'viridis', vmin=vmin, vmax=vmax) 
    set_ax(ax)
    outname = fits.replace(f".fits",f"_crop-{r_deg}.fits")
    print(outname)

    pf.writeto(outname, data, hdr, overwrite=True)
    print(f"The file has been restored as {outname}")


# %%
