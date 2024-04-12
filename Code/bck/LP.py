#!/usr/bin/env python
# %%

from astropy.wcs.wcs import WCS
from matplotlib.pyplot import plot, sca
import numpy as np
from astropy.io import fits as pf
from astropy import wcs
from matplotlib import pylab as plt
from reproject import reproject_interp

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from astropy import units as u
from copy import deepcopy
from regions import PixCoord
from regions import EllipseSkyRegion, EllipsePixelRegion
from astropy.coordinates import SkyCoord
import matplotlib.colors as colors
import math
import sys
import pyregion
sys.path.append("/Users/jing/academic/module/")
print(sys.path)
from myfits import *

config = {
    "font.family":'serif',
    "font.size": 22,
    "mathtext.fontset":'stix',
    "font.serif": ['SimSun'], # simsun字体中文版就是宋体
}
plt.rcParams.update(config)


plt.rcParams['font.size'] = 22


plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内

plt.rc('font',family='Times New Roman')


#-------   general ----------#
# plt.rcParams.update({'font.size': 18})
def set_ax(ax, num_ticks=5, minpad = 1 ):
    ra = ax.coords[0]
    dec = ax.coords[1]

    ra.set_axislabel('Right Ascension (J2000)',minpad = minpad)
    dec.set_axislabel('Declination (J2000)', minpad = minpad,)
    ra.set_major_formatter('dd:mm:ss.s')
    dec.set_major_formatter('dd:mm')
    ax.tick_params(axis='ra', colors='white')
    ax.tick_params(axis='dec', colors='white')
    ra.set_ticklabel(color='black', size=22, pad = 10)
    dec.set_ticklabel(color='black', size=22,)


    ra.display_minor_ticks(True)
    dec.display_minor_ticks(True)
    ra.set_minor_frequency(10)
    dec.set_minor_frequency(10)
    

    # ax.spines['top'].set_color('white') 
    # ax.spines['right'].set_color('white')
    # ax.spines['left'].set_color('black')


    # tick_spacing = 0.04
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(tick_spacing))
    
    
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_ticks(number = num_ticks)
    dec.set_ticks(number = num_ticks)
    dec.set_ticklabel_visible(True)
    dec.set_ticklabel(rotation='vertical', va='center')

I_hdr,I_data, I_w = load_fits_image('../Fits/G310_i_EMU_r250.fits')
PI_hdr, PI_data, PI_w = load_fits_image('../Fits/FDF_maxPI_crop.fits')
RM_hdr, RM_data, RM_w = load_fits_image('../Fits/FDF_peakRM_crop.fits')
hdu = pf.open('../Fits/FDF_peakRM_crop.fits')[0]

I_data = I_data
PI_data = PI_data
# PA_data = PA_data
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, projection=I_w)

vmin = 0
#np.percentile(I_data,20)
vmax = np.percentile(I_data,99.3)

im = ax.imshow(I_data, origin = 'lower',cmap='viridis', vmin = vmin, vmax = vmax)
show_cb(im, vmin, vmax)

set_ax(ax)
# sv_fig('kkbr_i')
# polarization

box_s = f"fk5;polygon(210.1910460,-63.4149374,210.1602719,-63.4234925,210.1440260,-63.4398608,210.1785609,-63.4432171,210.2160066,-63.4340989,210.2334541,-63.4167919)"
box = pyregion.parse(box_s)
data_box = box.get_mask(hdu=hdu,shape=RM_data.shape)
print(data_box)


maskcondition = (I_data < 25e-6 * 5) | (PI_data < 100e-6) | (data_box == False)

PI_data[maskcondition] = np.nan
# PI_data[np.where(Frac_data < 0.005)] = np.nan
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, projection=PI_w)
vmin =  np.nanpercentile(PI_data,0)
vmax =  np.nanpercentile(PI_data,99)
print(vmax)
im = ax.imshow(PI_data, origin='lower',cmap='viridis', vmin=vmin, vmax=vmax)
ax.patch.set_facecolor("mintcream")  
set_ax(ax)
show_ctr(ax, I_data, I_hdr, PI_hdr, 5e-5, num = 12, hpbw = 1.5, color = 'gray', linewidths = 1)
show_ctr(ax, PI_data, PI_hdr, PI_hdr, 7.56184e-06, num = 6,hpbw = 1.5, color = 'red', linewidths = 1)

c_PSR=('14:00:45.69',' -63:25:42.6')
c_PSR= SkyCoord(ra=c_PSR[0], dec=c_PSR[1], frame='fk5',unit=(u.hourangle, u.deg))
ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=300, marker='+', color='tomato' ,transform=ax.get_transform('world'),linewidth = 5)
print(c_PSR.ra.deg, c_PSR.dec.deg,)

show_cb(im, vmin, vmax, n=1e6, label = r'$\rm{\mu Jy ~beam^{-1}}$')
# sv_fig('../Output/PI')



fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, projection=PI_w)

Frac_data = PI_data/I_data

Frac_data[maskcondition] = np.nan
# Frac_data[np.where(data_box == False)] = np.nan
# Frac_data[np.where(Frac_data < 0)] = np.nan
vmin, vmax = 0.0, 0.05
im = ax.imshow(Frac_data, origin = 'lower',cmap = 'rainbow_r', vmin = vmin, vmax = vmax) 
show_cb(im, vmin, vmax, n = 100, label = r'%')
show_ctr(ax, PI_data, PI_hdr, PI_hdr, 7.56184e-06, num = 6,hpbw = 1.5, color = 'red', linewidths = 1)
set_ax(ax)
ax.patch.set_facecolor("aliceblue")  
# sv_fig('../Output/frac')


# RM_data = data_box * RM_data
# print(RM_data)
# RM_data[np.where(data_box == False)] = np.nan
RM_data[maskcondition] = np.nan

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, projection=RM_w)
vmin = -800
vmax = -500

im = ax.imshow(RM_data, origin = 'lower',cmap='Spectral', vmin = vmin, vmax = vmax)
show_cb(im, vmin, vmax, n=1, )
show_ctr(ax, PI_data, PI_hdr, PI_hdr, 7.56184e-06, num = 8,hpbw = 1.5, color = 'red', linewidths = 1)
set_ax(ax)
ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=300, marker='+', color='tomato' ,transform=ax.get_transform('world'), linewidth = 5)
ax.patch.set_facecolor("silver")  

# ax.set_xlim(30,210)
# ax.set_ylim(30,210)
sv_fig('../Output/RM')
%%



# %%
