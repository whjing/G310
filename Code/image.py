#%%
import os
import re
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
from glob import glob
from matplotlib.patches import Ellipse, Rectangle
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import interpolate
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator
import sys
import matplotlib.ticker as ticker 
from globalSetting import *
set_plot_settings()

def smooth(data_in, hpbw=2.):
    g = Gaussian2DKernel(hpbw)
    data_out = convolve(data_in, g, boundary='extend')
    return data_out

def show_ctr(ax, data_orig, hdr_ctr, hdr_img, rms, num = 6, hpbw = 2., color = 'tomato', linewidths = 2):
    data_s = smooth(data_orig, hpbw)
    levels=[]
    # for i in range(num): 
    #     level = rms* 5 * (2** (i))
    #     levels.append(level)
    # print(levels)
    levels = [0.00015, 0.00030, 0.00040, 0.0006, 0.0012, 0.0024000000000000002, 0.0048000000000000004, 0.009600000000000001, 0.019200000000000002, 0.038400000000000004, 0.07680000000000001, 0.15360000000000001, 0.30720000000000003, 0.6144000000000001, 1.2288000000000001]
    # levels = [0.00015, 0.00021, 0.0003, 0.0004, 0.0005, 0.0006, 0.0012, 0.0015, 0.0024000000000000002, 0.0048000000000000004, 0.009600000000000001, 0.019200000000000002, 0.038400000000000004, 0.07680000000000001, 0.15360000000000001, 0.30720000000000003, 0.6144000000000001, 1.2288000000000001]
    data_ctr, _ = reproject_interp(pf.PrimaryHDU(data_s, hdr_ctr), hdr_img)
    contour = ax.contour(data_ctr, levels=levels, colors = color, linewidths = linewidths)
    return contour


def show_cb(im, vmin, vmax, stp = 6, n = 1000, label = r'mJy beam$^{-1}$'):
    cb = plt.colorbar(im,  pad=0, fraction=0.0476, location='right',)#orientation='horizontal')#fraction=0.05,
    colors_range = np.linspace(vmin, vmax, stp)
    cb.set_ticks(colors_range)
    tick_label = []
    for l in range(stp):
        tick_label.append(round(colors_range[l] * n, 1))
    cb.set_ticklabels(tick_label)
    cb.ax.tick_params(labelsize=22)
    cb.set_label(label)
    return cb

def set_ax(ax, num_ticks=5, minpad = 1 , dec_display = True):
    ra = ax.coords[0]
    dec = ax.coords[1]
    
    ra.set_axislabel('Right Ascension (J2000)',minpad = minpad)
    if dec_display == True:
        dec.set_axislabel('Declination (J2000)', minpad = minpad,)
    else:
        dec.set_axislabel(' ', minpad = minpad,)
    ra.set_major_formatter('dd:mm:ss.s')
    dec.set_major_formatter('dd:mm')
    ax.tick_params(axis='ra', colors='white')
    ax.tick_params(axis='dec', colors='white')
    ra.set_ticklabel(color='black', size=15, pad = 10)
    dec.set_ticklabel(color='black', size=15,)


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

def sv_fig(name):
    sv = plt.savefig(name + '.png', bbox_inches='tight', dpi = 300)
    sv = plt.savefig(name + '.pdf', bbox_inches='tight')
    return sv



def beampatch(hdr, color = 'white'):
    pix_size1 = abs(hdr['CDELT1'])
    pix_size2 = abs(hdr['CDELT2'])
    p_pix = [0.005 / pix_size1, 0.005 / pix_size2]
    beam_size_x = hdr['BMAJ'] / pix_size1
    beam_size_y = hdr['BMIN'] / pix_size2
    print("++++++++++++++++++", beam_size_x, beam_size_y)
    beam = Ellipse((p_pix[0], p_pix[1]), beam_size_x, beam_size_y, angle=hdr['BPA'], facecolor=color, edgecolor=None, alpha=1, zorder=200, linewidth=2)

    res=np.sqrt(hdr['BMAJ'] *hdr['BMIN'])*3600
    print('resolution'+str(res))
    return beam

##################################
###### EMU data 
##################################

EMU_file = '../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-250.fits'
# hdr, data, wcs = load_fits_image(EMU_file)
with pf.open(EMU_file) as hdul:
    print(hdul)
    hdu = hdul[0]
    hdr = hdu.header
    data = hdu.data
    wcs = WCS(hdr)



rms = 3e-5

# log scale
fig = plt.figure(figsize=(30,10))
ax = fig.add_subplot(131, projection = wcs)

vmin = rms * 3 
vmax = np.nanpercentile(data,100)
data[np.where(data < vmin)] = vmin
im = ax.imshow(data,  cmap = 'jet', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower') #vmin = vmin, vmax = vmax)
set_ax(ax)
print(vmin, vmax)
# PSR 
c_PSR=('14:00:45.69',' -63:25:42.6')
c_PSR= SkyCoord(ra=c_PSR[0], dec=c_PSR[1], frame='fk5',unit=(u.hourangle, u.deg))
ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=500, marker='+', color='cyan' ,transform=ax.get_transform('world'))
print(c_PSR.ra.deg, c_PSR.dec.deg,)

beam = beampatch(hdr, color = 'white')
ax.add_patch(beam)
show_cb(im, vmin=vmin, vmax= vmax)

# sv_fig("../figures/ASKAP")

show_ctr(ax, data, hdr, hdr, rms, num = 14, hpbw = 2., color = 'white', linewidths = 2)
# ax.set_title("ASKAP image")

# sv_fig("../figures/ASKAP_ctr")
# plt.show()




ASKAP_hdr, ASKAP_data, ASKAP_wcs = hdr, data, wcs

#############################
#### Chandra 

# orig
Chandra_file = '../Fits/G310_i_Chandra_r250.fits'
with pf.open(Chandra_file) as hdul:
    hdu = hdul[0]
    hdr = hdu.header
    data = hdu.data
    wcs = WCS(hdr)

ax = fig.add_subplot(132, projection = wcs)
vmin = 5e-9
vmax = np.nanpercentile(data,99.9)
data[np.where(data < vmin)] = vmin
# in sqrt scale
# im = ax.imshow(data,  cmap = 'jet', norm=colors.PowerNorm(gamma = 0.5, vmin=vmin, vmax=vmax), origin='lower')
im = ax.imshow(data,  cmap = 'jet', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower')
set_ax(ax, dec_display=False)
print(vmin, vmax)
# show_cb(im, vmin=vmin, vmax= vmax, n= 1e6, label = '?')
show_ctr(ax, ASKAP_data, ASKAP_hdr, hdr, rms, num = 14, hpbw = 2., color = 'white', linewidths = 2)
ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=500, marker='+', color='cyan' ,transform=ax.get_transform('world'))
show_cb(im, vmin=vmin, vmax= vmax, n = 1e6,label=f"")
# beam = beampatch(hdr, color = 'yellow')
# ax.add_patch(beam)
# ax.set_title("Chandra image with ASKAP contour")
# sv_fig("../figures/Chandra_ctr")

# plt.show()



# smooth
ChandraSmooth_file = '../Fits/G310_i_Chandra_r250_smooth.fits'
with pf.open(ChandraSmooth_file) as hdul:
    print(hdul)
    hdu = hdul[0]
    hdr = hdu.header
    data = hdu.data
    wcs = WCS(hdr)

ax = fig.add_subplot(133, projection = wcs)
vmin = 9e-10
vmax = np.nanpercentile(data,99.9)
data[np.where(data < vmin)] = vmin
im = ax.imshow(data,  cmap = 'jet_r', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower') #vmin = vmin, vmax = vmax)
set_ax(ax, dec_display=False)
print(vmin, vmax)

show_ctr(ax, ASKAP_data, ASKAP_hdr, hdr, rms, num = 14, hpbw = 2., color = 'white', linewidths = 2)
beam = beampatch(hdr, color = 'white')
ax.add_patch(beam)
ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=500, marker='+', color='cyan' ,transform=ax.get_transform('world'))
# ax.set_title("Chandra image (convolved) with ASKAP contour")
# 
show_cb(im, vmin=vmin, vmax= vmax, n = 1e6, label=f"")

# sv_fig("../figures/image")
plt.show()




# %%
