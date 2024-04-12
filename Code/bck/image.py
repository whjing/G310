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
sys.path.append("~/academic/module/")



# from myfits import load_fits_image, crop



# Basic Setup


plt.rcParams['font.size'] = 20


plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内

plt.rc('font',family='Times New Roman') 

# Function

def smooth(data_in, hpbw=2.):
    g = Gaussian2DKernel(hpbw)
    data_out = convolve(data_in, g, boundary='extend')
    return data_out

def show_ctr(ax, data_orig, hdr_ctr, hdr_img, rms, num = 6, hpbw = 2., color = 'tomato', linewidths = 2):
    data_s = smooth(data_orig, hpbw)
    levels=[]
    for i in range(num): 
        level = rms* 5 * (2** (i))
        levels.append(level)
    data_ctr, _ = reproject_interp(pf.PrimaryHDU(data_s, hdr_ctr), hdr_img)
    contour = ax.contour(data_ctr, levels=levels, colors = color, linewidths = linewidths)
    return contour


def show_cb(im, vmin, vmax, stp = 7, n = 1000, label = r'mJy beam$^{-1}$'):
    cb = plt.colorbar(im,  pad=0, fraction=0.0476, location='right')#orientation='horizontal')fraction=0.05,
    colors_range = np.linspace(vmin, vmax, stp)
    cb.set_ticks(colors_range)
    tick_label = []
    for l in range(stp):
        tick_label.append(round(colors_range[l] * n, 1))
    cb.set_ticklabels(tick_label)
    cb.ax.tick_params(labelsize=22)
    cb.set_label(label)
    return cb

def set_ax(ax, num_ticks=5, minpad = 1 ):
    ra = ax.coords[0]
    dec = ax.coords[1]

    ra.set_axislabel('Right Ascension (J2000)',minpad = minpad)
    dec.set_axislabel('Declination (J2000)', minpad = minpad,)
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



def beampatch(hdr, color = 'yellow'):
    pix_size1 = abs(hdr['CDELT1'])
    pix_size2 = abs(hdr['CDELT2'])
    p_pix = [0.005 / pix_size1, 0.005 / pix_size2]
    beam_size_x = hdr['BMAJ'] / pix_size1
    beam_size_y = hdr['BMIN'] / pix_size2
    print("++++++++++++++++++", beam_size_x, beam_size_y)
    beam = Ellipse((p_pix[0], p_pix[1]), beam_size_x, beam_size_y, angle=hdr['BPA'], facecolor='none', edgecolor=color, alpha=1, zorder=200, linewidth=2)

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



# linner scale
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection = wcs)
vmin = np.nanpercentile(data, 1)
vmax = np.nanpercentile(data, 90)
im = ax.imshow(data, origin='lower', cmap = 'viridis', vmin=vmin, vmax=vmax)
beam = beampatch(hdr, color = 'yellow')
ax.add_patch(beam)
set_ax(ax)
# sv_fig(name)
plt.show()


rms = 3e-5

# log scale
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection = wcs)

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

beam = beampatch(hdr, color = 'yellow')
ax.add_patch(beam)
show_cb(im, vmin=vmin, vmax= vmax)

sv_fig("../Output/ASKAP")

show_ctr(ax, data, hdr, hdr, rms, num = 14, hpbw = 2., color = 'white', linewidths = 2)
ax.set_title("ASKAP image")

sv_fig("../Output/ASKAP_ctr")
plt.show()




ASKAP_hdr, ASKAP_data, ASKAP_wcs = hdr, data, wcs

#############################
#### Chandra 

# orig
Chandra_file = '../Fits/G310_i_Chandra_r250.fits'
with pf.open(Chandra_file) as hdul:
    print(hdul)
    hdu = hdul[0]
    hdr = hdu.header
    data = hdu.data
    wcs = WCS(hdr)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection = wcs)
vmin = 5e-9
vmax = np.nanpercentile(data,99.9)
data[np.where(data < vmin)] = vmin
im = ax.imshow(data,  cmap = 'gist_stern', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower') #vmin = vmin, vmax = vmax)
set_ax(ax)
print(vmin, vmax)
# show_cb(im, vmin=vmin, vmax= vmax, n= 1e6, label = '?')
show_ctr(ax, ASKAP_data, ASKAP_hdr, hdr, rms, num = 14, hpbw = 2., color = 'white', linewidths = 2)
ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=500, marker='+', color='cyan' ,transform=ax.get_transform('world'))
# beam = beampatch(hdr, color = 'yellow')
# ax.add_patch(beam)
ax.set_title("Chandra image with ASKAP contour")
sv_fig("../Output/Chandra_o_ctr")

plt.show()



# smooth
ChandraSmooth_file = '../Fits/G310_i_Chandra_r250_smooth.fits'
with pf.open(ChandraSmooth_file) as hdul:
    print(hdul)
    hdu = hdul[0]
    hdr = hdu.header
    data = hdu.data
    wcs = WCS(hdr)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection = wcs)
vmin = 9e-10
vmax = np.nanpercentile(data,99.9)
data[np.where(data < vmin)] = vmin
im = ax.imshow(data,  cmap = 'jet_r', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower') #vmin = vmin, vmax = vmax)
set_ax(ax)
print(vmin, vmax)
# show_cb(im, vmin=vmin, vmax= vmax, n= 1e6, label = '?')
show_ctr(ax, ASKAP_data, ASKAP_hdr, hdr, rms, num = 14, hpbw = 2., color = 'white', linewidths = 2)
beam = beampatch(hdr, color = 'yellow')
ax.add_patch(beam)
ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=500, marker='+', color='cyan' ,transform=ax.get_transform('world'))
ax.set_title("Chandra image (convolved) with ASKAP contour")
sv_fig("../Output/Chandra_smooth_o")

plt.show()


# %%
print(repr(hdr))

# %%
EMU_file = '../Fits/SpecIndex_src.fits'
# hdr, data, wcs = load_fits_image(EMU_file)
with pf.open(EMU_file) as hdul:
    print(hdul)
    hdu = hdul[0]
    hdr = hdu.header
    data = hdu.data
    wcs = WCS(hdr)



# linner scale
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection = wcs)
vmin = -1
vmax = 1
im = ax.imshow(data, origin='lower', cmap = 'RdBu_r', vmin=vmin, vmax=vmax)
beam = beampatch(hdr, color = 'yellow')
ax.add_patch(beam)
set_ax(ax)
show_ctr(ax, ASKAP_data, ASKAP_hdr, hdr, rms, num = 14, hpbw = 2., color = 'red', linewidths = 2)

show_cb(im, vmin=vmin, vmax= vmax, n =1)
sv_fig("../Output/SpecIndex_TT1")
plt.show()
# %%
