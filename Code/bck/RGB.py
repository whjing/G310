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
        level = rms* 3 * (2** (i))
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

file_EMU = '../Fits/G310_i_EMU_r250.fits'
# hdr, data, wcs = load_fits_image(EMU_file)
hdu_ASKAP = pf.open(file_EMU )[0]
hdr_ASKAP = hdu_ASKAP.header
data_ASKAP = hdu_ASKAP.data
wcs_ASKAP = WCS(hdr_ASKAP)


file_Chandra = '../Fits/G310_i_Chandra_r250.fits'
# hdul =  pf.open(Chandra_file)
hdu_Chandra = pf.open(file_Chandra)[0]
hdr_Chandra = hdu_Chandra.header
data_Chandra = hdu_Chandra.data
wcs_Chandra = WCS(hdr_Chandra)



# log scale
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection = wcs_Chandra)

vmin = 1.2e-4
vmax = np.nanpercentile(data_ASKAP,100)
data_ASKAP[np.where(data_ASKAP < vmin)] = vmin
data_ASKAP, _ = reproject_interp(hdu_ASKAP, hdr_Chandra)
im = ax.imshow(data_ASKAP,  cmap = 'Reds_r', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower', alpha = 1) #vmin = vmin, vmax = vmax)
set_ax(ax)
print(vmin, vmax)
# PSR 
# c_PSR=('14:00:45.69',' -63:25:42.6')
# c_PSR= SkyCoord(ra=c_PSR[0], dec=c_PSR[1], frame='fk5',unit=(u.hourangle, u.deg))
# ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=500, marker='+', color='cyan' ,transform=ax.get_transform('world'))
# print(c_PSR.ra.deg, c_PSR.dec.deg,)

# beam = beampatch(hdr, color = 'yellow')
# ax.add_patch(beam)
# show_cb(im, vmin=vmin, vmax= vmax)
# sv_fig("../Output/ASKAP")
# plt.show()


#############################
#### Chandra 

# orig


# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(111, projection = wcs)
vmin = 5e-9
vmax = np.nanpercentile(data_Chandra,99.9)
data_Chandra[np.where(data_Chandra < vmin)] = vmin

im = ax.imshow(data_Chandra,  cmap = 'Blues_r', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower', alpha = 0.7)#vmin = vmin, vmax = vmax)
set_ax(ax)
print(vmin, vmax)

plt.show()



# linner scale
# %%
