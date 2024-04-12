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
from matplotlib.patches import Annulus
import pyregion
from sklearn.preprocessing import StandardScaler
from matplotlib.ticker import ScalarFormatter, FuncFormatter
sys.path.append("~/academic/module/")



# from myfits import load_fits_image, crop



# Basic Setup
config = {
    "font.family":'serif',
    "font.size": 20,
    "mathtext.fontset":'stix',
    "font.serif": ['SimSun'], # simsun字体中文版就是宋体
}
plt.rcParams.update(config)


plt.rcParams['font.size']
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
    cb.ax.tick_params(labelsize=12)
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


def normalize_and_scale(data, data_err):
    ymin = np.min(data)
    ymax = np.max(data)
    normalized_data = [(y_i - ymin) / (ymax - ymin) + 1e-6 for y_i in data]
    # norm_data_err = 
    # scaled_data = [y_i * ymax * 1000 for y_i in normalized_data]
    return normalized_data
##################################
###### EMU data 
##################################


# open fits file 

EMU_file = '../Fits/G310_i_Chandra_r250_smooth.fits'

with pf.open(EMU_file) as hdul:
    print(hdul)
    hdu = hdul[0]
    hdr = hdu.header
    data = hdu.data
    wcs = WCS(hdr)

fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(111, projection = wcs)

# Mask Data
box_s = f'fk5;polygon(14:00:54.0800,-63:23:09.825,14:00:14.5310,-63:24:46.807,14:00:45.6449,-63:25:45.558)'
box = pyregion.parse(box_s)
data_box = box.get_mask(hdu=hdu,shape=data.shape)
data = data_box * data


vmin = 5e-11
vmax = np.nanpercentile(data,99.9)
# data[np.where(data < vmin)] = vmin
im = ax1.imshow(data,  cmap = 'jet', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower') #vmin = vmin, vmax = vmax)
set_ax(ax1)
show_cb(im, vmin=vmin, vmax= vmax)

c_0 = box_s.split(",")[-2:]
print(c_0)
ra = c_0[0]
dec = c_0[1].replace(")", "")
print(ra,dec)

ra = "210.1882728"
dec = "-63.4280232"
r = 120

shell_c = SkyCoord(ra, dec, frame='fk5',unit=(u.deg, u.deg))
x, y = wcs.wcs_world2pix(shell_c.ra.deg, shell_c.dec.deg, 0)
r_out = r/((3600*(abs(hdr['CDELT1']))))
r_list = np.linspace(0, r, 30)
r_list_ring = np.linspace(0, r_out, 30)
print("r_out",r_out)

data_old = np.zeros_like(data)
data_profile = []
data_err = []
r_old = 0


for radius in r_list:

    reg =f'fk5;circle({ra},{dec},{radius}")'
    outcircle = pyregion.parse(reg)
    
    # circel_reg = outcircle.as_imagecoord(hdr).get_mpl_patches_texts()
    # print(circel_reg)
    # ax1.add_patch(circel_reg)
    r_out = radius/((3600*(abs(hdr['CDELT1']))))
    color = 'red'
    if radius == 100:
        color = 'pink'

    ring = Circle((x,y), r_out, facecolor = 'none', edgecolor = color, alpha=1, zorder=200)
    ax1.add_patch(ring)
    
    
    data_reg = outcircle.get_mask(hdu=hdu,shape=data.shape)
    data_mask = data_reg * data
    data_shell = data_mask - data_old
    data_shell[np.where(data_shell == 0)] = np.nan

    # im = ax.imshow(data_shell,  cmap = 'jet', norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin='lower')
    # ax1.add_patch(data_reg)
    # print(data_shell)
    data_shell = data_shell[~np.isnan(data_shell)]
    # data_flux = np.sum(data_shell) / (np.pi * (r_out **2 - r_old **2))
    data_profile.append(np.mean(data_shell))
    # data_err.append(np.std(data_shell)/ np.sqrt(data_shell.size))
    
    data_old = data_mask 
    r_old = r_out
    # plt.show()

plt.show()


fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
print(r_list)
print(data_profile)
print(data_err)
ax.plot(r_list, data_profile, color='darkblue')
# ax.errorbar(r_list, data_profile, data_err, fmt='o',ecolor='r',color='b',elinewidth=2,capsize=4)
# ax.semilogy()
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x*1e9:.1f}"))

# ax.set_xlim(-0.2, 5)
y_major = np.arange(0, 0.015, 0.003)
y_minor = np.arange(0, 0.015, 0.001)
ax.tick_params(direction='in')
major_ticks = np.arange(0, 121, 10.)  # 从0到80，步长为20
ax.set_xticks(major_ticks)
minor_ticks = np.arange(0, 121, 2)
minor_ticks = np.setdiff1d(minor_ticks, major_ticks)
ax.set_xticks(minor_ticks, minor=True)
major_ticks = y_major  # 从0到80，步长为20
minor_ticks = y_minor 

# ax.set_yticks(major_ticks)
# ax.set_yticks(minor_ticks, minor=True)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
# ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.1f}'.format(val)))
ax.set_ylim(0, 5e-9)

ax.set_xlabel(f'R ($arcsec$)')
ax.set_ylabel(r'Mean Flux Density ($mJy beam^{-1} cm^{-2}$)')
plt.show()

# %%
