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
# from matplotlib.patches import Annulus
import pyregion
from sklearn.preprocessing import StandardScaler
from matplotlib.ticker import ScalarFormatter, FuncFormatter
# sys.path.append("~/academic/module/")



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
    ymin = np.nanmin(data)
    ymax = np.nanmax(data)
    dy = ymax - ymin 
    normalized_data = [(y_i - ymin) / (ymax - ymin) for y_i in data]
    norm_data_err = [(y_i_err) / (ymax - ymin) for y_i_err in data_err]
    return normalized_data, norm_data_err, dy


##################################
###### EMU data 
##################################


# open fits file 


def plot_profile(file, ra, dec, vmin, colorbar_times, cmap = "jet"):
    with pf.open(file) as hdul:
        print(hdul)
        hdu = hdul[0]
        hdr = hdu.header
        data = hdu.data
        wcs = WCS(hdr)

    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(111, projection = wcs)

    # Mask Data

    box_s = f'fk5;polygon(210.2649256,-63.4323334,210.1940403,-63.4267705,210.2184661,-63.3950644,210.1172636,-63.3946010,210.1156781,-63.4285364,210.1800642,-63.4303422,210.1455727,-63.4613691,210.2669987,-63.4613567,210.2659964,-63.4602409)'
    # box_s = f'fk5;polygon(14:00:54.0800,-63:23:09.825,14:00:14.5310,-63:24:46.807,14:00:45.6449,-63:25:45.558)'
    box = pyregion.parse(box_s)
    data_box = box.get_mask(hdu=hdu,shape=data.shape)
    data = data_box * data
    # data = np.where(data < vmin, vmin, data)
    # print(f"data min = {np.nanmin(data)}")
    vmax = np.nanpercentile(data, 98)
    
    im = ax1.imshow(data,  cmap = cmap, vmin = vmin, vmax = vmax, origin='lower') #vmin = vmin, vmax = vmax)norm=colors.LogNorm(
    set_ax(ax1)
    # show_cb(im, vmin = vmin, vmax = vmax, n = colorbar_times)


    # ra = "210.1882728"
    # dec = "-63.4280232"
    r = 120

    shell_c = SkyCoord(ra, dec, frame='fk5',unit=(u.deg, u.deg))
    x, y = wcs.wcs_world2pix(shell_c.ra.deg, shell_c.dec.deg, 0)
    r_list = np.linspace(0, r, 21)
    # r_list_ring = np.linspace(0, r_out, 30)

    data_old = np.zeros_like(data)
    data_profile = []
    data_err = []
    r_old = 0

    for radius in r_list:

        reg =f'fk5;circle({ra},{dec},{radius}")'
        outcircle = pyregion.parse(reg)
        

        r_out = radius/((3600*(abs(hdr['CDELT1']))))
        color = 'tomato'

        if (radius > 60) & (radius<80):
            color = 'magenta'

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
        data_flux = np.mean(data_shell)
        # np.sum(data_shell) / (np.pi * (r_out **2 - r_old **2))
        data_profile.append(data_flux)
        data_err.append(5 * np.std(data_shell)/ np.sqrt(data_shell.size))
        
        data_old = data_mask 
        r_old = r_out

    sv_fig(f"../figures/{ra_r}-{dec_r}_0_{file.split('/')[-1]}_image")

    plt.show()

    data_orig = np.array(list(zip(data_profile, data_err)))
    data_profile, data_err, data_max = normalize_and_scale(data_profile, data_err)
    data_norm = np.array(list(zip(data_profile, data_err)))
    print(f"data_orig shape {data_orig.shape}\ndata_norm shape {data_norm.shape}")

    data_all = np.hstack((data_orig, data_norm))
    print(f"data_all shape {data_all.shape}")
    
    return r_list, data_profile, data_err, data_max


def plot_2Dprofile(ra_r, dec_r):

    file_askap = '../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-250.fits'

    # file_chandra_smooth = '../Fits_4-subband-useful/hard_flux.img.smooth_bin-2.fits'
    file_chandra_smooth = '../Fits_4-subband-useful/hard_flux.img.smooth_bin-0.5.fits'
    RMS_ASKAP = 3e-05
    RMS_Chandra = 2.05395e-09
    RMS_Chandra_sm = 4.20539e-10

    r_list, data_profile_askap, data_err_askap, data_max_askap = plot_profile(file_askap, ra = ra_r, dec = dec_r, vmin = RMS_ASKAP, colorbar_times = 6 ,cmap = 'CMRmap_r')

    # r_list, data_profile_chandra, data_err_chandra, data_max_chandra = plot_profile(file_chandra, vmin = RMS_Chandra, colorbar_times = 12, cmap = "ocean_r") 
    # Chandra 的流量还要再确定一下
    r_list, data_profile_chandra_smooth, data_err_chandra_smooth, data_max_chandra_sm = plot_profile(file_chandra_smooth, ra = ra_r, dec = dec_r,vmin = RMS_Chandra_sm, colorbar_times = 12, cmap = "ocean_r") 


    def plot_profile_2D(ax):
        # ax.plot(r_list, data_profile_askap, color='darkred', label = "ASKAP")
        # ax.plot(r_list, data_profile_chandra, color='darkroyalblue', label = "Chandra")
        # 应该是6还是12
        ax.errorbar(r_list, data_profile_askap, data_err_askap,  fmt='o',ecolor='red',color='red',elinewidth=2,capsize=4, label = "ASKAP", alpha = 0.5) #xerr = 6,
        # ax.errorbar(r_list, data_profile_chandra, data_err_chandra, fmt='o',ecolor='r',color='darkroyalblue',elinewidth=2,capsize=4, label = "Chandra", alpha = 0.8)
        ax.errorbar(r_list, data_profile_chandra_smooth, data_err_chandra_smooth, fmt='o',ecolor='royalblue',color='royalblue',elinewidth=2,capsize=4, label = "Chandra convolved", alpha = 0.5)

        # ax.semilogy()
        ax.yaxis.set_major_formatter(ScalarFormatter())
        # ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))


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
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.1f}'.format(val)))
        # ax.set_ylim(0,)
        for radius in r_list:
            ax.axvline(x=radius, color='grey', linestyle='--', linewidth=1, alpha = 0.5)
        # ax.plot(line)

        ax.set_xlabel(r'R ($\rm{arcsec}$)')
        ax.set_ylabel(r'Normalized Flux Density (A.U)')#($mJy beam^{-1} cm^{-2}$)')


    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111)
    plot_profile_2D(ax)
    plt.legend(loc='lower right',bbox_to_anchor=[1,0.87], ncol = 2, shadow = False, fancybox= True) #loc="upper right",, 

    ax.set_xlim(-2, 125)
    ax_in = fig.add_axes([0.35, 0.4, 0.55, 0.38])
    plot_profile_2D(ax_in)
    # ax_in.set_xlim(43, 125)
    ax_in.set_xlim(35, 125)
    ax_in.set_ylim(0, 0.03)
    ax_in.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x*1000:.1f}"))
    ax_in.set_ylabel(r"yaxis $\times$ 1000")
    # ax_in.axhline(y  = RMS_ASKAP / data_max_askap, color = "red", linestyle="--")
    # ax_in.axhline(y  = RMS_Chandra / data_max_chandra, color = "darkroyalblue")
    # ax_in.axhline(y  = RMS_Chandra_sm / data_max_chandra_sm, color = "royalblue", linestyle="--")
    ax.set_title(f"{ra_r}  {dec_r}")

    sv_fig(f"../figures/{ra_r}-{dec_r}_profile")
    plt.show()

import numpy as np
import random
ras = 210.1882728
decs = -63.4280232
range = 0.004
raList = np.random.uniform(ras-range/2, ras+range/2,20)
decList = np.random.uniform(decs-range/2, decs+range/2,20)
# Combine raList and decList into a list of tuples

coordinates = list(zip(raList, decList))

# Shuffle the list of coordinates
random.shuffle(coordinates)

# Select 10 unique coordinates
unique_coordinates = random.sample(coordinates, 10)

# Output the unique coordinates
for coord in unique_coordinates:
    ra_r = str(coord[0])
    dec_r = str(coord[1])
    plot_2Dprofile(ra_r, dec_r)

# %%
