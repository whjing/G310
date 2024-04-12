#!/usr/bin/env python
# %%

from astropy.wcs.wcs import WCS
# from matplotlib.pyplot import plot, sca
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
from globalSetting import *

set_plot_settings()

#-------   general ----------#
# plt.rcParams.update({'font.size': 18})
def set_ax(ax, num_ticks=5, minpad = 1, declabel = "default" ):
    ra = ax.coords[0]
    dec = ax.coords[1]

    ra.set_axislabel('Right Ascension (J2000)',minpad = minpad)

    
    dec.set_axislabel('Declination (J2000)', minpad = minpad,)
    if declabel == None:
        dec.set_axislabel(' ', minpad = minpad,)
        
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
    

    ax.set_xlim(40,90)
    ax.set_ylim(37,87)

    ax.patch.set_facecolor("aliceblue") 


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



def read_fits(file_path):
    with pf.open(file_path) as hdul:
        header = hdul[0].header
        data = hdul[0].data 
        wcs = WCS(header)
    return header, data, wcs

def main():

    file_I = '../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-250.fits'
    file_PI = "../Fits_smallCube/test_ampPeakPIfitEff_crop-250.fits"
    file_PA = "../Fits_smallCube/test_polAngle0Fit_deg_crop-250.fits"
    file_RM = "../Fits_smallCube/test_phiPeakPIfit_rm2_crop-250.fits"

    # read fits
    hdr_I , data_I,  wcs_I  = read_fits(file_I)
    hdr_PI, data_PI, wcs_PI = read_fits(file_PI)
    hdr_PA, data_PA, wcs_PA = read_fits(file_PA)
    hdr_RM, data_RM, wcs_RM = read_fits(file_RM)

    # parameters
    box_s = f"fk5;polygon(210.1910460,-63.4149374,210.1602719,-63.4234925,210.1440260,-63.4398608,210.1785609,-63.4432171,210.2160066,-63.4340989,210.2334541,-63.4167919)"
    box = pyregion.parse(box_s)
    # hdu = fits.PrimaryHDU(data, header=header)
    # data_box = box.get_mask(hdu=hdu,shape = RM_data.shape)
    # print(data_box)

    maskcondition = (data_I < rms_I * 5) | (data_PI < 3 *rms_PI) #| (data_box == False)


    data_PI_img = data_PI.copy()

    data_PI_img[maskcondition] = np.nan
    data_PA[maskcondition] = np.nan
    # PI_data[np.where(Frac_data < 0.005)] = np.nan
    fig = plt.figure(figsize=(30,9))
    ax = plt.subplot(131, projection = wcs_PI)
    vmin =  np.nanpercentile(data_PI_img, 0)
    vmax =  np.nanpercentile(data_PI_img, 99)
    im = ax.imshow(data_PI_img, origin='lower',cmap='GnBu_r', vmin=vmin, vmax=vmax)#, alpha = 0.8)
    set_ax(ax)
    show_ctr(ax, data_I, hdr_I, hdr_PI, rms_I, num = 14, hpbw = 1.5, color = 'gray', linewidths = 1)
    # show_ctr(ax, PI_data, PI_hdr, PI_hdr, rms_PI, num = 6,hpbw = 1.5, color = 'red', linewidths = 1)

    scale = 0.005
    Nx = 2  #这两个值为偏振取值，越小对偏振的取值次数就会越频繁，反之亦然
    Ny = 2
    sigma = 3 * rms_PI #g315.PI.fits图像中的 rms=0.055 mJy，所以5rms=0.275 mJy
    M0, N0 = data_PA.shape
    count = 0

    for i in range(0,M0,Nx):
        for j in range(0,N0,Ny):
            if np.isnan(data_PI[i,j]) or data_PA[i,j] < sigma: continue    #没有偏振角或者偏振小于sigma，就跳出循环
            pa0 = np.radians(data_PA[i, j]) + np.pi/2.
            pi0 = np.sqrt(data_PI[i, j]* scale) / 2.
            l, b = wcs_PA.wcs_pix2world([[j+1,i+1]], 0)[0][0:2]

            y1 = b+pi0*np.cos(pa0)
            x1 = l+pi0*np.sin(pa0)/np.cos(np.radians(y1))
            y2 = b-pi0*np.cos(pa0)
            x2 = l-pi0*np.sin(pa0)/np.cos(np.radians(y2))

            px1, py1 = wcs_PA.wcs_world2pix([[x1, y1]], 0)[0][0:2]
            px2, py2 = wcs_PA.wcs_world2pix([[x2, y2]], 0)[0][0:2]

            plt.plot([px1,px2],[py1,py2],'tomato', linewidth = 4, alpha = 1)
            count += 1  #统计画的偏振条目数 

    # plot the PSR

    c_PSR=('14:00:45.69',' -63:25:42.6')
    c_PSR= SkyCoord(ra=c_PSR[0], dec=c_PSR[1], frame='fk5',unit=(u.hourangle, u.deg))
    ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=300, marker='+', color='cyan' ,transform=ax.get_transform('world'),linewidth = 5, zorder=200)
    print(c_PSR.ra.deg, c_PSR.dec.deg,)
    show_cb(im, vmin, vmax, n=1e6, label = r'$\rm{\mu Jy ~beam^{-1}}$')
    # sv_fig('../figure/PI')
                                                               

    ax = plt.subplot(132, projection = wcs_PI)

    Frac_data = data_PI/data_I

    Frac_data[maskcondition] = np.nan

    vmin, vmax = 0.01, 0.05
    im = ax.imshow(Frac_data, origin = 'lower',cmap = 'rainbow', vmin = vmin, vmax = vmax) 
    show_cb(im, vmin, vmax, n = 100, label = r'%')
    show_ctr(ax, data_I, hdr_I, hdr_PI, rms_I, num = 14, hpbw = 1.5, color = 'gray', linewidths = 1)
    # show_ctr(ax, data_PI, hdr_PI, hdr_PI, rms_PI, num = 6, hpbw = 1.5, color = 'tomato', linewidths = 2)
    ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=300, marker='+', color='cyan' ,transform=ax.get_transform('world'), linewidth = 5)
    set_ax(ax, declabel=None)


    #----------------#
    # RM 

    data_RM[maskcondition] = np.nan

    ax = plt.subplot(133, projection = wcs_RM)

    vmin = -750
    vmax = -600

    im = ax.imshow(data_RM, origin = 'lower',cmap='Spectral', vmin = vmin, vmax = vmax)
    show_cb(im, vmin, vmax, n=1, label = r"rad m$^{-2}$" )
    show_ctr(ax, data_I, hdr_I, hdr_PI, rms_I, num = 14, hpbw = 1.5, color = 'gray', linewidths = 1)
    ax.scatter(c_PSR.ra.deg, c_PSR.dec.deg, s=300, marker='+', color='cyan' ,transform=ax.get_transform('world'), linewidth = 5)
    set_ax(ax, declabel=None)
    
    # sv_fig("../figure/PI-RM")

    sv_fig('../figures/PI_all')
    

if __name__ == "__main__":
    main()

# %%
