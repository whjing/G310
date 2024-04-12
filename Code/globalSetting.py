import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pf
from astropy.wcs import WCS
from matplotlib import pylab as plt
from matplotlib.colors import LogNorm
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord
from matplotlib.patches import Ellipse, Rectangle
from regions import EllipseSkyRegion, EllipsePixelRegion
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
import matplotlib.ticker as ticker
from copy import deepcopy
from astropy import units as u
import os,sys,re

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from reproject import reproject_interp
from astropy import units as u
from astropy import constants as const
# ++++++++ Parameters +++++++++
rms_I  = 2.5e-5
rms_PI = 5.3e-5
rms_X =  4.356845238665e-10
print(f"+++ Here are global parameters +++")
print(f"rms_I = {rms_I}")
print(f"rms_PI = {rms_I}")
print(f"rms_X = {rms_X}")
print(f"++++++++++++++++++++++++++++++++++")

# set the plot format
def set_plot_settings():
    config = {
        "font.family": 'serif',
        "font.size": 22,
        "mathtext.fontset": 'stix',
        "font.serif": ['SimSun'],  # simsun字体中文版就是宋体
    }
    plt.rcParams.update(config)
    
    plt.rcParams['font.size'] = 22
    
    plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
    plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内
    
    plt.rc('font', family='Times New Roman')#!/usr/bin/env python3





def load_fits_image(filename):

    # Read the header and image data from the file
    header=pf.getheader(filename)
    data = pf.getdata(filename)
    naxis = len(data.shape)

    # Strip unused dimensions from the data array
    if naxis == 2:
        xydata = data.copy()
        del data
    elif naxis == 3:
        xydata = data[0].copy()
        del data
    elif naxis == 4:
        xydata = data[0][0].copy()
        del data
    elif naxis == 5:
        xydata = data[0][0][0].copy()
        del data
    else:
        print("Data array contains %s axes" % naxis)
        print("This script supports up to 5 axes only.")
        sys.exit(1)

    # Strip unused dimensions from the header
    stripHeadKeys = ['NAXIS','CRVAL','CRPIX','CDELT','CTYPE','CROTA',
                     'CD1_','CD2_','CUNIT','PC1_','PC2_','PC3_','PC4_']
    
    stripHeadKeys1 = ['PC1_','PC2_','PC3_','PC4_','PC01_','PC02_','PC03_','PC04_']
    if naxis >= 2:
        for i in range(3,6):
            for key in stripHeadKeys:
                if (key+str(i)) in header: del header[key+str(i)]
        for j in range(1,6):            
            for key1 in stripHeadKeys1:
                if (key1+str(j)) in header: del header[key1+str(j)]
                if (key1+ '0' + str(j)) in header: del header[key1+ '0' + str(j)]
        header['NAXIS'] = 2

    # Determine the coordinate type and the corresponding keyword index
    # Make a note in the header
    ra_regx = re.compile('^RA')
    dec_regx = re.compile('^DEC')
    glon_regx = re.compile('^GLON')
    glat_regx = re.compile('^GLAT')
    if 'CTYPE1' in header:
        if ra_regx.match(header['CTYPE1']) or glon_regx.match(header['CTYPE1']):
            for i in range(int(header['NAXIS'])):
                keyword = "CTYPE"+str(i+1)
                if ra_regx.match(header[keyword]): coord_type="EQU"; x_indx=i+1
                if dec_regx.match(header[keyword]): y_indx=i+1
                if glon_regx.match(header[keyword]): coord_type="GAL"; x_indx=i+1
                if glat_regx.match(header[keyword]): y_indx=i+1
            if not x_indx or not y_indx:
                print("Failed to find Equatorial or Galactic axis coordinate types.")
                del data; del header
                sys.exit(1)
            else:
                header['XINDEX'] = x_indx
                header['YINDEX'] = y_indx
        else:
            print('Does not have "RA-DEC" and "GLON-GLAT" coordinates!!!')
    else:
        print('key "CTYPE1" does not, please check the header!!!')
    # Convert AIPS clean-beam types to standard BMAJ, BMIN, BPA
    try:
        bmaj = header['CLEANBMJ']
        bmin = header['CLEANBMN']
        bpa = header['CLEANBPA']
        header.update('BMAJ',bmaj)
        header.update('BMIN',bmin)
        header.update('BPA',bpa)
    except Exception:
        print("No AIPS style beam keywords found.")

    # Check for PIXSCAL keywords and write to CDELT standard
    try:
        xdelt=(-1)*(header['PIXSCAL'+str(x_indx)])/3600.0
        ydelt=(header['PIXSCAL'+str(y_indx)])/3600.0
        header['CDELT'+str(x_indx)] = xdelt
        header['CDELT'+str(y_indx)] = ydelt
    except Exception:
        pass


    wcs = WCS(header)

    return [header,xydata,wcs]


def show_img(data, wcs, cb = True ):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = wcs)
    vmin, vmax = vmin, vmax
    im = ax.imshow(data, origin='lower', cmap = 'viridis', vmin=vmin, vmax=vmax)  
    if cb == True:
        show_cb(im, vmin, vmax) 
    set_ax(ax)   
    
    return ax


def smooth(data_in, hpbw=2.):
    g = Gaussian2DKernel(hpbw)
    data_out = convolve(data_in, g, boundary='extend')
    return data_out


def show_ctr(ax, data_orig, hdr_ctr, hdr_img, rms, num = 6, hpbw = 2., color = 'tomato', linewidths = 2):
    data_s = smooth(data_orig, hpbw)
    levels=[]
    for i in range(num): 
        level = rms* 3 * (2** (i/1.5))
        levels.append(level)
    print(levels)
    data_ctr, _ = reproject_interp(pf.PrimaryHDU(data_s, hdr_ctr), hdr_img)
    contour = ax.contour(data_ctr, levels=levels, colors = color, linewidths = linewidths)
    return contour



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

def set_ax(ax, num_ticks=5, minpad = 1 ):
    ra = ax.coords[0]
    dec = ax.coords[1]
    ra.set_axislabel('Right Ascension (J2000)',minpad = minpad)
    dec.set_axislabel('Declination (J2000)', minpad = minpad,)
    
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_ticks(number = num_ticks)
    dec.set_ticks(number = num_ticks)
    dec.set_ticklabel_visible(True)
    #dec.set_ticklabel(rotation='vertical', va='center')



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


def sv_fig(name):
    sv = plt.savefig(name + '.png', bbox_inches='tight', dpi = 300)
    sv = plt.savefig(name + '.pdf', bbox_inches='tight')
    return sv

def cpt_EF(E, freq):
    h = const.h
    if E == True:
        freq = (E/h).to(u.Hz)
        return freq
    elif freq == True:
        E = (h * freq).to(u.keV)
        return E
        