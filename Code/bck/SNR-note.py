# coding:utf-8  #定义coding则告诉系统脚本是何编码格式
#！/usr/bin/env #定义#！,会去找指定路径下的python解释器
#%%
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from matplotlib.patches import Circle
import numpy as np
from astropy.io import fits as pf
from astropy.wcs import WCS
from matplotlib import pylab as plt
from matplotlib.colors import LogNorm
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord
from matplotlib.patches import Ellipse, Rectangle
#from regions import EllipseSkyRegion, EllipsePixelRegion
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
from glob import glob


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


def sv_fig(name):
    sv = plt.savefig(name + '.png', bbox_inches='tight', dpi = 1000)
    sv = plt.savefig(name + '.pdf', bbox_inches='tight', dpi = 1000)
    return sv


def show_cb(im, vmin, vmax, stp = 7, n = 1000, label = r'mJy beam$^{-1}$'):
    cb = plt.colorbar(im,  pad=0.1, fraction=0.04, location='right')#orientation='horizontal')fraction=0.05,
    colors_range = np.linspace(vmin, vmax, stp)
    cb.set_ticks(colors_range)
    tick_label = []
    for l in range(stp):
        tick_label.append(round(colors_range[l] * n, 3))
    cb.set_ticklabels(tick_label)
    cb.ax.tick_params(labelsize=12)
    cb.set_label(label)
    return cb
#########################################

args = sys.argv[1:]
fname = args[0]
fname = '/Users/jing/EMU_dataset/image.i.EMU_0208-09B.SB59253.cont.taylor.0.restored.conv.fits'

hdr, data, wcs = load_fits_image(fname)
print('load')

for p1,p2 in [(10, 99), (5, 90)]:
    print(p1,p2)

    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': wcs}, figsize=(12,12))

    overlay = ax.get_coords_overlay('galactic')
    overlay.grid(color='yellow', ls='dotted', alpha = 0.5)
    overlay[0].set_axislabel('Galactic Longitude')
    overlay[1].set_axislabel('Galactic Latitude')


    # SNRcat
    df_SNR = pd.read_csv('~/academic/Cat/SNRcat.csv', sep=';', skiprows=2)

    for index, row in df_SNR.iterrows():
        coord = SkyCoord(row['J2000_ra (hh:mm:ss)'], row['J2000_dec (dd:mm:ss)'], unit=(u.hourangle, u.deg), frame='icrs')
        x, y = wcs.wcs_world2pix(coord.ra.deg, coord.dec.deg, 0)
        r = 0.5 * row['size_coarse (arcmin)']/((60*(abs(hdr['CDELT2']))))
        if row['J2000_from'] == 'Green J coord':
            color = 'lightblue'
        else:
            color = 'limegreen'
        circle = Circle((x,y), radius = r, facecolor='none', edgecolor=color, alpha=1, zorder=200)
        ax.add_patch(circle)
        ax.annotate(row['G'], (x,y), xytext=(r/100,r/100), textcoords='offset points',color= color)



    df_HII = pd.read_csv('/Users/jing/academic/Cat/wise_hii_V2.2.csv', sep=',')#, skiprows=2)

    # df_HII = df_HII[(df_HII['Catalog'] == 'K') | (df_HII['Catalog'] == 'C') | (df_HII['Catalog'] == 'G')]

    for index, row in df_HII.iterrows():
        coord = SkyCoord(row['GLong<br>(deg.)'], row['GLat<br>(deg.)'], unit=(u.deg, u.deg), frame='galactic').transform_to('icrs')
        # print(coord)
        x, y = wcs.wcs_world2pix(coord.ra.deg, coord.dec.deg, 0)
        # print(x,y)
        r = 0.5 * row['Radius<br>(arcsec.)']/((3600*(abs(hdr['CDELT2']))))
        if row['Catalog'] == 'K':
            color = 'red'
        elif row['Catalog'] == 'C':
            color = 'pink'
        elif row['Catalog'] == 'Q':
            color = 'magenta'
        else:
            color = 'grey'
        circle = Circle((x,y), radius = r, facecolor='none', edgecolor=color, alpha=1, zorder=200)
        ax.add_patch(circle)
        # ax.annotate(row['WISE Name'], (x,y), xytext=(r/100,r/100), textcoords='offset points',color= color)



    df_ATNF = pd.read_csv('~/academic/Cat/ATNF.txt', sep=';', usecols=(2,3,5), na_values='NAN')
    df_ATNF_Young = df_ATNF#df_ATNF[df_ATNF['AGE'] < 2e4]
    ATNF = SkyCoord(df_ATNF_Young['RAJ'], df_ATNF_Young['DECJ'], unit=(u.hourangle, u.deg), frame='icrs')
    ax.scatter(ATNF.ra.deg, ATNF.dec.deg, marker='+', color='cyan', transform=ax.get_transform('world'), label='ATNF')

    
    vmin = np.nanpercentile(data,p1)
    vmax = np.nanpercentile(data,p2)
    im = ax.imshow(data,vmin=vmin,vmax=vmax, interpolation='spline36', cmap='jet', origin='lower')
    ax.coords.grid(color = 'grey', linestyle = 'dotted', alpha = 0.5)

    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    plt.tight_layout()
    show_cb(im, vmin, vmax,)
    plt.title(fname.split('/')[-1].split('.')[2:4], y = -0.1)

    outname = fname.replace('.fits', '')+ '_' + str(p1) + '-' + str(p2)#.replace('RawData', 'Output')
    sv_fig(outname) 
    print('saved')


# %%
