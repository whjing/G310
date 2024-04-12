#%%


from re import T
from traceback import print_tb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
#from colourspace import maps
import copy
from astropy.modeling import models, fitting
from astropy.coordinates import SkyCoord
from regions import EllipseSkyRegion, EllipsePixelRegion
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
import matplotlib.ticker as ticker
from astropy.io import fits as pf
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from copy import deepcopy
import os,sys,re
import glob
import pandas as pd
from pandas import Series,DataFrame
from reproject import reproject_interp
from astropy import units as u




def spectral_map(a, x_hz, x, y):
    '''spectral map
    a: data narray
    x_hz: 
    '''
    
    spec_map_data = np.zeros((y,x))
    
    for i in range(y):
        for j in range(x):

            tmp_data = a[:,i,j]
            if  np.isnan(tmp_data).any() == False:
                fit_data = np.log10(tmp_data[tmp_data>0.])
                x_fit = np.log10(x_hz[tmp_data>0.])
                 
                model = models.Linear1D()
                fitter = fitting.LinearLSQFitter()
                best_fit = fitter(model, x_fit, fit_data)
                spec_map_data[i,j] = best_fit.parameters[0]
                
                #check the data fit
                if i==110 and j==79:
                    fig = plt.figure(figsize=(12,9))
                    ax = fig.add_subplot(111)
                    model = models.Linear1D()
                    fitter = fitting.LinearLSQFitter()
                    best_fit = fitter(model, x_fit, fit_data)
                    ax.plot(x_fit, fit_data,'k.',alpha=1 ,label='data',markersize=12)
                    ax.plot(x_fit, best_fit(x_fit), color='r', linewidth=1.5,label = 'fit')
                    print(best_fit.parameters[0])
                
            else:
                spec_map_data[i,j] = np.nan
    return spec_map_data



def main():
    #In put three dimension cub fits file",
    #path='/Users/wst/fwq/wsclean/10463/beam14/casa-soomth/'
    file_name = 'kkbr_cube_1GHz'
    file = path + file_name + '.fits'
    hdu = pf.open(file)
    hdr = pf.getheader(file)
    # x_hz = np.fromfile(path + 'freq.txt', sep='\n')
    # hdrx_hz = ('[1752519369.13, 2562910079.96]')
    x_hz = np.array(eval(hdr['FREQ']))
    print(x_hz)
    data = pf.getdata(file)
    axis_x = hdu[0].header['NAXIS1']
    axis_y = hdu[0].header['NAXIS2']
    spectral_map_data = spectral_map(data, x_hz, axis_x, axis_y)

    pf.writeto(path + 'spectral_index_2dimension.fits', spectral_map_data, hdr, overwrite=True)
    
    #draw map
    ##############################################################
    file1 = path + 'spectral_index_2dimension' + '.fits'
    hdu1 = pf.open(file1)
    hdr1 = pf.getheader(file1)
    wcs = WCS(hdr1)
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_subplot(111,projection=wcs,slices=('x','y'))
    vmax = 1
    vmin = -1
    im = ax.imshow(spectral_map_data, vmin=vmin, vmax=vmax , cmap = 'viridis',origin='lower')
    #设置colorbar的位置和大小
    norm = ImageNormalize(spectral_map_data, interval=MinMaxInterval(), stretch=SqrtStretch())
    cbar = plt.colorbar(im, ax=ax,pad=0.007)
    cbar.ax.set_yticklabels(cbar.ax.get_yticks(), fontsize=18)
    #cbar.set_label(r'mJy/beam',fontsize=22,labelpad=15)
    fmt = ticker.FormatStrFormatter('%.2f')  #颜色棒的值保留1位小数
    cbar.ax.yaxis.set_major_formatter(fmt)
    #plt.savefig(path + file_name + '.png', dpi = 300)
    plt.show()
    ##############################################################

if __name__ == "__main__":
    main()











