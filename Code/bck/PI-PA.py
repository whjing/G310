#%%
# Todo list:
# 1. Just display the center part.



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

import sys
import pyregion
import os, sys

# Set my module
# os.chdir(sys.path[0])
from Code.globalSetting import *
set_plot_settings()

#-------   general ----------#

# plt.rcParams['font.size'] = 22
# plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
# plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内

# plt.rc('font',family='Times New Roman')

#需要对PI进行比例的缩放，具体缩放程度看PI和I的数值差距,没看懂

I_hdr,I_data, I_w = load_fits_image('../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-250.fits')
PI_hdr, PI_data, PI_w = load_fits_image('../Fits_smallCube/test_ampPeakPIfitEff_crop-250.fits')
PA_hdr, PA_data, PA_w = load_fits_image('../Fits_smallCube/test_polAngle0Fit_deg_crop-250.fits')

# intensity 

I_data = I_data
PI_data = PI_data
PA_data = PA_data
fig = plt.figure()
ax = plt.subplot(111, projection=I_w)

vmin = np.nanpercentile(I_data,10)
#np.percentile(I_data,20)
vmax = np.nanpercentile(I_data,99.3)

im = ax.imshow(I_data, origin = 'lower', cmap='viridis', vmin = vmin, vmax = vmax)
show_cb(im, vmin, vmax)

set_ax(ax)
show_ctr(ax, I_data, I_hdr, PI_hdr, 0.01, num = 8,hpbw = 1.5, color = 'limegreen', linewidths = 1)
# sv_fig('kkbr_i')
# polarization
# ax.set_xlim(30,210)
# ax.set_ylim(30,210)

# fig = plt.figure(figsize=(10,10))

Frac_data = PI_data/I_data
PI_data[np.where(Frac_data > 0.1)] = np.nan
PI_data[np.where(Frac_data < 0.005)] = np.nan
fig = plt.figure()
ax = plt.subplot(111, projection=PI_w)
print(f"max {np.nanmax(PI_data)}")
vmin = np.nanpercentile(PI_data,10)
vmax = np.nanpercentile(PI_data,99.3)
im = ax.imshow(PI_data, origin='lower',cmap='spring_r', vmin=vmin, vmax=vmax) #norm=colors.LogNorm(vmin=vmin, vmax=vmax))cubehelix_r
ax.patch.set_facecolor("mintcream")
set_ax(ax)
show_ctr(ax, I_data, I_hdr, PI_hdr, 2.5e-5, num = 8, hpbw = 1.5, color = 'limegreen', linewidths = 1)


scale = 0.05
Nx = 4  #这两个值为偏振取值，越小对偏振的取值次数就会越频繁，反之亦然
Ny = 4
sigma = 0.00015 #g315.PI.fits图像中的 rms=0.055 mJy，所以5rms=0.275 mJy
M0, N0 = PA_data.shape
count = 0

for i in range(0,M0,Nx):
    for j in range(0,N0,Ny):
        if np.isnan(PA_data[i,j]) or PI_data[i,j] < sigma: continue    #没有偏振角或者偏振小于sigma，就跳出循环
        pa0 = np.radians(PA_data[i, j]) + np.pi/2.
        pi0 = np.sqrt(PI_data[i, j]* scale) / 2.
        l, b = PA_w.wcs_pix2world([[j+1,i+1]], 0)[0][0:2]

        y1 = b+pi0*np.cos(pa0)
        x1 = l+pi0*np.sin(pa0)/np.cos(np.radians(y1))
        y2 = b-pi0*np.cos(pa0)
        x2 = l-pi0*np.sin(pa0)/np.cos(np.radians(y2))

        px1, py1 = I_w.wcs_world2pix([[x1, y1]], 0)[0][0:2]
        px2, py2 = I_w.wcs_world2pix([[x2, y2]], 0)[0][0:2]

        plt.plot([px1,px2],[py1,py2],'cyan', linewidth = 1)
        count += 1  #统计画的偏振条目数 

show_cb(im, vmin, vmax, n=1000, label=r'mJy beam$^(-1)$')
sv_fig("PI")

# ax.set_xlim(30,210)
# ax.set_ylim(30,210)
fig = plt.figure()
ax = plt.subplot(111, projection=PI_w)

Frac_data = PI_data/I_data
Frac_data[np.where(Frac_data > 0.1)] = np.nan
vmin, vmax = 0.01, 0.1
im = ax.imshow(Frac_data, origin = 'lower',cmap = 'magma_r', vmin = vmin, vmax = vmax) 
show_cb(im, vmin, vmax, n = 100, label = r'%')
show_ctr(ax, I_data, I_hdr, PI_hdr, 0.01, num = 8, hpbw = 1.5, color = 'limegreen', linewidths = 1)
set_ax(ax)


# ax.set_xlim(30,210)
# ax.set_ylim(30,210)

# %%
