#%%

import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.io import fits as pf
from astropy.wcs import WCS
import numpy as np
from astropy.convolution import convolve,Gaussian2DKernel
from reproject import reproject_interp

import os
import sys
from scipy import odr
from scipy import signal
from astropy.modeling import models, fitting


grid = plt.GridSpec(2, 3, wspace=0.6, hspace=0.2)
# fig = plt.figure(figsize=(10, 6),dpi=72)

def snr_plot(name,file,position,title='unknown'):
    hdu = pf.open(name)[0]
    hdr = hdu.header
    wcs = WCS(hdr)
    data = hdu.data
    
    data[np.isnan(data)] = 0
    vmin = np.percentile(data,5)
    vmax = np.percentile(data,98.5)
    
    ax = plt.subplot(position,projection = wcs)
    im = ax.imshow(data, vmin=vmin, vmax=vmax, interpolation='spline36', cmap='viridis', origin='lower')
    
    hdu_x = pf.open(file)[0]
    data_x, footprint = reproject_interp(hdu_x, hdr)
    ax.contour(data_x, levels=0, colors='r')
    
    ax.set_title(title)
    
    # ax.set_xlabel('RA')
    # ax.set_ylabel('Dec')
    ra = ax.coords[0]
    dec = ax.coords[1]
    
    ra.set_ticks_visible(False)
    ra.set_ticklabel_visible(False)
    dec.set_ticks_visible(False)
    dec.set_ticklabel_visible(False)
    
    return im

#t-tplot

def cal_alpha(data1,data2):
    n = data1.size
    x = np.zeros(n) # 长度为n的数组，每个值都是0
    y = np.zeros(n)
    kk = 0
    for i in range(n):
        if np.isnan(data1[i]) or np.isnan(data2[i]): continue #如果有为nan的值，则跳过进行下一次迭代
        x[kk] = data1[i]
        y[kk] = data2[i]
        kk += 1

    if kk<10: return np.nan,np.nan,np.nan,np.nan

    p_init = models.Polynomial1D(1)
    fit_p = fitting.LinearLSQFitter()
    p = fit_p(p_init, x[:kk], y[:kk])
    a, b = p.c0.value, p.c1.value #c0 截距 c1 斜率
    return a, b

freq1 = 872/1.e3  # in GHz   
freq2 = 1016/1.e3  # in GHz
def ttplot(file,file1,file2,position):

    #设置频率
    
	
	#对指定区域进行拟合，例如本文就是设定mask大于0.5区域是目标源
    mask = pf.getdata(file)  

	#导入两个频率下的温度数据
    data1 = pf.getdata(file1)[mask>0.5]
    data2 = pf.getdata(file2)[mask>0.5]
	
	#以下不需要修改
    a1, b1 = cal_alpha(data1, data2)
    a2, b2 = cal_alpha(data2, data1)
    alpha1 = np.log10(b1)/np.log10(freq2/freq1)
    alpha2 = np.log10(b2)/np.log10(freq1/freq2)
    print(alpha1, alpha2)
    alpha_1 = (alpha1 + alpha2) / 2.
    d_alpha_1 = abs(alpha2 - alpha1)
    a1_p = a1 * 1.
    b1_p = b1 * 1.
    a2_p = -a2/b2
    b2_p = 1./b2
	
	#画图
    ax = fig.add_subplot(position)
    ax.plot(data1, data2, 'o')
    ax.plot(data1, a1_p + b1_p * data1)
    ax.plot(data1, a2_p + b2_p * data1)
    ss = r'Spectral index = $-$' + '%4.2f' % abs(alpha_1) + r'$\pm$' + '%4.2f' % d_alpha_1
    ax.text(0.1, 0.9, ss, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    ax.set_xlabel('J/beam ( at %5.3f GHz)' % freq1)
    ax.set_ylabel('J/beam ( at %5.3f GHz)' % freq2)
    # ax.set_title("S-S plot")
    return ax



file = "../Fits/G310_i_EMU_0-144_PWN.mask.fits"
file1 = "../Fits/G310_i_EMU_0-144.fits"
file2 = "../Fits/G310_i_EMU_144-288.fits"

fig = plt.figure(figsize=(12, 6), dpi=300)

# file = ('../fits/mask_' + aim[i] + '.fits')
# file1 = ('../fits/SUMSS_pysmooth.fits')
position = grid[0, 0]
title = " %5.3f GHz" % freq1
snr_plot(file1,file,position,title)

# file2 = ('../fits/4850MHz_cut1.fits')
position = grid[1, 0]
title = " %5.3f GHz" % freq2
snr_plot(file2,file,position,title)

position = grid[0:, 1:]
ttplot(file,file1,file2,position)

# plt.show()为什么加上这个就保存不出来图片了
output = ('../Output/T-Tplot.pdf')
plt.savefig(output, bbox_inches='tight')
    
plt.show()

# %%
import myfits
# %%
