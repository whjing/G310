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
from astropy import units as u
from globalSetting import *
plt.rcParams.update({'font.size': 22})
plt.rcParams['lines.linewidth'] = 3





"""
Usage: input File has to be the same data shape
"""


def smooth(data_in, hpbw=4.):
    """
    smoothing
    """
    g = Gaussian2DKernel(hpbw)

    data_out = convolve(data_in, g, boundary='extend')

    return data_out


grid = plt.GridSpec(2, 4, wspace=0.4, hspace=0.2)




def snr_plot(name,file,position,vmin=20,vmax=90,cmap='viridis',title='unknown'):
    
    hdu = pf.open(name)[0]
    hdr = hdu.header
    wcs = WCS(hdr)
    data = hdu.data
    
    data[np.isnan(data)] = 0
    vmin = np.percentile(data,vmin)
    vmax = np.percentile(data,vmax)
    
    ax = plt.subplot(position,projection = wcs)
    im = ax.imshow(data, vmin=vmin, vmax=vmax, interpolation='spline36', cmap=cmap, origin='lower')


    hdu_x = pf.open(file)[0]
    data_x, footprint = reproject_interp(hdu_x, hdr)
    # data_x = smooth(data_x,1.1)
    ax.contour(data_x, levels=0, colors='magenta')
    
    

    ra = ax.coords[0]
    dec = ax.coords[1]
    
    ra.set_ticks_visible(False)
    ra.set_ticklabel_visible(False)
    dec.set_ticks_visible(False)
    dec.set_ticklabel_visible(False)


    ax.set_title(title,fontsize=18,y=-0.15)
    
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    
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


def ttplot(file,file1,file2,f1_freq,f2_freq,position):

    #设置频率
      # in GHz
    freq1 = f1_freq
    freq2 = f2_freq
	#对指定区域进行拟合，例如本文就是设定mask大于0.5区域是目标源
    mask = pf.getdata(file)  
    print(f'mask shape {mask.shape}')


	#导入两个频率下的温度数据
    data1 = pf.getdata(file1)[mask>0.5]
    data2 = pf.getdata(file2)[mask>0.5]
    print(f'before select {data1.min(), data2.min()}')

    mask_condition = (data1 > 3 * rms_I) & (data2 > 3 * rms_X)
    # mask_condition = (data1.squeeze() > 3 * rms_I) & (data2.squeeze() > 3 * rms_X)
    data1 = data1[mask_condition][::5]
    data2 = data2[mask_condition][::5]
    print(f'{data1.shape}, {data2.shape}')   
    print(f'after select {data1.min(), data2.min()}') 
    data1 = data1 * 1e3

    ax = plt.subplot(position)

    import numpy as np

    # 进行线性拟合
    coefficients = np.polyfit(data1, data2, 1)
    a, b = coefficients
    
    # 构造包含拟合公式的字符串
    fit_equation = "y = {:.2e}x + {:.2e}".format(a, b)

    
    # 绘制线性拟合直线，并使用label显示拟合公式
    ax.plot(data1, a * data1 + b, label=fit_equation)
    
    # 添加图例以显示拟合公式
    ax.legend()

    ax.plot(data1, data2, 'o')

    ax.set_xlabel('mJy/beam (at %5.2f MHz)' % freq1.value)
    ax.set_ylabel('flux/pixel (at %1.0f keV)' % freq2.value)

    ax.set_title("radio - x-ray")
    return ax
# aim = ['bottom','left','right','all']

# for i in range(len(aim)):




freq1 = 943*u.MHz 
freq2 = 2.5*u.keV  
maskList = ['SNR', 'shell', 'shell-1', 'PWN', 'shell-Chandra', 'PWN-Chandra', 'shell-p1', 'shell-p2']

for maskFile in maskList:
    fig = plt.figure(figsize=(16, 7),dpi=300)
    file = (f'../Fits_4-subband-useful/{maskFile}.mask_crop-250.fits')
    print(file)

    file1 = ('../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-250.fits')
    position = grid[0, 0]
    snr_plot(file1,file,position,vmin=40, vmax=90,title='ASKAP')


    
    file2 = ('../Fits/G310_i_Chandra_r250_smooth.fits')

    position = grid[1, 0]
    snr_plot(file2,file,position,vmin=20, vmax=95,cmap='hot',title='Chandra')

    
    position = grid[0:, 1:]
    ttplot(file,file1,file2,freq1,freq2,position)
    
    # plt.suptitle(i)
        
        # plt.show()为什么加上这个就保存不出来图片了
    output = (f'../Output/relation_{maskFile}.pdf')
    plt.savefig(output, bbox_inches='tight')
        
    plt.show()
# %%
