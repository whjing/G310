#%%  #####!/usr/bin/env python
from astropy.io import fits as pf
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.wcs as wcs
from astropy.coordinates import Angle
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from astropy.wcs import WCS
from matplotlib.patches import Circle
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator


def smooth(data_in, hpbw=4.):
    """
    smoothing
    """
    g = Gaussian2DKernel(hpbw)
    data_out = convolve(data_in, g, boundary='extend')
    return data_out


plt.rcParams.update({'font.size': 22})
plt.rcParams['lines.linewidth'] = 3


def snr_plot(name,position,vmin=20,vmax=90,cmap='viridis',title='unknown'):
    
    hdu = pf.open(name)[0]
    hdr = hdu.header
    wcs = WCS(hdr)
    data = hdu.data
    
    data[np.isnan(data)] = 0
    vmin = np.percentile(data,vmin)
    vmax = np.percentile(data,vmax)
    
    ax = plt.subplot(position,projection = wcs)
    im = ax.imshow(data,  cmap=cmap, origin='lower',vmin=vmin, vmax=vmax, )


    ra_vals = [c0.ra.deg, c1.ra.deg, c2.ra.deg]
    dec_vals = [c0.dec.deg, c1.dec.deg, c2.dec.deg]
    
    # 绘制两条线
    plt.plot([ra_vals[0], ra_vals[1]], [dec_vals[0], dec_vals[1]], marker='o', color='pink', transform=ax.get_transform('world'))
    plt.plot([ra_vals[0], ra_vals[2]], [dec_vals[0], dec_vals[2]], marker='o', color='orange', transform=ax.get_transform('world'))
    

    ra = ax.coords[0]
    dec = ax.coords[1]
    
    ra.set_ticks_visible(False)
    ra.set_ticklabel_visible(False)
    dec.set_ticks_visible(False)
    dec.set_ticklabel_visible(False)


    ax.set_title(title,fontsize=18,y=-0.1)
    
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    

    return im


file1 = '../Fits/G310_i_EMU_r250.fits'
file2 = '../Fits/G310_i_Chandra_r250_smooth.fits'

with pf.open(file1) as hdul:
    hdr_askap = hdul['PRIMARY'].header
    data_askap = hdul['PRIMARY'].data
    wcs_askap = WCS(hdr_askap)

with pf.open(file2) as hdul:
    hdr_chdr = hdul['PRIMARY'].header
    data_chdr = hdul['PRIMARY'].data
    wcs_chdr = WCS(hdr_chdr)


pt0 = ('14:00:45', '-63:25:45.844')
pt1 = ('14:00:31.2569', '-63:25:45.782')
# pt1 = ('14:00:28.6276','-63:24:55.718')
# pt2 = ('14:00:44.8','-63:23:39.727')
# pt2 = ('14:00:52.8217', '-63:23:48.794')

pt2 = ('14:00:50.8221', '-63:24:06.983')
# pt2 = ('14:00:43.7745 -63:24:14.924')
c0 = SkyCoord(ra=pt0[0], dec=pt0[1], frame='fk5', unit=(u.hourangle, u.deg))
c1 = SkyCoord(ra=pt1[0], dec=pt1[1], frame='fk5', unit=(u.hourangle, u.deg))
c2 = SkyCoord(ra=pt2[0], dec=pt2[1], frame='fk5', unit=(u.hourangle, u.deg))

sep1 = c0.separation(c1).deg
sep2 = c0.separation(c2).deg
ang1 = c0.position_angle(c1).value
ang2 = c0.position_angle(c2).value + 2.*np.pi

print(np.degrees(sep1), np.degrees(sep2))
print(np.degrees(ang1), np.degrees(ang2))

d_sep = 0.001 # 0.5 arcmin
sep = max(sep1, sep2)
Npts = int(sep / d_sep) - 1
x = (np.arange(Npts)+0.5) * d_sep

d_ang = np.radians(1.)
Npf = int((ang2 - ang1) / d_ang)
ang_all = ang1 + np.arange(Npf)*d_ang

hdr_list=[hdr_askap,hdr_chdr]
data_list=[data_askap,data_chdr]

Num_inf = len(hdr_list)

yList = []
yminList = []
ymaxList = []
erryList = []


for j in range(Num_inf):

    prof_data = np.zeros((Npts, Npf))

    for i in range(Npf):
        hdr = hdr_list[j]
        data = data_list[j]
        w = wcs.WCS(hdr)

        cc = c0.directional_offset_by(ang_all[i], Angle(x*u.deg))
        
        N = cc.size
        radec = np.zeros((N,2))
        radec[:,0] = cc.ra.deg
        radec[:,1] = cc.dec.deg
        pix = w.wcs_world2pix(radec,0)
    

        y = np.zeros_like(x)
        for k in range(y.size):
            y[k] = data[int(pix[k][1]), int(pix[k][0])]

        prof_data[:,i] = y

    err_y = np.zeros_like(x)
    for k in range(y.size):
        y[k] = np.mean(prof_data[k,:])
        err_y[k] = np.std(prof_data[k,:]) / np.sqrt(err_y.size)

    yList.append(y)
    erryList.append(err_y)


fig = plt.figure(figsize=(24, 9))

grid = plt.GridSpec(2, 4, wspace=0.1, hspace=0.1)

def process_data(index,x=x, yList=yList, erryList=erryList, ):
    xp = x[6:] * 3600
    yp = yList[index][6:]
    eyp = erryList[index][6:]
    return xp, yp, eyp


xp1, yp1, eyp1 = process_data(0)
xp2, yp2, eyp2 = process_data(1) 

print(yp2)

def normalization(data, err):
    _range = np.max(data) - np.min(data)/10
    data_out =((data - np.min(data)) / _range)
    err = err / _range
    
    return data_out, err

yp1_p, eyp1 = normalization(yp1, eyp1)
yp2_p, eyp2 = normalization(yp2, eyp2)
# yp2_p, eyp2 = normalization(yp2, eyp2)[0] * 4, normalization(yp2, eyp2)[1] * 4
print
# yp2_p, _range = normalization(yp2)
print(np.min(yp1))

annu =sep*3600
r1_shell = 50
r2_shell = 95
color_cir = 'magenta'
line_cir='--'



position = grid[0, 0]
snr_plot(file1,position,vmin=40, vmax=90, title='ASKAP')


position = grid[1, 0]
snr_plot(file2,position,vmin=20, vmax=93, cmap='hot',title='Chandra')

ax1 = fig.add_subplot(grid[0:, 1:])
ax2 = ax1.twinx()#核心
# ax1.plot(xp1,yp1,color='tab:blue',label='ASKAP data')
# ax1.plot(xp1,yp1_p,color='tab:blue',)
# yp1_l/og = np.log(yp1)
# print(eyp1)
# eyp1_log = np.log(eyp1)

# ax1.plot(xp1,yp1_p,color='tab:blue',)
# ax1.plot(xp1,yp1_log,color='tab:blue',)
# ax1.plot(xp1,yp1max,color='blue',label='ASKAP data max')
# eyp1 = eyp1/_range
# ax1.set_xlabel(r'angular distance (arcsec)')

# # ax2.plot(xp2,yp2,color="tab:red")#alpha=1.0)
# ax1.plot(xp2,yp2_p,color='tab:red',label='Chandra')
ax1.errorbar(xp1, yp1_p, yerr = eyp1, fmt='o', color='b', capsize=5, label='Data with Error Bars')
# ax1.set_xlim(-0.5,120)
# ax1.tick_params(axis='y')
# ax1.set_ylabel(r'normalised intensity (a.u.)')
# x_major_locator = MultipleLocator(10)
# ax1.xaxis.set_major_locator(x_major_locator)

# ax1.set_xscale("log", nonpositive='clip')
# ax1.set_yscale("log", nonpositive='clip')
# ax2.plot(xp1,yp1*1e3,color='tab:blue',label='Chandra',alpha=0)
# ax2.set_ylabel(r'radio intensity ($mJy$)')

ax1.errorbar(xp2, yp2_p, yerr =  eyp2, fmt='o', color='tomato', capsize=5, label='Data with Error Bars')
# ax3 = fig.add_axes([0.4695, 0.5, 0.4305, 0.38]) # inside axes
# # ax4 = ax3.twinx()
# # ax.inset_axes([0.6, 0.6, 0.35, 0.35])


# ax3.plot(xp1,yp1_p,color="tab:blue",label='ASKAP(mean)',)
# ax3.errorbar(xp1, yp1_p, yerr=eyp1, fmt='o', color='b', capsize=5, label='Data with Error Bars')
# ax3.plot(xp2,yp2_p,color="tab:red",label='Chandra(mean)',)#绿色看右边y轴。
# ax3.set_xlim(30,120)
# ax3.set_ylim(-.001, 0.04)
# # ax3.set_yticklabels([])
# # ax3.set_title('zoom in',y=-0.3, fontsize=20)
ax1.set_ylim(0,0.3)
# # axvline
# ax1.axvline(r1_shell,color=color_cir,linestyle=line_cir,linewidth=1)
# ax1.axvline(r2_shell,color=color_cir,linestyle=line_cir,linewidth=1)
# # ax1.axvline(annu,color='yellow',)

# ax3.axvline(r1_shell,color=color_cir,linestyle=line_cir,linewidth=1)
# ax3.axvline(r2_shell,color=color_cir,linestyle=line_cir,linewidth=1)
# # ax3.axvline(annu,color='magenta')


# ax1.text(5,1.0,'PSR')
# ax1.text(40,0.2,'PWN')
# ax1.text(70,0.2,'shell')
# ax1.text(100,0.2,'background')
# ax3.legend(loc='upper right')

# plt.savefig('plotProf_fan_mean_zoom.pdf',format='pdf', bbox_inches='tight')

plt.show()
# %%
