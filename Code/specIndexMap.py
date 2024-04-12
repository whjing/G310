
#%%

import numpy as np
from astropy.io import fits as pf
from astropy.wcs import WCS
from matplotlib import pylab as plt
from matplotlib.colors import LogNorm
from reproject import reproject_interp
from copy import deepcopy
from astropy.coordinates import SkyCoord
from regions import EllipseSkyRegion, EllipsePixelRegion
from astropy import units as u
from astropy.modeling import models, fitting
from astroquery.vizier import Vizier
import scipy.optimize
from scipy import optimize
from glob import glob
import sys
from globalSetting import *
import warnings

# Suppress all warnings
warnings.filterwarnings("ignore")


#%%
def make_cube(fileNameList, freq_all = [], writeto = False):
    """
    Make a cube file:
    1. f
    """
    data_all = []
    fileList = sorted(glob(fileNameList))
    print("File list is: \n" + '\n'.join(fileList))
    for fileNameNow in fileList:
        print(f"Now we are working on the file: \n {fileNameNow}")
        
        hdr_orig = pf.getheader(fileNameNow)
        header, data, wcs = load_fits_image(fileNameNow)
        data_all.append(data)
        
        if 'CRVAL3' in hdr_orig:
            freq = hdr_orig['CRVAL3']
            print(fileNameNow, freq)
            freq_all.append(freq)
        elif 'freq' in hdr_orig:
            freq = hdr_orig['freq']
            freq_all.append(freq)
        else:
            print('Please provide the freqlist')

        # Draw map
        ###########################
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(1,1,1, projection = wcs)
        vmin, vmax = np.nanpercentile(data, 0), np.nanpercentile(data,99.5)
        im = ax.imshow(data, origin='lower', cmap = 'viridis', vmin=vmin, vmax=vmax)  
        set_ax(ax)
        print(freq_all)
        show_cb(im, vmin, vmax) 
        ax.set_title(fileNameNow)
        plt.show()
        ##########################
    # Write new fits file to 3 dimension cubes
    header, data, _ = load_fits_image(fileList[0])
    header['FREQLIST'] = str(freq_all)
    # data_all = np.array(data_all)
    hdu = pf.ImageHDU(data = data_all, header = header)
    header, data, wcs = hdu.header, hdu.data, WCS(hdu.header)
    if writeto == True:
        hdu.writeto('test.fits')
    return header, data, wcs


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
                if i==360 and j==360:
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


def plot_specIndexMap(fileNameList, freqData):
    '''
    f = '~/academic/fits/6cm-90cm_cub.fits'
    header = pf.getheader(f)
    data = pf.getdata(f)
    # 
    '''
    print(f"+++ 1. Make cube file. +++")
    header, data, wcs = make_cube(fileNameList) 

    print(f"+++ 2. Open the frequency file. +++")
    x_hz = freqData
    print(f"FreqList in cube is {x_hz}")

    print(f"+++ 3. Combine the spectral map data. +++")

    axis_x = header['NAXIS1']
    axis_y = header['NAXIS2']
    spectral_map_data = spectral_map(data, x_hz, axis_x, axis_y)

    # make new hdu
    hdu = pf.PrimaryHDU(data = spectral_map_data, header = header)
    header, spec_map_data, wcs = hdu.header, hdu.data, WCS(hdu.header)
    hdu.writeto('../Fits/specIndex_4-spws.fits', overwrite= True)

    # mask specific data
    spectral_map_data[np.where (spectral_map_data > 1)] = np.nan
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_subplot(111, projection = wcs, )#slices=('x','y',))
    vmin, vmax = -0.8, 0.8
    im = ax.imshow(spectral_map_data, vmin=vmin, vmax=vmax , cmap = 'Spectral',origin='lower')
    ax.plot(360, 360,'ro')

    show_cb(im, vmin, vmax, n = 1, label = ' ')
    # show_ctr(ax, data[0], header, header, rms = 0.008)

    # plt.savefig(path + file_name + '.png', dpi = 300)
    plt.show()
    ##############################################################




filePath = "../Fits_4-subband-useful"
fileNameCom = "image.i.EMU_1356-64.SB53310.beam27-28-33.spw-*.taylor.0.restored.total_crop-0.4.fits"
fileNameList = f"{filePath}/{fileNameCom}"
freqData = np.loadtxt("../Data/freq_4.dat", delimiter=",", dtype=np.float64, usecols=range(1))

print(f"freqData is {freqData}")
plot_specIndexMap(fileNameList, freqData)

# %%



# Basic Setup


plt.rcParams['font.size'] = 20


plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内

plt.rc('font',family='Times New Roman') 


file = "../Fits/specIndex_4-spws_crop-0.06.fits"
file_I = "../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.spw-1.taylor.0.restored.total_crop-0.4_crop-0.06.fits"

rms = 3e-5
hdr, data, wcs = load_fits_image(file)
hdr_I, data_I, wcs_I = load_fits_image(file_I)
fig = plt.figure(figsize=(12,9), dpi=100)
ax = fig.add_subplot(111, projection = wcs, )#slices=('x','y',))
vmin, vmax = -0.8, 0
#(data > 0) | (data < vmin)|
data[ np.where((data_I < rms * 30))] = np.nan
im = ax.imshow(data, vmin=vmin, vmax=vmax , cmap = 'Spectral',origin='lower')
show_cb(im, vmin, vmax, n = 1, label = ' ')
show_ctr(ax, data_I, hdr_I, hdr, rms, num  = 12, color = "grey")
set_ax(ax)
sv_fig("../figures/specIndexMap_PWN")
plt.show()


#%%



file_Mask_ERS1= "../Fits_4-subband-useful/ERS-1.mask.fits"
file_Mask_SNR= "../Fits_4-subband-useful/SNR.mask.fits"
file_I = "../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.spw-1.taylor.0.restored.total_crop-0.4.fits"

hdr_mask_ERS1, data_mask_ERS1, wcs_mask_ERS1 = load_fits_image(file_Mask_ERS1)
hdr_mask_SNR, data_mask_SNR, wcs_mask_SNR = load_fits_image(file_Mask_SNR)
hdr_I, data_I, wcs_I = load_fits_image(file_I)
# data_mask = np.where(np.isnan(data_mask_SNR), data_mask_ERS1, data_mask_SNR)
data = data_mask_ERS1 * data_I

fig = plt.figure(figsize=(12,9), dpi=100)
ax = fig.add_subplot(111, projection = wcs_I, )#slices=('x','y',))
vmin, vmax = np.nanpercentile(data_I, 5), np.nanpercentile(data_I, 95)
im = ax.imshow(data, vmin=vmin, vmax=vmax , cmap = 'Spectral',origin='lower')
show_cb(im, vmin, vmax, n = 1, label = ' ')
# show_ctr(ax, data_I, hdr_I, hdr, rms, num  = 12, color = "grey")
set_ax(ax)
sv_fig("../Output/specIndexMask_ERS1")

# %%
