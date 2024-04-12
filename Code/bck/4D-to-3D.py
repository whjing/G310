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
from copy import deepcopy
from glob import glob

fileList = glob("../RawImage/image.i.cube.cutout.total.fits")
for file in fileList:
    try:

        hdr = pf.getheader(file)
        hdr0 = deepcopy(hdr)
        data = pf.getdata(file)[:,0,:,:]
        # modify the hdr0
        print(repr(hdr0))
        del hdr0["NAXIS3"]
        del hdr0["CTYPE3"]
        del hdr0["CRVAL3"]
        del hdr0["CDELT3"]
        del hdr0["CRPIX3"]
        del hdr0["CUNIT3"]

        hdr0["NAXIS3"] = hdr0["NAXIS4"]
        hdr0["CTYPE3"] = hdr0["CTYPE4"]
        hdr0["CRVAL3"] = hdr0["CRVAL4"]
        hdr0["CDELT3"] = hdr0["CDELT4"]
        hdr0["CRPIX3"] = hdr0["CRPIX4"]
        hdr0["CUNIT3"] = hdr0["CUNIT4"]

        del hdr0["NAXIS4"]
        del hdr0["CTYPE4"]
        del hdr0["CRVAL4"]
        del hdr0["CDELT4"]
        del hdr0["CRPIX4"]
        del hdr0["CUNIT4"]

        hdu = pf.PrimaryHDU(data, hdr0)
        # hdu = 
        print(repr(hdu.header))
        hdu.writeto(file.replace(".fits", "_3D.fits"), overwrite=True)
    except Exception as e:
        print(e)






#%%

file  = "../Fits/FDF_tot_dirty.fits"
hdr = pf.getheader(file)
wcs = WCS(hdr)
print(repr(hdr))
data = pf.getdata(file)[:, 0, 610, 593]
data_img = pf.getdata(file)[40, 0, :, :]
print(data_img.shape)
RM = hdr['CRVAL4'] + (np.arange(hdr['NAXIS4'])+1 - hdr['CRPIX4'])*hdr['CDELT4']
print(freq)
# print(RM_arr)
plt.plot(RM, data)
fig = plt.figure()
ax = fig.subfigures(111)
im = ax.imshow(data_img)
# %%
