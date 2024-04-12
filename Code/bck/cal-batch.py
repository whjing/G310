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

from subprocess import run




maskNameList = ["SNR", "PWN"]#, "ERS1"]
fileList = sorted(glob("../Fits_for-4-subband/image.i.EMU_1356-64.SB53310.beam27.spw-*.taylor.0.restored.sm_crop-0.4.fits", ))
bkgFile = (f"../Fits_for-4-subband/bkg_mask.fits")
bkgData = (f"../Fits_for-4-subband/bkg_mask.w.dat")


for maskName in maskNameList:
    for file in fileList:
        print(file)
        maskFile = f"../Fits_for-4-subband/{maskName}_mask.fits"
        maskData = f"../Fits_for-4-subband/{maskName}_mask.w.dat"
        run([
            "calIntFlux_polygon_simple.py",
            file,
            maskData,
            bkgData,
            maskFile,
            bkgFile,
            
        ],)



# %%
