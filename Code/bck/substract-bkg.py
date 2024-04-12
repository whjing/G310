#%%
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
import pyregion
sys.path.append("/Users/jing/academic/module/")
print(sys.path)
from myfits import *

file_I = '../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-0.1.fits'
file_bkg = "../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-0.1_bkg.fits"

hdr_I, data_I, wcs_I = load_fits_image(file_I)
hdr_bkg, data_bkg, wcs_bkg = load_fits_image(file_bkg)

data = data_I - data_bkg


pf.writeto("test.fits", data, hdr_I, overwrite=True)
# %%
