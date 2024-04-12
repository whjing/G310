# %%
import matplotlib.pyplot as plt
from astropy.io import fits as pf
import numpy as np
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from reproject import reproject_interp
import sys
from copy import deepcopy



def regrid(file_o, file_t):
    hdul = pf.open(file_o)
    hdu_o = hdul["PRIMARY"]
    hdr_o = hdu_o.header
    data_o = hdu_o.data

    with pf.open(file_t) as hdul:
        data_t = hdul["PRIMARY"].data
        hdr_t = hdul["PRIMARY"].header

    hdr_n = deepcopy(hdr_o)
    hdr_n["NAXIS1"] = hdr_t["NAXIS1"]
    hdr_n["NAXIS2"] = hdr_t["NAXIS2"]
    hdr_n["CRVAL1"] = hdr_t["CRVAL1"]
    hdr_n["CRVAL2"] = hdr_t["CRVAL2"]
    hdr_n["CDELT1"] = hdr_t["CDELT1"]
    hdr_n["CDELT2"] = hdr_t["CDELT2"]
    hdr_n["CRPIX1"] = hdr_t["CRPIX1"]
    hdr_n["CRPIX2"] = hdr_t["CRPIX2"]
    try:
        hdr_n["CROTA1"] = hdr_t["CROTA1"]
        hdr_n["CROTA2"] = hdr_t["CROTA2"]
    except:
        print("No need to rotate. ")
        pass


    data_n, footprint = reproject_interp(hdu_o, hdr_n)
    hdul.close()
    pf.writeto("test.fits", data_n, hdr_n, overwrite = True)
    return hdr_n, data_n



def smooth(data_in, hpbw_o, hpbw_n, hdr_n):
    """
    smoothing, all angles are in arcmin

    hpbw_o and hpbw_n in arcsec

    """
    if hpbw_n < hpbw_o:
        print("New beam less than the old one, no smoothing!")
        sys.exit()

    else:
        hpbw = np.sqrt(hpbw_n ** 2 - hpbw_o ** 2)  # in arcmin
        
        grid = abs(hdr_n["CDELT1"]) * 3600

        g = Gaussian2DKernel(hpbw / (2.0 * np.sqrt(2.0 * np.log(2.0))) / grid)

        data_out = convolve(data_in, g, boundary="extend")

        return data_out
    


def make_smoothedFits(file_o, file_t, hpbw_n, hpbw_o, outFile=None, plot = False):
    """
    
    """
    hdr_n, data_n = regrid(file_o, file_t)

    # 射电的数据还要再乘上一个因子，这个因子意思目前还不知道
    # data = data_o * (hpbw_n**2) / (hpbw_o**2)
    data_sm = smooth(data_n, hpbw_o, hpbw_n, hdr_n)
    if outFile != None:
        pf.writeto(outFile, data_sm, hdr_n, overwrite=True)

    """
    想画4个图，1 模板 2 原始图像 3 网格后 4 平滑后
    """
    if plot !=False:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True)
        im1 = ax1.imshow(data_n)
        im2 = ax2.imshow(data_sm)
        plt.show()

# file_o = "../Chandra_process/process/flux_1/hard_flux.img"
file_o = "../Chandra_process/process/flux_2/hard_flux.img"
file_t = "../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-250.fits"

hpbw_n = 12
hpbw_o = 0.5

# ciao bin per 2arcsec
# hpbw_o = 2

outFileName = "../Fits_4-subband-useful/hard_flux.img.smooth_bin-2-hpbw-0.5.fits"


make_smoothedFits(file_o, file_t, hpbw_n, hpbw_o, outFile=outFileName, plot = True)
# if __name__ == "__main__":
#     main()

# %%

from astropy.io import fits as pf
file_o = "../Chandra_process/process/flux_1/hard_flux.img_crop-250.fits"
file_t = "../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-250.fits"

newimage =  pf.hcongrid(pf.open(file_o)[0].data, pf.open(file_o)[0].header, pf.open(file_t)[0].header)
# %%
from astropy.io import fits as pf
import FITS_tools
import matplotlib.pyplot as plt
from astropy.wcs import WCS

file_o = "../Chandra_process/process/flux_1/hard_flux.img_crop-250.fits"
file_t = "../Fits_4-subband-useful/image.i.EMU_1356-64.SB53310.beam27-28-33.taylor.0.restored_crop-250.fits"

data_o = pf.open(file_o)[0].data
header_o =  pf.open(file_o)[0].header
wcs = WCS(header_o)
ctype = wcs.equinox
print(ctype)

# header_t =  pf.open(file_t)[0].header

# newimage =  FITS_tools.hcongrid.hcongrid(data_o, header_o, header_t)


# fig = plt.figure()
# ax = fig.add_subplot(111)
# im = ax.imshow(newimage ,  cmap = 'jet_r')
# plt.show()
# %%
