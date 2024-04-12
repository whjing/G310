#%%
import os
import re
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from glob import glob
from astropy.io import fits as pf
import numpy as np
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from reproject import reproject_interp
from matplotlib.patches import Circle
from glob import glob
from matplotlib.patches import Ellipse, Rectangle
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import interpolate
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator
import sys
import matplotlib.ticker as ticker 
sys.path.append("~/academic/module/")


def load_fits_image(filename):

    # Read the header and image data from the file
    header=pf.getheader(filename)
    data = pf.getdata(filename)
    naxis = len(data.shape)

    # Strip unused dimensions from the data array
    if naxis == 2:
        xydata = data.copy()
        del data
    elif naxis == 3:
        xydata = data[0].copy()
        del data
    elif naxis == 4:
        xydata = data[0][0].copy()
        del data
    elif naxis == 5:
        xydata = data[0][0][0].copy()
        del data
    else:
        print("Data array contains %s axes" % naxis)
        print("This script supports up to 5 axes only.")
        sys.exit(1)

    # Strip unused dimensions from the header
    stripHeadKeys = ['NAXIS','CRVAL','CRPIX','CDELT','CTYPE','CROTA',
                     'CD1_','CD2_','CUNIT','PC1_','PC2_','PC3_','PC4_']
    
    stripHeadKeys1 = ['PC1_','PC2_','PC3_','PC4_','PC01_','PC02_','PC03_','PC04_']
    if naxis >= 2:
        for i in range(3,6):
            for key in stripHeadKeys:
                if (key+str(i)) in header: del header[key+str(i)]
        for j in range(1,6):            
            for key1 in stripHeadKeys1:
                if (key1+str(j)) in header: del header[key1+str(j)]
                if (key1+ '0' + str(j)) in header: del header[key1+ '0' + str(j)]
        header['NAXIS'] = 2

    # Determine the coordinate type and the corresponding keyword index
    # Make a note in the header
    ra_regx = re.compile('^RA')
    dec_regx = re.compile('^DEC')
    glon_regx = re.compile('^GLON')
    glat_regx = re.compile('^GLAT')
    if 'CTYPE1' in header:
        if ra_regx.match(header['CTYPE1']) or glon_regx.match(header['CTYPE1']):
            for i in range(int(header['NAXIS'])):
                keyword = "CTYPE"+str(i+1)
                if ra_regx.match(header[keyword]): coord_type="EQU"; x_indx=i+1
                if dec_regx.match(header[keyword]): y_indx=i+1
                if glon_regx.match(header[keyword]): coord_type="GAL"; x_indx=i+1
                if glat_regx.match(header[keyword]): y_indx=i+1
            if not x_indx or not y_indx:
                print("Failed to find Equatorial or Galactic axis coordinate types.")
                del data; del header
                sys.exit(1)
            else:
                header['XINDEX'] = x_indx
                header['YINDEX'] = y_indx
        else:
            print('Does not have "RA-DEC" and "GLON-GLAT" coordinates!!!')
    else:
        print('key "CTYPE1" does not, please check the header!!!')
    # Convert AIPS clean-beam types to standard BMAJ, BMIN, BPA
    try:
        bmaj = header['CLEANBMJ']
        bmin = header['CLEANBMN']
        bpa = header['CLEANBPA']
        header.update('BMAJ',bmaj)
        header.update('BMIN',bmin)
        header.update('BPA',bpa)
    except Exception:
        print("No AIPS style beam keywords found.")

    # Check for PIXSCAL keywords and write to CDELT standard
    try:
        xdelt=(-1)*(header['PIXSCAL'+str(x_indx)])/3600.0
        ydelt=(header['PIXSCAL'+str(y_indx)])/3600.0
        header['CDELT'+str(x_indx)] = xdelt
        header['CDELT'+str(y_indx)] = ydelt
    except Exception:
        pass


    wcs = WCS(header)

    return [header,xydata,wcs]


# from myfits import load_fits_image, crop



# Basic Setup

plt.rcParams['font.size']
plt.rcParams['font.size'] = 20


plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内

plt.rc('font',family='Times New Roman') 

EMU_file = '../RawImage/image.i.EMU_1356-64.SB53310.cont.taylor.0.restored.conv.fits'
hdr, data, wcs = load_fits_image(EMU_file)



Taylor1_file = '../RawImage/image.i.EMU_1356-64.SB53310.cont.taylor.1.restored.conv.fits'

hdr_1, data_1, wcs_1 = load_fits_image(Taylor1_file)

data_spec = data_1 / data 

# linner scale
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection = wcs)
vmin = -1#np.nanpercentile(data, 1)
vmax = 1 #np.nanpercentile(data, 90)
# im = ax.imshow(data_spec, origin='lower', cmap = 'Spectral', vmin=vmin, vmax=vmax)
pf.writeto("../RawImage/SpecIndex.fits",data_spec, hdr)

# plt.savefig("SpecIndex_total" + '.png', bbox_inches='tight', dpi = 500)
plt.show()
# set_ax(ax)



# %%
import os
import re
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from glob import glob
from astropy.io import fits as pf
import numpy as np
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from reproject import reproject_interp
from matplotlib.patches import Circle
from glob import glob
from matplotlib.patches import Ellipse, Rectangle
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import interpolate
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator
import sys
import matplotlib.ticker as ticker 
sys.path.append("~/academic/module/")



# from myfits import load_fits_image, crop



# Basic Setup


plt.rcParams['font.size'] = 22


plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内

plt.rc('font',family='Times New Roman') 

# Function

def smooth(data_in, hpbw=2.):
    g = Gaussian2DKernel(hpbw)
    data_out = convolve(data_in, g, boundary='extend')
    return data_out

def show_ctr(ax, data_orig, hdr_ctr, hdr_img, rms, num = 6, hpbw = 2., color = 'tomato', linewidths = 2):
    data_s = smooth(data_orig, hpbw)
    levels=[]
    for i in range(num): 
        level = rms* 3 * (2** (i))
        levels.append(level)
    data_ctr, _ = reproject_interp(pf.PrimaryHDU(data_s, hdr_ctr), hdr_img)
    contour = ax.contour(data_ctr, levels=levels, colors = color, linewidths = linewidths)
    return contour


def show_cb(im, vmin, vmax, stp = 7, n = 1000, label = r'mJy beam$^{-1}$'):
    cb = plt.colorbar(im,  pad=0, fraction=0.0476, location='right')#orientation='horizontal')fraction=0.05,
    colors_range = np.linspace(vmin, vmax, stp)
    cb.set_ticks(colors_range)
    tick_label = []
    for l in range(stp):
        tick_label.append(round(colors_range[l] * n, 1))
    cb.set_ticklabels(tick_label)
    cb.ax.tick_params(labelsize=22)
    cb.set_label(label)
    return cb

def set_ax(ax, num_ticks=5, minpad = 1 ):
    ra = ax.coords[0]
    dec = ax.coords[1]

    ra.set_axislabel('Right Ascension (J2000)',minpad = minpad)
    dec.set_axislabel('Declination (J2000)', minpad = minpad,)
    ra.set_major_formatter('dd:mm:ss.s')
    dec.set_major_formatter('dd:mm')
    ax.tick_params(axis='ra', colors='white')
    ax.tick_params(axis='dec', colors='white')
    ra.set_ticklabel(color='black', size=15, pad = 10)
    dec.set_ticklabel(color='black', size=15,)


    ra.display_minor_ticks(True)
    dec.display_minor_ticks(True)
    ra.set_minor_frequency(10)
    dec.set_minor_frequency(10)
    

    # ax.spines['top'].set_color('white') 
    # ax.spines['right'].set_color('white')
    # ax.spines['left'].set_color('black')


    # tick_spacing = 0.04
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(tick_spacing))
    
    
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_ticks(number = num_ticks)
    dec.set_ticks(number = num_ticks)
    dec.set_ticklabel_visible(True)
    dec.set_ticklabel(rotation='vertical', va='center')

def sv_fig(name):
    sv = plt.savefig(name + '.png', bbox_inches='tight', dpi = 300)
    sv = plt.savefig(name + '.pdf', bbox_inches='tight')
    return sv



def beampatch(hdr, color = 'yellow'):
    pix_size1 = abs(hdr['CDELT1'])
    pix_size2 = abs(hdr['CDELT2'])
    p_pix = [0.005 / pix_size1, 0.005 / pix_size2]
    beam_size_x = hdr['BMAJ'] / pix_size1
    beam_size_y = hdr['BMIN'] / pix_size2
    beam = Ellipse((p_pix[0], p_pix[1]), beam_size_x, beam_size_y, angle=hdr['BPA'], facecolor='none', edgecolor=color, alpha=1, zorder=200, linewidth=2)
    res=np.sqrt(hdr['BMAJ'] *hdr['BMIN'])*3600
    print('resolution'+str(res))
    return beam

##################################
###### EMU data 
##################################




# %%
