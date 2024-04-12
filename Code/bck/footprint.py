#%%
import pandas as pd
import sys
sys.path.append("/Users/jing/academic/module")
from astropy.io import fits as pf
from myfits import load_fits_image, crop, set_ax
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.patches import Circle


def shell(ra,dec,r,color="red"):
    shell_c = SkyCoord(ra,dec, frame='fk5',unit=(u.hourangle, u.deg))
    r = r/((3600*(abs(hdr['CDELT2']))))
    x, y = wcs.wcs_world2pix(shell_c.ra.deg, shell_c.dec.deg, 0)
    print(x,y)
    shell_image = Circle((x,y), r, facecolor='none', edgecolor=color, alpha=1, zorder=200, linewidth = 2)
    ax.add_patch(shell_image)
    return ax
def sv_fig(name):
    sv = plt.savefig(name + '.png', bbox_inches='tight', dpi = 300)
    sv = plt.savefig(name + '.pdf', bbox_inches='tight')
    return sv


footprint_file = "../MetaData/calibration-metadata-processing-logs-SB53310_2023-12-14-150635/metadata/footprintOutput-sb53310-src1-EMU_1356-64.txt"

fits_file = "../RawImage/image.i.EMU_1356-64.SB53310.cont.taylor.0.restored.conv.fits"


hdr, data, wcs = load_fits_image(fits_file)


freq = 0.94
xlam = 0.3 / freq
D = 22.
beam_size = xlam/D*180./np.pi
# print(beam)
theta = xlam/(np.sqrt(3.)*D)*180./np.pi
r = beam_size /(abs(hdr['CDELT1']))
print(r)

fig = plt.figure(figsize = (10, 10))
ax = fig.add_subplot(1,1,1, projection = wcs)
vmin = np.nanpercentile(data, 1)
vmax = np.nanpercentile(data, 99)
im = ax.imshow(data, origin='lower', cmap = 'summer', vmin=vmin, vmax=vmax, alpha = 0.8) 

with open(footprint_file, "r") as f:
    lines = f.read().split("\n")[:-1]
    for line in lines:
        print(line)
        beam_name = line[0:2]
        coord = line.split()[-1].split(",")
        ra = coord[0]
        dec = coord[1]
        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
        ra = coord.ra
        dec = coord.dec
        x, y = wcs.wcs_world2pix(ra, dec, 0)
        # beam_size = (hdr['CDELT1'] + )/2

        if beam_name in ["27", "28", "33"]:
            color = 'red'
            color_anno = "black"
            zorder = 100
        else:
            color = "pink"
            color_anno = "white"
            zorder = 0
        circle = Circle((x, y), radius =r, facecolor='none', edgecolor = color, linewidth = 2, zorder = zorder)
        ax.add_patch(circle)
        ax.annotate(beam_name, (x,y), xytext=(1, 1), textcoords='offset points',color = color_anno, fontsize = 22, zorder = 1000)

pt0 = ('14:00:45', '-63:25:45')
c0 = SkyCoord(ra=pt0[0], dec=pt0[1], frame='fk5',unit=(u.hourangle, u.deg))
r = 200
shell(c0.ra, c0.dec, r, color = "cyan", )
print(data.shape)

set_ax(ax)

sv_fig('../Output/footprint-all')

# ax.set_xlim(5000, 9000)
# ax.set_ylim(7800, 11800)
# # sv_fig('../Output/footprint-zoom')
# %%
