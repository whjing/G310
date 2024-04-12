#!/usr/bin/env python
#%%
from astropy.io import fits as pf
from reproject import reproject_interp
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.wcs import WCS
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from astropy.visualization import SqrtStretch
from astropy.visualization import ZScaleInterval
import numpy as np
from astropy.coordinates import SkyCoord

plt.rcParams.update({'font.size': 22})

def smooth(data_in, hpbw=2.):
    """
    smoothing
    """
    g = Gaussian2DKernel(hpbw)

    data_out = convolve(data_in, g, boundary='extend')

    return data_out


def plotRGB(r_file, g_file, b_file, w_file, sm_r = 1., sm_g=1., sm_b=1., lo=10., up=90., pol=False, figname=None):
    hdr = pf.getheader(w_file)
    wcs = WCS(hdr)

    if r_file == w_file:
        data_r_regrid = pf.getdata(r_file)
        if pol:
            mask1 = pf.getdata('mask.P.bad.1.fits')
            mask2 = pf.getdata('mask.P.bad.2.fits')
            mask = mask1 + mask2
            data_r_regrid[mask > 0.5] = np.nan

    else:
        hdu_r = pf.open(r_file)[0]
        hdu_r.data = smooth(hdu_r.data, hpbw=sm_r)
        data_r_regrid, footprint = reproject_interp(hdu_r, hdr)

    if g_file == w_file:
        data_g_regrid = pf.getdata(g_file)
        if pol:
            mask1 = pf.getdata('mask.P.bad.1.fits')
            mask2 = pf.getdata('mask.P.bad.2.fits')
            mask = mask1 + mask2
            data_g_regrid[mask > 0.5] = np.nan
    else:
        hdu_g = pf.open(g_file)[0]
        hdu_g.data = smooth(hdu_g.data, hpbw=sm_g)
        data_g_regrid, footprint = reproject_interp(hdu_g, hdr)

    if b_file == w_file:
        data_b_regrid = pf.getdata(b_file)
        if pol:
            mask1 = pf.getdata('mask.P.bad.1.fits')
            mask2 = pf.getdata('mask.P.bad.2.fits')
            mask = mask1 + mask2
            data_b_regrid[mask > 0.5] = np.nan
    else:
        hdu_b = pf.open(b_file)[0]
        hdu_b.data = smooth(hdu_b.data, hpbw=sm_b)
        data_b_regrid, footprint = reproject_interp(hdu_b, hdr)

    forCasting = np.float_()
    r = np.array(data_r_regrid, forCasting)
    g = np.array(data_g_regrid, forCasting)
    b = np.array(data_b_regrid, forCasting)
    r_tmp = r[~np.isnan(r)]
    r[np.isnan(r)] = r_tmp.min() 
    g_tmp = g[~np.isnan(g)]
    g[np.isnan(g)] = g_tmp.min() 
    b_tmp = b[~np.isnan(b)]
    b[np.isnan(b)] = b_tmp.min() 

    stretch = SqrtStretch() + ZScaleInterval()
    r = stretch(r)
    g = stretch(g)
    b = stretch(b)

    lo_val, up_val = np.percentile(np.hstack((r.flatten(), g.flatten(), b.flatten())), (lo, up))  # Get the value of lower and upper 0.5% of all pixels

    stretch_val = up_val - lo_val
    rgb = make_lupton_rgb(r, g, b, minimum=lo_val, stretch=stretch_val, Q=0)

    if figname == None:
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(2,2,1,projection=wcs)
        ax.imshow(r, origin='lower')
        ax = fig.add_subplot(2,2,2,projection=wcs)
        ax.imshow(g, origin='lower')
        ax = fig.add_subplot(2,2,3,projection=wcs)
        ax.imshow(b, origin='lower')
        ax = fig.add_subplot(2,2,4,projection=wcs)
        ax.imshow(rgb, origin='lower')
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        plt.show()
    else:
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(1,1,1, projection=wcs)
        ax.imshow(rgb, origin='lower')
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')

        # c = SkyCoord(ra='20h58m52.8s', dec='31d38m52.00s', frame='fk5')
        # xt, yt = wcs.wcs_world2pix(c.ra.deg, c.dec.deg, 0)
        # ax.text(xt, yt, 'NGC 6992/6995',c='white', va='center', rotation=60., fontweight='bold')

        # c = SkyCoord(ra='20h44m24.0s', dec='30d26m36.00s', frame='fk5')
        # xt, yt = wcs.wcs_world2pix(c.ra.deg, c.dec.deg, 0)
        # ax.text(xt, yt, 'NGC 6960',c='white', va='center', rotation=-90., fontweight='bold')

        # c = SkyCoord(ra='20h46m55.0s', dec='32d01m22.00s', frame='fk5')
        # xt, yt = wcs.wcs_world2pix(c.ra.deg, c.dec.deg, 0)
        # ax.text(xt, yt, "Pickering's triangle",c='white', va='center', rotation=-45.,fontweight='bold')

        plt.tight_layout()
        plt.savefig(figname, bbox_inches='tight')
        plt.close()


def main():
    
    
    w_file = 'CygnusLoop_I_med_beam_info.fits'

    """
    r_file = 'cygnusloop_uv.fits'
    g_file = 'CygnusLoop_I_med_beam_info.fits'
    b_file = 'cygnusloop_Xray.fits'
    figname='cyg_uv_x_radio.pdf'
    plotRGB(r_file, g_file, b_file, w_file, lo=5, up=95, figname=figname)
    """
    
    """
    r_file = 'cygnusloop_wise_22.fits'
    g_file = 'CygnusLoop_I_med_beam_info.fits'
    b_file = 'cygnusloop_353.fits'
    figname='cyg_wise22_radio_353.pdf'
    plotRGB(r_file, g_file, b_file, w_file, sm_r=2., sm_b=0.5, lo=5, up=95, figname=figname)
    """

    """
    r_file = 'cygnusloop_wise_22.fits'
    g_file = 'CygnusLoop_I_med_beam_info.fits'
    b_file = 'cygnusloop_iras_100.fits'
    figname='cyg_wise_22_radio_100.pdf'
    plotRGB(r_file, g_file, b_file, w_file, sm_r=2., sm_b=0.5, lo=5, up=95, figname=figname)
    """
    
    """
    r_file = 'cygnusloop_dss2_r.fits'
    g_file = 'CygnusLoop_I_med_beam_info.fits'
    b_file = 'cygnusloop_Xray.fits'
    figname='cyg_dss_x_radio.pdf'
    plotRGB(r_file, g_file, b_file, w_file, sm_r=4., lo=5, up=95, figname=figname)
    """
    
    
    
    # w_file = 'CygnusLoop_p_cln_flag_0_clip_edit.fits'
    w_file = '../Fits/G310_i_Chandra_r250.fits'
    pol = False

    """
    r_file = 'cygnusloop_wise_22.fits'
    g_file = 'CygnusLoop_p_cln_flag_0_clip_edit.fits'
    b_file = 'cygnusloop_353.fits'
    figname='cyg_wise22_p_353.pdf'
    plotRGB(r_file, g_file, b_file, w_file, sm_r=4., lo=5, up=95, pol=pol, figname=figname)
    """

    # r_file = 'cygnusloop_uv.fits'
    r_file = '../Fits/G310_i_Chandra_r250_smooth.fits'
    # g_file = 'CygnusLoop_p_cln_flag_0_clip_edit.fits'
    g_file = '../Fits/G310_i_EMU_r250.fits'
    # b_file = 'cygnusloop_Xray.fits'
    b_file = '../Fits/G310_i_EMU_r250.fits'
    figname='../Output/RGB.png'
    plotRGB(r_file, g_file, b_file, w_file, lo=5, up=95, pol=pol, figname=figname)
    
    """
    r_file = 'cygnusloop_dss2_r.fits'
    g_file = 'CygnusLoop_p_cln_flag_0_clip_edit.fits'
    b_file = 'cygnusloop_Xray.fits'
    figname='cyg_dss_x_p.pdf'
    plotRGB(r_file, g_file, b_file, w_file, sm_r=4., lo=5, up=95, pol=pol, figname=figname)
    """

if __name__ == "__main__":
    main()


# %%
