#%%
from astropy.io import fits as pf
file = "../ref/components2020v2.fits"
hdul = pf.open(file)
# hdul.info()
data = hdul[1].data
# data = pf.getdata(f ie)
# %%
