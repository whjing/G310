#%%
import numpy as np
from astropy.io import fits as pf


file = "/Users/jing/G310/MetaData/calibration-metadata-processing-logs-SB53310_2023-12-14-150635/scratch/askaprt/askapops/science-processing/53310/LinmosBeamImages/akpb.iquv.closepack36.54.943MHz.SB51811.taylor.fits"

data = pf.getdata(file)

print(data.shape)

# %%
