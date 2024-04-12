#%%
from astropy.io import fits as pf
import numpy as np




f = "../RawImage/image.q.cube.cutout.fits"
hdr = pf.getheader(f)
hdr
freq_stp = int(hdr['CDELT4'])
freq_num = int(hdr['NAXIS4'])
freq_range =  freq_num * freq_stp 




data = pf.getdata(f)[:,0,0,0]
# print(data)
freq_start = int(hdr['CRVAL4'])
freq_end = freq_start + freq_range
freq_array = np.linspace(freq_start, freq_end, freq_num )
# print(freq_array.shape)
# freq_array_txt = str(freq_array)
# np.savetxt("../Data/freq.test.dat", freq_array)
# with open ("../Data/freq.test.dat", "w") as f:
    # f.write(freq_array_txt)

# %%
print(repr(hdr))
freq = hdr['CRVAL4'] + (np.arange(freq_num)+1 - hdr['CRPIX4'])*hdr['CDELT4']

freq.tofile('../Data/freq.dat',sep='\n')
# %%
