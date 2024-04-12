#%%
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits

# Load your image data
# Assuming you have FITS files for each color channel, replace 'red.fits' and 'blue.fits' with your filenames
red_data = fits.getdata('../Fits/G310_i_EMU_r250.fits')
blue_data = fits.getdata('../Fits/G310_i_Chandra_r250_smooth.fits')

# Normalize the data if necessary
red_data /= np.max(red_data)
blue_data /= 0.8 *np.max(blue_data)

# Create a green channel (for a simple example)
# You can manipulate your data to create the green channel according to your requirements
green_data = (red_data + blue_data) / 2

# Combine the three channels into a single RGB image
rgb_image = make_lupton_rgb(red_data, green_data, blue_data, stretch=0.01)

# Display the resulting image
plt.imshow(rgb_image, origin='lower')
plt.axis('off')
plt.show()

# %%
