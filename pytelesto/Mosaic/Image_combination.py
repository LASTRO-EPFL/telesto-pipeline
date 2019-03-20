"""
IMAGE_COMBINATION /Combination of all the images of a same portion of the
sky but with different exposure time, into one single images (mean of the
images pixel by pixel).
"""

##### Library #########################################################

import numpy as np
import astropy.io as astro
import glob
import astroalign as aa


from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn

##### Teatment of the images #########################################

def image_combination(k):
    name_all_images = glob.glob('.\February 12 2019\patch%d*.fit' % k)

    image_file0 = get_pkg_data_filename(name_all_images[0])
    image_data0 = astro.fits.getdata(image_file0, ext=0)

    for i in range(1,len(name_all_images)-1):
        image_file = get_pkg_data_filename(name_all_images[i])
        image_data = astro.fits.getdata(image_file, ext=0)
        image_data0 = (image_data0 + image_data)/2

    hdu = fits.PrimaryHDU(image_data0)
    hdul = fits.HDUList([hdu])
    hdul.writeto('patch%d_combined.fit' % k)

    return image_data0

image1 = image_combination(1)
image2 = image_combination(2)
image3 = image_combination(3)
image4 = image_combination(4)

