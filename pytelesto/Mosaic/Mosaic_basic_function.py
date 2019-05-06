"""
MOSAIC_BASIC_FUNCTION: basic function used in the creation
of the Mosaic.
"""

from astropy.io import fits


def save_image(image, name, save):
    """Save image in a new fits file to be open with ds9.

    Do not return anything as it is only a succession of actions.

    Args:
        image (array-like): image that we want to save in a new
            file

        name (string): name of the future file
    """
    if save:
        name = name + '.fit'
        hdu = fits.PrimaryHDU(image)
        hdul = fits.HDUList([hdu])
        hdul.writeto(name, overwrite=True)