"""
MOSAIC_OF_IMAGES2 : from 4 images that superpose partly, creates a fifth,
bigger, image that is a combination of the previous ones.
Importance of the PLACE of the four images in the mosaic. Name them such
that :
  - image_data1 = top right of the mosaic
  - image_data2 = top left of the mosaic
  - image_data3 = bottom left of the mosaic
  - image_data4 = bottom right of the mosaic
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

##### Data importation ###############################################

# import image_combination.py
"""
Idealy, import image_combination into this program so everything run in
a single time, but compilation time of image_combination being big, we
directly import the final combined images into the program.
"""

image_file1 = get_pkg_data_filename('.\success_images\patch1_combined.fit')
image_data1 = astro.fits.getdata(image_file1, ext=0)

image_file2 = get_pkg_data_filename('.\success_images\patch2_combined.fit')
image_data2 = astro.fits.getdata(image_file2, ext=0)

image_file3 = get_pkg_data_filename('.\success_images\patch3_combined.fit')
image_data3 = astro.fits.getdata(image_file3, ext=0)

image_file4 = get_pkg_data_filename('.\success_images\patch4_combined.fit')
image_data4 = astro.fits.getdata(image_file4, ext=0)

##### Functions ##################################################

def save_image(image,name):
    """Save image in a new fits file to be open with ds9

    Do not return anything as it is only succession of actions.

    Args:
        image (array-like): image that we want to save in a new
            file
        name (string): name of the future file
    """
    name = name + '.fit'
    hdu = fits.PrimaryHDU(image)
    hdul = fits.HDUList([hdu])
    hdul.writeto(name)

def big_image(image, place):
    """Add pixels set to defaults around initial image

    Return the new image (rectangle described below), placed in a much
    bigger rectangle that have approximately the shape of the final
    mosaic. And call save_image to save the new image.

    Args:
        image (array-like): initial image that we want to align
        place (string) : 'top_right' , 'top_left', 'bottom_right'
           'bottom_left', depending on where should approximately
            placed the image in the mosaic

    Returns:
        image_bigger : initial image with bigger size, all new pixels
            being set to default.

    """
    image_bigger = np.copy(image_with_bigger_size)

    image[:, 0:cut_boarder] = horizontal_boarder
    image[:, N[1] - cut_boarder:N[1]] = horizontal_boarder
    image[0:cut_boarder, :] = vertical_boarder
    image[(N[0] - cut_boarder):N[0], :] = vertical_boarder

    if place == 'top_right':
        image_bigger[(N[0] - 1 + int(add_edges / 2)):(2 * N[0] - 1 + int(add_edges / 2)),
        N[1] - 1 + int(add_edges / 2):(2 * N[1] - 1 + int(add_edges / 2))] = image
    if place == 'top_left':
        image_bigger[(N[0] - 1 + int(add_edges / 2)):(2 * N[0] - 1 + int(add_edges / 2)),
        int(add_edges / 2) - 1:(N[1] - 1 + int(add_edges / 2))] = image
    if place == 'bottom_left':
        image_bigger[int(add_edges / 2) - 1:(N[0] - 1 + int(add_edges / 2)),
        int(add_edges / 2) - 1:(N[1] - 1 + int(add_edges / 2))] = image
    if place == 'bottom_right':
        image_bigger[int(add_edges / 2) - 1:(N[0] - 1 + int(add_edges / 2)),
        N[1] - 1 + int(add_edges / 2):(2 * N[1] - 1 + int(add_edges / 2))] = image

    place = 'image_' + place
    save_image(image_bigger, place)

    return image_bigger


def alignment(source, target, name):
    """ Align 2 images at a time.

    Return the combination of source and target images
    aligned with each other. Firstly, align source on
    target. Then re-add target image on align image.
    Call save_image to save the final image.

    Args:
        source (array-like): A numpy array of the source image to be
            transformed, its pixels will match the ones of the target
            image
        target (array-like): A numpy array of the target (destination)
            image.
        name (string): Name of the future file, that is a saving of the
            final image.

    Returns:
        image_align : The new image, combination of the source and
            target image.

    """
    image_align = aa.register(source, target)

    for i in range(len(image_align)):
        for j in range(len(image_align[i])):
            if (target[i, j] != default) and (image_align[i, j] != default):
                image_align[i, j] = np.mean([image_align[i, j], target[i, j]])
            if (target[i, j] != default) and (image_align[i, j] == default):
                image_align[i, j] = target[i, j]

    name = 'image_align' + name
    save_image(image_align, name)

    return image_align


def cut_image(image):
    """Cut the full bands set to default in the image

    Return the first cut of the image. The bands fully equal to
    default along the total width or total length of the image
    are removed.

    Args:
        image (array-like): Mosaic of the images, sticks on
            a bigger rectangle with all pixels set to default.

    Return:
        image (array-like): Initial image for which the horizontal
            and vertical bands with value default have been
            removed.

    """
    M = image.shape
    i = 0
    # Starting from the left boarder, we count the number of
    # vertical line that have only default values in order to
    # remove them later. Same thing done starting from the right
    # boarder, top and bottom boarder.
    while (image[i, :] == vertical_default).all():
        i = i + 1

    j = M[0]-1
    while (image[j, :] == vertical_default).all():
        j = j - 1

    k = 0
    while (image[:, k] == horizontal_default).all():
        k = k + 1

    l = M[1]-1
    while (image[:, l] == horizontal_default).all():
        l = l - 1

    image_cut = np.copy(image[i:j, k:l])
    save_image(image_cut, 'image_cut')
    return image_cut


def final_cut(image):
    """Final purification of the image with all default values
        removed

    Return the mosaic image with the right boarders, all bands
    that still contain partly default values are removed.

    Args:
         image (array-like): mosaic of image that has been
            firstly cut with function cut image. The only
            bands that remain are the one partly equal to
            default.

    Return:
        image (array-like): Return the final image, a
            combination of 4 images that superpose partly
            and that has been reframed correctly.

    """

    # Starting from the bottom-left corner we count the number
    # of pixels on a diagonal of a square with default values.
    # We search for the maximal diagonal possible for the square
    # to be only composed of pixels set to default. Then, we
    # check the direction of the band to be removed. If its
    # vertical, we removed it and then redo the same process,
    # in case there is also an horizontal band to removed.
    if image[0, 0] == default:
        index = 1
        while image[index, index] == default:
            index = index + 1
        if image[index + 1, index] != default:
            image = image[:, index:-1]
            save_image(image, '1')
            if image[0, 0] == default:
                index = 1
                while image[index, index] == default:
                    index = index + 1
                image = image[index:-1, :]
                save_image(image, '2')
        else:
            image = image[index:-1, :]
            save_image(image, '3')
            if image[0, 0] == default:
                index = 1
                while image[index, index] == default:
                    index = index + 1
                image = image[:, index:-1]
                save_image(image, '4')

    if image[image.shape[0]-1, 0] == default:
        index1 = image.shape[0] - 2
        index2 = 1
        while image[index1, index2] == default:
            index1 = index1 - 1
            index2 = index2 + 1
        if image[index1 - 1, index2] != default:
            image = image[:, index2:-1]
            save_image(image, '5')
            if image[image.shape[0]-1, 0] == default:
                index1 = image.shape[0] - 2
                index2 = 1
                while image[index1, index2] == default:
                    index1 = index1 - 1
                    index2 = index2 + 1
                image = image[0:index1, :]
                save_image(image, '6')

        else:
            image = image[0:index1, :]
            save_image(image, '7')
            if image[image.shape[0]-1, 0] == default:
                index1 = image.shape[0] - 2
                index2 = 1
                while image[index1, index2] == default:
                    index1 = index1 - 1
                    index2 = index2 + 1
                image = image[:, index2:-1]
                save_image(image, '8')

    if image[0, image.shape[1]-1] == default:
        index1 = 1
        index2 = image.shape[1] - 2
        while image[index1, index2] == default:
            index1 = index1 + 1
            index2 = index2 - 1
        if image[index1 + 1, index2] != default:
            image = image[:, 0:index2]
            save_image(image, '9')
            if image[0, image.shape[1]-1] == default:
                index1 = 1
                index2 = image.shape[1] - 2
                while image[index1, index2] == default:
                    index1 = index1 + 1
                    index2 = index2 - 1
                image = image[index1:-1, :]
                save_image(image, '10')
        else:
            image = image[index1:-1 :]
            save_image(image, '11')
            if image[image.shape[0]-1, 0] == default:
                index1 = 1
                index2 = image.shape[1] - 2
                while image[index1, index2] == default:
                    index1 = index1 + 1
                    index2 = index2 - 1
                image = image[:, 0:index2]
                save_image(image, '12')

    if image[image.shape[0]-1, image.shape[1]-1] == default:
        index1 = image.shape[0] - 2
        index2 = image.shape[1] - 2
        while image[index1, index2] == default:
            index1 = index1 - 1
            index2 = index2 - 1
        if image[index1 - 1, index2] != default:
            image = image[:, 0:index2]
            save_image(image, '13')
            if image[image.shape[0]-1, image.shape[1]-1] == default:
                index1 = image.shape[0] - 2
                index2 = image.shape[1] - 2
                while image[index1, index2] == default:
                    index1 = index1 - 1
                    index2 = index2 - 1
                image = image[0:index1, :]
                save_image(image, '14')
        else:
            image = image[0:index1, :]
            save_image(image, '15')
            if image[image.shape[0]-1, image.shape[1]-1] == default:
                index1 = image.shape[0] - 2
                index2 = image.shape[1] - 2
                while image[index1, index2] == default:
                    index1 = index1 - 1
                    index2 = index2 - 1
                image = image[:, 0: index2]
                save_image(image, '16')

        save_image(image, 'mosaic')

    return image


# def final_cut(image):

    # while image[0, 0] == default:
    #     image = image[1:-1, 1:-1]
    # save_image(image,'1')
    # while image[image.shape[0]-1, 0] == default:
    #     image = image[0:image.shape[0]-2, 1:-1]
    # save_image(image,'2')
    # while image[0, image.shape[1]-1] == default:
    #     image = image[1:-1, 0:image.shape[1]-2]
    # save_image(image,'3')
    # while image[image.shape[0]-1, image.shape[1]-1] == default:
    #     image = image[0:image.shape[0]-2, 0:image.shape[1]-2]
    # save_image(image, 'final')
    # return image


##### Main ################################################################

"""
Firstly, each images is positioned in a rectangle of sides 2 times
the ones of the image, adding also a band (add_edges) all around
the rectangle such that the initial image doesn't touch the boarder
of the new rectangle where it is placed. This ensure that the
alignment of the 2 images can be done even if the 'source' one
has an offset in both x and y direction. The new pixels, surrounding
the initial image are fixed to a default value (must be much bigger
than 0 so it works).
"""

N = image_data1.shape  # Suppose all images have the same shape
add_edges = 400
default = 1000  # amélioration: idée, pour que ce ne soit pas vu comme du bruit, utiliser la moyenne des valeurs des pixels
                # de l'image comme défault
cut_boarder = 20  # The boarder of initial images are too noisy, so we cut them a bit to decrease the noise bands in
                  # between all images when thay will be aligned.

image_with_bigger_size = np.ones((2 * N[0] + add_edges, 2 * N[1] + add_edges), dtype='f')
n = image_with_bigger_size.shape

for i in range(len(image_with_bigger_size)):
    for j in range(len(image_with_bigger_size[i])):
        image_with_bigger_size[i, j] = default

# To cut the boarder of the initial images, replace it with a rectangular boarder with all pixels
# set to default:
horizontal_boarder = np.zeros((N[0], cut_boarder), dtype='f')
vertical_boarder = np.zeros((cut_boarder, N[1]), dtype='f')
for i in range(len(horizontal_boarder)):
    for j in range(len(horizontal_boarder[i])):
        horizontal_boarder[i, j] = default

for i in range(len(vertical_boarder)):
    for j in range(len(vertical_boarder[i])):
        vertical_boarder[i, j] = default

# not necessary : can all be implemented so each function is calling the other
image_top_right = big_image(image_data1, 'top_right')
image_top_left = big_image(image_data2, 'top_left')
image_bottom_left = big_image(image_data3, 'bottom_left')
image_bottom_right = big_image(image_data4, 'bottom_right')

image_align1_2 = alignment(image_top_left, image_top_right, '1_2')
image_align1_2_3 = alignment(image_bottom_left, image_align1_2, '1_2_3')
image_align1_2_3_4 = alignment(image_bottom_right, image_align1_2_3, '1_2_3_4')

N = image_align1_2_3_4.shape
horizontal_default = np.zeros((N[0], 1), dtype='f')
vertical_default = np.zeros((1, N[1]), dtype='f')

for i in range(len(horizontal_default)):  # method can be reduced, using previous horizontal band created
        horizontal_default[i] = default

for i in range(len(vertical_default)):
        vertical_default[i] = default

image_cut_intermediate = cut_image(image_align1_2_3_4)  # not necessary
final_mosaic = final_cut(image_cut_intermediate)  # not necessary

