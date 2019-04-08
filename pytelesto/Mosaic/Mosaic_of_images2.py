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
#import astropy.io as astro
#import glob
import matplotlib.pyplot as plt
#import astropy.io.fits as pyfits
import astroalign as aa

#from astropy.visualization import astropy_mpl_style
#from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
#from astropy.wcs import WCS
#from astropy.io import ascii
#from astropy.table import Table, Column, MaskedColumn

##### Data importation ###############################################

# import image_combination.py
"""
Idealy, import image_combination into this program so everything run in
a single time, but compilation time of image_combination being big, we
directly import the final combined images into the program.
"""

image_data1 = fits.open('.\February 12 2019\patch1.00009674.UCAC3 191_45351.fit')[0].data
image_data2 = fits.open('.\February 12 2019\patch2.00009686.UCAC3 191_44161.fit')[0].data
image_data3 = fits.open('.\February 12 2019\patch3.00009698.Tycho 154_318.fit')[0].data
image_data4 = fits.open('.\February 12 2019\patch4.00009711.Tycho 154_1116.fit')[0].data

##### Functions ##################################################

disp_func = lambda x: np.log10(x)


def save_image(image, name):
    """Save image in a new fits file to be open with ds9

    Do not return anything as it is only succession of actions.

    Args:
        image (array-like): image that we want to save in a new
            file
        name (string): name of the future file
    """
    if save_all_steps:
        name = name + '.fit'
        hdu = fits.PrimaryHDU(image)
        hdul = fits.HDUList([hdu])
        hdul.writeto(name, overwrite=True)


def fusion_of_2_images(im1, im2, rows, columns):
    """Fusion (replacement) of a smaller one in a bigger one

    Return the bigger image with a defined part replaced by the
    smaller one.

    Args:
        im1 (array-like): the biggest image

        im2 (array-like): the smallest image

        rows (range): the row-range where im2 should corresponds
            to im1.

        columns (range)the column-range where im2 should corresponds
            to im1.

    Return:
        im1 (array-like): the bigger image fusioned with the smaller
            one.
    """

    im1[rows[0]:rows[len(rows)-1]+1, columns[0]:columns[len(columns)-1]+1] = im2
    return im1


def cut_distorsion_of_boarder(image, h_boarder, v_boarder):
    image[:, 0:cut_boarder] = h_boarder
    image[:, N[1] - cut_boarder:N[1]] = h_boarder
    image[0:cut_boarder, :] = v_boarder
    image[(N[0] - cut_boarder):N[0], :] = v_boarder
    return image


def adding_frame_to_image(image, place):
    """Add pixels set to defaults around initial image

    Return a masked array of the new image (rectangle described below),
    placed in a much bigger rectangle that have approximately the shape
    of the final mosaic; and return also a mask of the big frame but with
    the area corresponding to the initial image a little smaller (useful
    for after). (Option: save the new image with save_image).

    Args:
        image (array-like): initial image that we want to align
        place (string) : 'top_right' , 'top_left', 'bottom_right'
           'bottom_left', depending on where should approximately
            placed the image in the mosaic

    Returns:
        new_image (MaskedArray) : initial image with bigger size, all new
            pixels being set to default. Mask set to False where the new
            pixels have been added.

        smaller_mask (array-like): mask of 0 and 1, the 1 correspond to
            the area of the image exept a small portion (crop variable)
            on the boarder. Wil be used to remove distorsions later on.

    """
    image_with_frame = np.copy(frame)

    mask_of_normal_size = np.copy(mask_normal_size)
    mask_with_frame = np.copy(frame_mask)

    smaller_mask = mask_with_frame * 0

    # Cut the distorsion on the boarder of the initial image:
    image = cut_distorsion_of_boarder(image, horizontal_boarder*default, vertical_boarder*default)
    save_image(image+mask_of_normal_size, 'test')
    mask_of_normal_size = cut_distorsion_of_boarder(mask_of_normal_size, horizontal_boarder, vertical_boarder)
    save_image(image+mask_of_normal_size, 'test_2')

    if place == 'top_right':
        rows = range(N[0] - 1 + int(add_frame / 2), 2 * N[0] - 1 + int(add_frame / 2))
        columns = range(N[1] - 1 + int(add_frame / 2), 2 * N[1] - 1 + int(add_frame / 2))

        image_with_frame = fusion_of_2_images(image_with_frame, image, rows, columns)

        mask_with_frame = fusion_of_2_images(mask_with_frame, mask_of_normal_size, rows, columns)

    if place == 'top_left':
        rows = range(N[0] - 1 + int(add_frame / 2), 2 * N[0] - 1 + int(add_frame / 2))
        columns = range(int(add_frame / 2) - 1, N[1] - 1 + int(add_frame / 2))

        image_with_frame = fusion_of_2_images(image_with_frame, image, rows, columns)

        mask_with_frame = fusion_of_2_images(mask_with_frame, mask_of_normal_size, rows, columns)

        smaller_mask[(N[0] - 1 + int(add_frame/2) + crop):(2 * N[0] - 1 + int(add_frame/2)-crop),
                     int(add_frame/2) - 1 + crop:(N[1] - 1 + int(add_frame/2) - crop)] = mask_of_ones_smaller_size

    if place == 'bottom_left':
        rows = range(int(add_frame / 2) - 1, N[0] - 1 + int(add_frame / 2))
        columns = range(int(add_frame / 2) - 1, N[1] - 1 + int(add_frame / 2))

        image_with_frame = fusion_of_2_images(image_with_frame, image, rows, columns)

        mask_with_frame = fusion_of_2_images(mask_with_frame, mask_of_normal_size, rows, columns)

        smaller_mask[int(add_frame / 2) - 1 + crop:(N[0] - 1 + int(add_frame / 2) - crop),
                     int(add_frame / 2) - 1 + crop:(N[1] - 1 + int(add_frame / 2) - crop)] = mask_of_ones_smaller_size

    if place == 'bottom_right':
        rows = range(int(add_frame / 2) - 1, N[0] - 1 + int(add_frame / 2))
        columns = range(N[1] - 1 + int(add_frame / 2), 2 * N[1] - 1 + int(add_frame / 2))

        image_with_frame = fusion_of_2_images(image_with_frame, image, rows, columns)

        mask_with_frame = fusion_of_2_images(mask_with_frame, mask_of_normal_size, rows, columns)

        smaller_mask[int(add_frame / 2) - 1 + crop:(N[0] - 1 + int(add_frame / 2)-crop),
                     N[1] - 1 + int(add_frame / 2) + crop:(2 * N[1] - 1 + int(add_frame / 2) - crop)] = mask_of_ones_smaller_size

    place = 'image_with_frame_' + place
    place_of_mask = place + '_mask'
    save_image(image_with_frame, place)
    save_image(mask_with_frame, place_of_mask)

    new_image = np.ma.array(image_with_frame, mask=mask_with_frame)

    return new_image, smaller_mask


def combine_image_with_mask(image, mask):
    """

    """
    common_area = np.where(mask == 1)
    new_im = np.copy(frame)
    new_im[common_area] = image[common_area]

    new_mask = frame_mask
    new_mask[common_area] = mask[common_area]*0
    # im = image.data
    # M = im.shape
    # for i in range(M[0]):
    #     for j in range(M[1]):
    #         if mask[i, j] == 0 and im[i, j] != default:
    #             im[i, j] = default
    # save_image(im-image, 'test')
    new_im = np.ma.array(new_im, mask=new_mask)
    return new_im


def alignment(source, target, smaller_mask, name):
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
    transf, (source_list, target_list) = aa.find_transform(source, target)
    # smaller_mask = aa.apply_transform(transf, smaller_mask, target)
    # save_image(smaller_mask, 'smaller_mask_' + name)

    aligned_image = aa.apply_transform(transf, source, target)

    save_image(aligned_image, 'test4' + name)
    #image = combine_image_with_mask(aligned_image.data, smaller_mask)

    #negative_pixels = np.where(image < 0)
    #image[negative_pixels] = image[negative_pixels]*default
    #save_image(image[negative_pixels], 'baba' + name)
    #print(image[negative_pixels].shape)
    #save_image(image.data, 'test' + name)

    mask_image = aa.apply_transform(transf, source.mask, target)

    save_image(mask_image, 'aloha' + name)
    dark_area = np.where(mask_image == 0)


    # print(dark_area)
    # new_x = []
    # new_y = []
    # n = dark_area[0][0]
    # i = 1
    # first_line = dark_area[0][0]
    # last_line = dark_area[0][-1]
    # while i < (len(dark_area[0])-1):
    #     if dark_area[0][i-1] == dark_area[0][i] & dark_area[0][i+1] == dark_area[0][i]: # & dark_area[0][i] != last_line:
    #         new_x.append(dark_area[0][i])
    #         new_y.append(dark_area[1][i])
    #
    #     i += 1
    #
    # im = image[(new_x, new_y)]

    #for x in dark_area[0]:
        #while x = n:



    # first_x = dark_area[0]
    # x_index = set(dark_area[0])
    # y_index = set(dark_area[1])
    # width = len(x_index)
    # hight = len(y_index)
    # im = np.zeros((width, hight), dtype='f')
    # im = np.reshape(image[dark_area], (width, hight))
    #save_image(im, 'mmm' + name)
    new = np.copy(frame_mask)*0
    new[dark_area] = aligned_image[dark_area]  #image
    #new[(new_x, new_y)] = im
    save_image(new, 'new' + name)

    # im[0, :] = im[0, :] * 0
    # im[M[0], :] = im[M[0], :] * 0
    # im[:, 0] = im[:, 0] * 0
    # im[:, M[1]] = im[:, M[1]] * 0
    # save_image(im, 'mmm1' + name)
    # image[dark_area] = im
    #save_image(im, 'mmm2' + name)

    #new_image = np.maximum(image.data, target.data)
    #save_image(target.data, 'essai1' + name)
    new_image = new + target.data #image.date
    #save_image(new_image, 'essai'+ name)
    common_area = np.where((new > default) & (target.data > default)) #image.data
    #print(common_area)
    new_image[common_area] = new_image[common_area]/2.#target.data[common_area]
    save_image(new_image, 'test2')

    mask_source = np.where(mask_image == 0)
    mask_target = np.where(target.mask == 0)
    new_image_mask = np.copy(frame_mask)
    new_image_mask[mask_source] = mask_image[mask_source]
    new_image_mask[mask_target] = target.mask[mask_target]
    #new_image_mask = image.mask & target.mask

    mosaic = np.ma.array(new_image, mask=new_image_mask)

    name = 'image_align' + name
    save_image(new_image, name)

    return mosaic


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

    # default_area_vertical = np.where(image == vertical_default)
    # print(default_area_vertical)
    # default_area_horizontal = np.where(image == horizontal_default)
    #
    # while i < N[0]:
    #     i = i + 1
    #
    # j = M[0]-1
    # while (image[j, :] == vertical_default).all():
    #     j = j - 1
    #
    # k = 0
    # while (image[:, k] == horizontal_default).all():
    #     k = k + 1
    #
    # l = M[1]-1
    # while (image[:, l] == horizontal_default).all():
    #     l = l - 1
    #
    # image_cut = np.copy(image[i:j, k:l])



    #image_cut = image[default_area_vertical, default_area_horizontal]

    while np.all(image[i, :] == vertical_default):
        i = i + 1

    j = M[0]-1
    while np.all(image[j, :] == vertical_default):
        j = j - 1

    k = 0
    while np.all(image[:, k] == horizontal_default):
        k = k + 1

    m = M[1]-1
    while np.all(image[:, m] == horizontal_default):
        m = m - 1

    image_cut = np.copy(image[i:j, k:m])
    save_image(image_cut, 'image_cut')
    return image_cut


# def final_cut(image):
#     """Final purification of the image with all default values
#         removed
#
#     Return the mosaic image with the right boarders, all bands
#     that still contain partly default values are removed.
#
#     Args:
#          image (array-like): mosaic of image that has been
#             firstly cut with function cut image. The only
#             bands that remain are the one partly equal to
#             default.
#
#     Return:
#         image (array-like): Return the final image, a
#             combination of 4 images that superpose partly
#             and that has been reframed correctly.
#
#     """
#
#     # Starting from the bottom-left corner we count the number
#     # of pixels on a diagonal of a square with default values.
#     # We search for the maximal diagonal possible for the square
#     # to be only composed of pixels set to default. Then, we
#     # check the direction of the band to be removed. If its
#     # vertical, we removed it and then redo the same process,
#     # in case there is also an horizontal band to removed.
#     if image[0, 0] == default:
#         index = 1
#         while image[index, index] == default:
#             index = index + 1
#         if image[index + 1, index] != default:
#             image = image[:, index:-1]
#             save_image(image, '1')
#             if image[0, 0] == default:
#                 index = 1
#                 while image[index, index] == default:
#                     index = index + 1
#                 image = image[index:-1, :]
#                 save_image(image, '2')
#         else:
#             image = image[index:-1, :]
#             save_image(image, '3')
#             if image[0, 0] == default:
#                 index = 1
#                 while image[index, index] == default:
#                     index = index + 1
#                 image = image[:, index:-1]
#                 save_image(image, '4')
#
#     if image[image.shape[0]-1, 0] == default:
#         index1 = image.shape[0] - 2
#         index2 = 1
#         while image[index1, index2] == default:
#             index1 = index1 - 1
#             index2 = index2 + 1
#         if image[index1 - 1, index2] != default:
#             image = image[:, index2:-1]
#             save_image(image, '5')
#             if image[image.shape[0]-1, 0] == default:
#                 index1 = image.shape[0] - 2
#                 index2 = 1
#                 while image[index1, index2] == default:
#                     index1 = index1 - 1
#                     index2 = index2 + 1
#                 image = image[0:index1, :]
#                 save_image(image, '6')
#
#         else:
#             image = image[0:index1, :]
#             save_image(image, '7')
#             if image[image.shape[0]-1, 0] == default:
#                 index1 = image.shape[0] - 2
#                 index2 = 1
#                 while image[index1, index2] == default:
#                     index1 = index1 - 1
#                     index2 = index2 + 1
#                 image = image[:, index2:-1]
#                 save_image(image, '8')
#
#     if image[0, image.shape[1]-1] == default:
#         index1 = 1
#         index2 = image.shape[1] - 2
#         while image[index1, index2] == default:
#             index1 = index1 + 1
#             index2 = index2 - 1
#         if image[index1 + 1, index2] != default:
#             image = image[:, 0:index2]
#             save_image(image, '9')
#             if image[0, image.shape[1]-1] == default:
#                 index1 = 1
#                 index2 = image.shape[1] - 2
#                 while image[index1, index2] == default:
#                     index1 = index1 + 1
#                     index2 = index2 - 1
#                 image = image[index1:-1, :]
#                 save_image(image, '10')
#         else:
#             image = image[index1:-1 :]
#             save_image(image, '11')
#             if image[image.shape[0]-1, 0] == default:
#                 index1 = 1
#                 index2 = image.shape[1] - 2
#                 while image[index1, index2] == default:
#                     index1 = index1 + 1
#                     index2 = index2 - 1
#                 image = image[:, 0:index2]
#                 save_image(image, '12')
#
#     if image[image.shape[0]-1, image.shape[1]-1] == default:
#         index1 = image.shape[0] - 2
#         index2 = image.shape[1] - 2
#         while image[index1, index2] == default:
#             index1 = index1 - 1
#             index2 = index2 - 1
#         if image[index1 - 1, index2] != default:
#             image = image[:, 0:index2]
#             save_image(image, '13')
#             if image[image.shape[0]-1, image.shape[1]-1] == default:
#                 index1 = image.shape[0] - 2
#                 index2 = image.shape[1] - 2
#                 while image[index1, index2] == default:
#                     index1 = index1 - 1
#                     index2 = index2 - 1
#                 image = image[0:index1, :]
#                 save_image(image, '14')
#         else:
#             image = image[0:index1, :]
#             save_image(image, '15')
#             if image[image.shape[0]-1, image.shape[1]-1] == default:
#                 index1 = image.shape[0] - 2
#                 index2 = image.shape[1] - 2
#                 while image[index1, index2] == default:
#                     index1 = index1 - 1
#                     index2 = index2 - 1
#                 image = image[:, 0: index2]
#                 save_image(image, '16')
#
#         save_image(image, 'mosaic')
#
#     return image

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
save_all_steps = True
N = image_data1.shape  # Suppose all images have the same shape
print(N)
add_frame = 400
crop = 20
default = 0  # amélioration: idée, pour que ce ne soit pas vu comme du bruit, utiliser la moyenne des valeurs des pixels
                # de l'image comme défault
cut_boarder = 50  # The boarder of initial images are too noisy, so we cut them a bit to decrease the noise bands in
                  # between all images when thay will be aligned.

frame = np.ones((2 * N[0] + add_frame, 2 * N[1] + add_frame), dtype='f')*default  # image_with_bigger_size

frame_mask = np.ones((2 * N[0] + add_frame, 2 * N[1] + add_frame), dtype='int')  # image_with_bigger_size_mask

mask_normal_size = np.zeros((N[0], N[1]), dtype='int')  # image normal size

mask_of_ones_smaller_size = np.ones((N[0]-crop*2, N[1]-crop*2), dtype='int')

# To cut the boarder of the initial images, replace it with a rectangular boarder with all pixels
# set to default:
horizontal_boarder = np.ones((N[0], cut_boarder), dtype='f')
vertical_boarder = np.ones((cut_boarder, N[1]), dtype='f')

image_top_right, smaller_mask_top_right = adding_frame_to_image(image_data1, 'top_right')
image_top_left, smaller_mask_top_left = adding_frame_to_image(image_data2, 'top_left')
image_bottom_left, smaller_mask_bottom_left = adding_frame_to_image(image_data3, 'bottom_left')
image_bottom_right, smaller_mask_bottom_right = adding_frame_to_image(image_data4, 'bottom_right')


fig, axes = plt.subplots(2, 2, figsize=(10, 10))
ax = axes[0, 0]
ax.imshow(disp_func(image_top_left), origin='lower', cmap='gist_stern')
ax = axes[0, 1]
ax.imshow(disp_func(image_top_right), origin='lower', cmap='gist_stern')
ax = axes[1, 0]
ax.imshow(disp_func(image_bottom_left), origin='lower', cmap='gist_stern')
ax = axes[1, 1]
ax.imshow(disp_func(image_bottom_right), origin='lower', cmap='gist_stern')

plt.show()

combined_image_of_tops = alignment(image_top_left, image_top_right,smaller_mask_top_left, '_tops')
combined_image_of_tops_and_bottom = alignment(image_bottom_left, combined_image_of_tops, smaller_mask_bottom_left, '_tops_and_bottom')
combined_image_of_tops_and_bottoms = alignment(image_bottom_right, combined_image_of_tops_and_bottom, smaller_mask_bottom_right, '_tops_and_bottoms')

fig, axes = plt.subplots(2, 2, figsize=(10, 10))
ax = axes[0, 0]
ax.imshow(disp_func(combined_image_of_tops), origin='lower', cmap='gist_stern')
ax = axes[0, 1]
ax.imshow(disp_func(combined_image_of_tops_and_bottom), origin='lower', cmap='gist_stern')
ax = axes[1, 0]
ax.imshow(disp_func(combined_image_of_tops_and_bottoms), origin='lower', cmap='gist_stern')

plt.show()

horizontal_default = np.zeros((N[0], 1), dtype='f') * default
vertical_default = np.zeros((1, N[1]), dtype='f') * default

final_mosaic = cut_image(combined_image_of_tops_and_bottoms.data)
#final_mosaic = final_cut(image_cut_intermediate)  # not necessary

