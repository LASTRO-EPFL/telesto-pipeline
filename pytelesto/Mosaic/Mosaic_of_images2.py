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
import matplotlib.pyplot as plt
import astroalign as aa
from astropy.io import fits

import sep
sep.set_extract_pixstack(1e6)


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

    im1[rows[0]:rows[-1]+1, columns[0]:columns[-1]+1] = im2
    return im1


def cut_distorsion_of_boarder(image, h_boarder, v_boarder):
    """Cut all sides of the image by a band of width 'cut_boarder'

    Return the same image, cropped. The pixels 'removed are in reality
    set to 0, which is the value of the pixels in the added frame of
    the image.

    Args:
        image (array-like): the image we want to crop.

        h_boarder (array-like): the horizontal vector of 0 to replaced
            the horizontal sides of the image (i.e crop them).

        v_boarder (array-like): the vertical vector of 0 to replaced
            the vertical sides of the image (i.e crop them).

    Return:
        image (array-like): the image cropped.

    """
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

    # Cut the distorsion on the boarder of the initial image:
    image = cut_distorsion_of_boarder(image, 0, 0)
    mask_of_normal_size = cut_distorsion_of_boarder(mask_of_normal_size, 1, 1)

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

    if place == 'bottom_left':
        rows = range(int(add_frame / 2) - 1, N[0] - 1 + int(add_frame / 2))
        columns = range(int(add_frame / 2) - 1, N[1] - 1 + int(add_frame / 2))

        image_with_frame = fusion_of_2_images(image_with_frame, image, rows, columns)

        mask_with_frame = fusion_of_2_images(mask_with_frame, mask_of_normal_size, rows, columns)

    if place == 'bottom_right':
        rows = range(int(add_frame / 2) - 1, N[0] - 1 + int(add_frame / 2))
        columns = range(N[1] - 1 + int(add_frame / 2), 2 * N[1] - 1 + int(add_frame / 2))

        image_with_frame = fusion_of_2_images(image_with_frame, image, rows, columns)

        mask_with_frame = fusion_of_2_images(mask_with_frame, mask_of_normal_size, rows, columns)

    place = '1_image_with_frame_' + place
    place_of_mask = '1_mask_with_frame_' + place
    save_image(image_with_frame, place)
    save_image(mask_with_frame, place_of_mask)

    new_image = np.ma.array(image_with_frame, mask=mask_with_frame)

    return new_image


def alignment(source, target, name):
    """ Aligns 2 images at a time.

    Return the combination of source and target images
    aligned with each other. Firstly, aligns source on
    target. Then re-add target image on align image.
    Keeps only the pixels of the aligned image that
    have a false value (0) in the associated mask.
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
        mosaic (masked array) : A mask array of the new image, combination of
            the source and target image.

    """
    transf, (source_list, target_list) = aa.find_transform(source, target)

    aligned_image = aa.apply_transform(transf, source, target)
    mask_aligned_image = aa.apply_transform(transf, source.mask, target)

    save_image(aligned_image, '2_aligned_image' + name)
    save_image(mask_aligned_image, '2_mask_aligned_image' + name)

    image_area = np.where(mask_aligned_image == 0)

    new_aligned_image = np.copy(frame_mask)*0
    new_aligned_image[image_area] = aligned_image[image_area]

    save_image(new_aligned_image, '3_aligned_image_purified' + name)

    new_aligned_image_combined = new_aligned_image + target.data
    save_image(new_aligned_image_combined, '4_test' +name)
    common_area = np.where((new_aligned_image > 0) & (target.data > 0))

    new_aligned_image_combined[common_area] = new_aligned_image_combined[common_area]/2.

    save_image(new_aligned_image_combined, '4_aligned_image_combined' + name)

    mask_source = np.where(mask_aligned_image == 0)
    mask_target = np.where(target.mask == 0)
    new_mask_aligned_image = np.copy(frame_mask)
    new_mask_aligned_image[mask_source] = mask_aligned_image[mask_source]
    new_mask_aligned_image[mask_target] = target.mask[mask_target]

    save_image(new_mask_aligned_image, '4_mask_aligned_image_combined' + name)

    mosaic = np.ma.array(new_aligned_image_combined, mask=new_mask_aligned_image)

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
    # Starting from the left boarder, we count the number of
    # vertical line that have only default values in order to
    # remove them later. Same thing done starting from the right
    # boarder, top and bottom boarder.


    # horizontal_dark_area = np.where(image != horizontal_default)
    # horizontal_dark_area =
    # print(len(horizontal_dark_area))
    # new_image = np.zeros((len(horizontal_dark_area[0]), len(horizontal_dark_area[1])), dtype='f')
    # new_image[horizontal_dark_area] = image[horizontal_dark_area]
    # save_image(new_image, '5_final_mosaic_test1')
    # vertical_dark_area = np.where(image != vertical_default)
    # new_image2 = np.zeros((len(vertical_dark_area[0]), len(vertical_dark_area[1])), dtype='f')
    # print(vertical_dark_area)
    # new_image2[vertical_dark_area] = new_image[vertical_dark_area]
    # save_image(new_image2, '5_final_mosaic_test2')
    # return new_image2

    i = 0
    save_image(image[i, :], 'aloha')
    save_image(image[:, 0], 'aloha2')
    save_image(horizontal_default, 'h')
    save_image(vertical_default, 'v')
    while (i < (M[0]/3.)) & np.all(image[i, :] == horizontal_default):
        i = i + 1
    print('a')
    j = M[0]-1
    while (j > (M[0]/3.)) & np.all(image[j, :] == horizontal_default):
        j = j - 1
    print('b')
    k = 0
    while (k < (M[1]/3.)) & np.all(image[:, k] == vertical_default):
        k = k + 1
    print('c')
    m = M[1]-1
    while (m > (M[1]/3.)) & np.all(image[:, m] == vertical_default):
        m = m - 1
    print('d')
    image_cut = np.copy(image[i:j, k:m])
    save_image(image_cut, '5_final_mosaic')
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
save_all_steps = input('Do you want to save step by step  the treatment of the images in a file? (Yes/No) ')
if save_all_steps == 'Yes' or save_all_steps == 'yes' or save_all_steps == 'y' or save_all_steps == 'Y':
    save_all_steps = True
elif save_all_steps == 'No' or save_all_steps == 'no' or save_all_steps == 'n' or save_all_steps == 'N':
    save_all_steps = False
else:
    raise Exception('Invalid entry, you should enter: Yes or No')

N = image_data1.shape  # Suppose all images have the same shape
cut_boarder = 50  # The boarder of initial images are too noisy, so we cut them a bit to decrease the noise bands in
                  # between all images when thay will be aligned.
add_frame = 400

frame = np.zeros((2 * N[0] + add_frame, 2 * N[1] + add_frame), dtype='f') # image_with_bigger_size

frame_mask = np.ones((2 * N[0] + add_frame, 2 * N[1] + add_frame), dtype='int')  # image_with_bigger_size_mask

mask_normal_size = np.zeros((N[0], N[1]), dtype='int')  # image normal size

# To cut the boarder of the initial images, replace it with a rectangular boarder with all pixels
# set to default:
horizontal_boarder = np.ones((N[0], cut_boarder), dtype='f')
vertical_boarder = np.ones((cut_boarder, N[1]), dtype='f')

image_top_right = adding_frame_to_image(image_data1, 'top_right')
image_top_left = adding_frame_to_image(image_data2, 'top_left')
image_bottom_left = adding_frame_to_image(image_data3, 'bottom_left')
image_bottom_right = adding_frame_to_image(image_data4, 'bottom_right')

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

combined_image_of_tops = alignment(image_top_left, image_top_right, '_tops')
combined_image_of_tops_and_bottom = alignment(image_bottom_left, combined_image_of_tops, '_tops_and_bottom')
combined_image_of_tops_and_bottoms = alignment(image_bottom_right, combined_image_of_tops_and_bottom, '_tops_and_bottoms')

fig, axes = plt.subplots(2, 2, figsize=(10, 10))
ax = axes[0, 0]
ax.imshow(disp_func(combined_image_of_tops), origin='lower', cmap='flag')
ax = axes[0, 1]
ax.imshow(disp_func(combined_image_of_tops_and_bottom), origin='lower', cmap='flag')
ax = axes[1, 0]
ax.imshow(disp_func(combined_image_of_tops_and_bottoms), origin='lower', cmap='flag')
plt.show()

M = combined_image_of_tops_and_bottoms.data.shape
horizontal_default = np.zeros((M[0], 1), dtype='f')
vertical_default = np.zeros((1, M[1]), dtype='f')

final_mosaic = cut_image(combined_image_of_tops_and_bottoms.data)

figure_2 = plt.figure()
plt.imshow(disp_func(final_mosaic), cmap='flag')
#plt.legend()
plt.xlabel('pixels')
plt.ylabel('pixels')
plt.grid(True)
plt.show()

satisfied = False
while not satisfied:
    height = input('What height do you want for your image? (in number of pixels)')
    height = int(height)
    width = input('What width do you want for your image? (in number of pixels)')
    width = int(width)

    final_mosaic_crop = final_mosaic[int(M[0]/2.-height/2.):int(M[0]/2+(height/2.)), int(M[1]/2.-width/2.):int(M[1]/2+(width/2.))]

    figure_3 = plt.figure()
    plt.imshow(disp_func(final_mosaic_crop), cmap='flag')
    plt.xlabel('pixels')
    plt.ylabel('pixels')
    plt.grid(True)
    plt.show()

    satisfied = input('Are you satisfied with the shape of this image? (Yes/No)')
    if satisfied == 'Yes' or satisfied == 'yes' or satisfied == 'y' or satisfied == 'Y':
        satisfied = True
        save_image(final_mosaic_crop, '6_final_mosaic_crop')
    elif satisfied == 'No' or satisfied == 'no' or satisfied == 'n' or satisfied == 'N':
        satisfied = False
    else:
        raise Exception('Invalid entry, tou should enter: Yes or No')

#final_mosaic = final_cut(image_cut_intermediate)  # not necessary

