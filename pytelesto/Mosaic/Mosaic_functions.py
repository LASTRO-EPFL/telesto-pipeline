"""
MOSAIC_FUNCTIONS: functions used for the alignment of the images of
the mosaic.
"""

from Mosaic_basic_function import save_image
import astroalign as aa
import numpy as np
import sep
sep.set_extract_pixstack(1e8)

# To plot the images:
disp_func = lambda x: np.log10(x)


def fusion_of_2_images(im1, im2, rows, columns):
    """Merging (placement) of a smaller image in a bigger one.

    Return the bigger image with a defined part replaced by the
    smaller one.

    Args:
        im1 (array-like): the biggest image

        im2 (array-like): the smallest image

        rows (range): the row-range where im2 should corresponds
            to im1.

        columns (range): the column-range where im2 should corresponds
            to im1.

    Return:
        im1 (array-like): the bigger image merged with the smaller
            one.
    """
    im1[rows[0]:rows[-1]+1, columns[0]:columns[-1]+1] = im2
    return im1


def cut_distorsion_of_boarder(image, index, cut_boarder, N):
    """Cut all sides of the initial image by a band of width 'cut_boarder'

    Return the same image, cropped. The pixels 'removed' are in reality
    set to 0, which is the value of the pixels in the added frame of
    the image.

    Args:
        image (array-like): the image we want to crop.

        index (int): value for the pixels cropped i.e value used in the
            frame.

        cut_boarder (int): width for the band of index numbers used
            to 'cut' the boarders if initial image.

        N (int): size of the initial image.

    Return:
        image (array-like): the image 'cropped'.
    """
    image[:, 0:cut_boarder] = index
    image[:, N[1] - cut_boarder:N[1]] = index
    image[0:cut_boarder, :] = index
    image[(N[0] - cut_boarder):N[0], :] = index
    return image


def adding_frame_to_image(image, place, frame, mask_normal_size,
                          frame_mask, cut_boarder, add_frame, N, w, h, save):
    """Creates a masked array of initial image, placed in a bigger frame of 0.

    Return a masked array of the initial image placed in a much bigger rectangle
    (frame of 0) that have approximately the shape of the final mosaic.
    (Option: save the new images with save).

    Args:
        image (array-like): initial image that we want to align/

        place (string) : 'top_right' , 'top_left', 'bottom_right'
           'bottom_left', depending on where should be approximately
            placed the image in the mosaic

        frame (array-like): big rectangle of 0, where will be placed at
            position 'place' the initial image.

        mask_normal_size (array-like): rectangle of 0 with the same shape
            as initial image.

        frame_mask (array-like): big rectangle of 1 with the same shape as
            'frame'. 'mask_normal_size' will be later placed at position
            in frame_mask.

        cut_boarder (int): width for the band used to 'cut' the boarders
            if initial image.

        add_frame (int): minimal number of pixels we want from the boader
            of initial image to the boarder of the frame.

        N (int): size of the initial image.

        w (double): the width of the frame will be of 'w + add_frame'.

        h (double): the height of the frame will be of 'w + add_frame'.

        save (bool): option to activate if you want to save all the steps
            of the process.

    Returns:
        new_image (MaskedArray) : initial image with bigger size, all new
            pixels being set to 0. Mask is set to False where the new
            pixels have been added.
    """
    image_with_frame = np.copy(frame)
    mask_of_normal_size = np.copy(mask_normal_size)
    mask_with_frame = np.copy(frame_mask)

    # Cuts the distorsion on the boarder of the initial image:
    image = cut_distorsion_of_boarder(image, 0, cut_boarder, N)
    mask_of_normal_size = cut_distorsion_of_boarder(mask_of_normal_size, 1, cut_boarder, N)

    # Positions the initial images in the frame according to their position in the mosaic:
    if place == 'top_right':
        rows = range(int(0.5*N[0]) - 1 + int(add_frame / 2), w - 1 + int(add_frame / 2))
        columns = range(int(0.5*N[1]) - 1 + int(add_frame / 2), h - 1 + int(add_frame / 2))

        image_with_frame = fusion_of_2_images(image_with_frame, image, rows, columns)

        mask_with_frame = fusion_of_2_images(mask_with_frame, mask_of_normal_size, rows, columns)

    if place == 'top_left':
        rows = range(int(0.5*N[0]) - 1 + int(add_frame / 2), w - 1 + int(add_frame / 2))
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
        columns = range(int(0.5*N[1]) - 1 + int(add_frame / 2), h - 1 + int(add_frame / 2))

        image_with_frame = fusion_of_2_images(image_with_frame, image, rows, columns)

        mask_with_frame = fusion_of_2_images(mask_with_frame, mask_of_normal_size, rows, columns)

    place = '1_image_with_frame_' + place
    place_of_mask = '1_mask_with_frame_' + place
    save_image(image_with_frame, place, save)
    save_image(mask_with_frame, place_of_mask, save)

    new_image = np.ma.array(image_with_frame, mask=mask_with_frame)

    return new_image


def alignment(source, target, name, frame, save):
    """ Aligns 2 images at a time.

    Return the combination of source and target images
    aligned with each other, presented as a Masked Array.
    Firstly, aligns source on target. Then re-add target
    image on align image. Keeps only the pixels of the a
    ligned image that have a false value (0) in the
    associated mask.
    (Option: save the new images with save_image).

    Args:
        source (masked array): A masked array of the source image to be
            transformed, its pixels will match the ones of the target
            image

        target (masked array): A nmask array of the target (destination)
            image.

        name (string): Name of the future file, that is a saving of the
            final image.

        frame (array-like): big rectangle of 0, frame of the initial image.

        save (bool): option to activate if you want to save all the steps
           of the process.

    Returns:
        mosaic (masked array) : A mask array of the new image, combination of
            the source and target image.

    """
    transf, (source_list, target_list) = aa.find_transform(source, target)
    aligned_image = aa.apply_transform(transf, source, target)
    mask_aligned_image = aa.apply_transform(transf, source.mask, target)

    save_image(aligned_image, '2_aligned_image' + name, save)
    save_image(mask_aligned_image, '2_mask_aligned_image' + name, save)

    image_area = np.where(mask_aligned_image == 0)

    new_aligned_image = np.copy(frame)
    new_aligned_image[image_area] = aligned_image[image_area]

    save_image(new_aligned_image, '3_aligned_image_purified' + name, save)

    new_aligned_image_combined = new_aligned_image + target.data
    common_area = np.where((new_aligned_image > 0) & (target.data > 0))

    new_aligned_image_combined[common_area] = new_aligned_image_combined[common_area]/2.  # np.maximum

    save_image(new_aligned_image_combined, '4_aligned_image_combined' + name, save)

    mask_target = np.where(target.mask == 0)
    mask_aligned_image[mask_target] = 0

    save_image(mask_aligned_image, '4_mask_aligned_image_combined' + name, save)

    mosaic = np.ma.array(new_aligned_image_combined, mask=mask_aligned_image)

    return mosaic

def cut_image(image, save):
    """
    Cropped the full bands of 0 in the final mosaic.

    Return the final mosaic where useless band of 0
    have been removed.

    Args:
        image (masked array): final mosaic to crop.

        save (bool): option to activate if you want to save all the steps
           of the process.

    Return:
        cut_image (array-like): final mosaic cropped.
    """
    image_area = np.where(image.mask == 0)
    cut_image = image.data[min(image_area[0]): max(image_area[0]), min(image_area[1]): max(image_area[1])]
    save_image(cut_image, '5_test', save)

    return cut_image