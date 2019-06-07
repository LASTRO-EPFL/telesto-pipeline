"""
MOSAIC_FUNCTIONS: functions used for the alignment of the images of
the mosaic.
"""

from Mosaic_basic_function import save_image
import astroalign as aa
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import cv2
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

    place_of_image = '1_image_with_frame_' + place
    place_of_mask = '1_mask_with_frame_' + place
    save_image(image_with_frame, place_of_image, save)
    save_image(mask_with_frame, place_of_mask, save)

    new_image = np.ma.array(image_with_frame, mask=mask_with_frame)

    return new_image


def circle_sources(aligned_source, source, target, p, pos_source, pos_target, name, save):
    """
    To circle the matching stars used to find the transformation. This function
    can be suppressed from from the pipelines, it was just useful for the report.
    """
    map = 'gist_ncar'
    print(len(pos_target))
    colors = ['r', 'g', 'b', 'floralwhite', 'cyan', 'w', 'm', 'chartreuse', 'lightpink', 'fuchsia', 'salmon', 'olive',
              'chocolate', 'darkblue', 'darkred', 'darkviolet', 'royalblue']

    if len(colors) > len(pos_target):
        n = len(pos_target)
    else:
        n = len(colors)

    source_bis = disp_func(source)
    target_bis = disp_func(target)
    aligned_source_bis = disp_func(aligned_source)

    fig, axes = plt.subplots()
    axes.imshow(source_bis, cmap=map, interpolation='none', origin='lower')
    axes.axis('off')
    axes.set_title("Source Image")
    for (xp, yp), c in zip(pos_source[:n], colors):
        circ = plt.Circle((xp, yp), 8, fill=False, edgecolor=c, linewidth=2)
        axes.add_patch(circ)
        source_bis = cv2.circle(source_bis, (int(xp), int(yp)), radius=int(10), color=(0, 255, 0),
                                thickness=4)
    plt.show()

    fig, axes = plt.subplots()
    axes.imshow(target_bis, cmap=map, interpolation='none', origin='lower')
    axes.axis('off')
    axes.set_title("Target Image")
    for (xp, yp), c in zip(pos_target[:n], colors):
        circ = plt.Circle((xp, yp), 8 * p.scale, fill=False, edgecolor=c, linewidth=2)
        axes.add_patch(circ)
        target_bis = cv2.circle(target_bis, (int(xp), int(yp)), radius=int(10 * p.scale), color=(0, 255, 0), thickness=4)
    plt.show()

    fig, axes = plt.subplots()
    axes.imshow(aligned_source_bis, cmap=map, interpolation='none', origin='lower')
    axes.axis('off')
    axes.set_title("Source Image aligned with Target")
    for (xp, yp), c in zip(pos_target[:n], colors):
        circ = plt.Circle((xp, yp), 8 * p.scale, fill=False, edgecolor=c, linewidth=2)
        axes.add_patch(circ)
        aligned_source_bis = cv2.circle(aligned_source_bis, (int(xp), int(yp)), radius=int(10 * p.scale), color=(
                                        0, 255, 0), thickness=4)

    plt.show()

    mpimg.imsave('2_bis_source' + name, source_bis, cmap='gist_ncar', origin='lower')
    mpimg.imsave('2_bis_target' + name,  target_bis, cmap='gist_ncar', origin='lower')
    mpimg.imsave('2_bis_aligned_source' + name,  aligned_source_bis, cmap='gist_ncar', origin='lower')


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
    transf, (pos_source, pos_target) = aa.find_transform(source, target)
    aligned_image = aa.apply_transform(transf, source, target)
    #circle_sources(aligned_image, source.data, target.data, transf, pos_source, pos_target, name, save)   # can be
    # suppressed from the code

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


# def alignment(source, target, name, frame, save):
#     """ Aligns 2 images at a time with method 2 (see report).
        # This func can be suppressed from the code.
#
#     Return the combination of source and target images
#     aligned with each other, presented as a Masked Array.
#     Firstly, aligns source on target. Then re-add target
#     image on align image. Keeps only the pixels of the a
#     aligned image that have a false value (0) in the
#     associated mask.
#     (Option: save the new images with save_image).
#
#     Args:
#         source (masked array): A masked array of the source image to be
#             transformed, its pixels will match the ones of the target
#             image
#
#         target (masked array): A mask array of the target (destination)
#             image.
#
#         name (string): Name of the future file, that is a saving of the
#             final image.
#
#         frame (array-like): big rectangle of 0, frame of the initial image.
#
#         save (bool): option to activate if you want to save all the steps
#            of the process.
#
#     Returns:
#         mosaic (masked array) : A mask array of the new image, combination of
#             the source and target image.
#
#     """
#     transf, (pos_source, pos_target) = aa.find_transform(source, target)
#     aligned_image = aa.apply_transform(transf, source, target)
#     circle_sources(aligned_image, source.data, target.data, transf, pos_source, pos_target, name, save)
#
#     mask_aligned_image = aa.apply_transform(transf, source.mask, target)
#
#     save_image(aligned_image, '2_aligned_image' + name, save)
#     save_image(mask_aligned_image, '2_mask_aligned_image' + name, save)
#
#     image_area = np.where(mask_aligned_image == 0)
#
#     new_aligned_image = np.copy(frame)
#     new_aligned_image[image_area] = aligned_image[image_area]
#
#     save_image(new_aligned_image, '3_aligned_image_purified' + name, save)
#
#     mosaic = np.ma.array(new_aligned_image, mask=mask_aligned_image)
#
#     return mosaic


def combined_images(im1, im2, name, save):
    """
    To combine images, using method 2 of alignment (see report).
    This method can be suppressed from the code also, so it is
    not detailed.
    """
    new_aligned_image_combined = im1.data + im2.data
    common_area = np.where((im1.data > 0) & (im2.data > 0))

    new_aligned_image_combined[common_area] = new_aligned_image_combined[common_area]/2.  # np.maximum

    save_image(new_aligned_image_combined, '4_aligned_image_combined' + name, save)

    mask_target = np.where(im2.mask == 0)
    im1.mask[mask_target] = 0

    mosaic = np.ma.array(new_aligned_image_combined, mask=im1.mask)

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


def histogram(image, save):
    """
    Correct the range of luminosity to reduce the noise on the image.

    Return the image which interval of luminosity is reduced. Setting
    the values outside the interval to the value of  the closest interval's
    boundary.

    Args:
        image (array-like): final mosaic to correct in luminosity.

        save (bool): option to activate if you want to save all the steps
           of the process.

    Return:
        image_corrected (array-like): final mosaic corrected.
    """
    # The main values of the interesting luminosity are in this range:
    area_covered = np.where(image < (1./5.)*image.max())

    hist, bins, _ = plt.hist(image[area_covered], bins=7000)

    plt.figure()
    plt.plot(np.log(hist))  # log scale
    plt.xlabel('Flux')
    plt.grid(True)
    plt.show()

    maximum = np.where(hist == max(hist))
    print('You can rescale the luminosity range of your image to [max(hist)-a; max(hist)+b], choosing a and b. ')
    a = input('What lower distance a from the maximum of the histogram do you want? (int)')
    a = int(a)
    b = input('What upper distance b from the maximum of the histogram do you want? (int)')
    b = int(b)
    inter = [bins[maximum] - a, bins[maximum] + b]

    image_corrected = np.copy(image)
    air1 = np.where(image < inter[0])
    air2 = np.where(image > inter[1])
    image_corrected[air1] = inter[0]
    image_corrected[air2] = inter[1]

    fig, axes = plt.subplots(2, 1)
    ax = axes[0]
    plt.title('Image not corrected')
    ax.imshow(disp_func(image), origin='lower', cmap='gist_ncar')
    ax = axes[1]
    plt.title('Image corrected')
    ax.imshow(disp_func(image_corrected), origin='lower', cmap='gist_ncar')
    plt.show()

    save_image(image_corrected, '7_final_mosaic_corrected_in_luminosity', save)

    return image_corrected