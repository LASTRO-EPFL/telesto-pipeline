"""
IMAGE_COMBINATION_FUNCTIONS: functions used for the combination
and the alignment of the patches.
"""

import astropy.io as astro
import glob
import astroalign as aa
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import numpy as np
from Mosaic_basic_function import save_image


def noise_calculation(image, k, save_steps):
    """
    Calculate the rms noise of an image, based on
    the pixels values in a small square of the image.

    Return the root mean squared value of the small
    square in the image.

    Args:
        image (array-like): image of which we want to
            measure the noise.

        k (int): number of the image (1, 2, 3 or 4)
            depending on what part of the mosaic the
            image corresponds to (top, bottom, right,
            left).

        save_steps (string): True if you want to save
            transitory steps.

    Return:
        rms (double): root mean squared value of  a
            small square in the image.
    """
    N = image.shape

    # Size of the small square:
    l = 45

    # We take the square on the corner of the image
    # the most distant from the center of the nebula:
    if k == 1:
        square = image[N[0]-l: N[0], N[1]-l: N[1]]

    if k == 2:
        square = image[N[1] - l: N[1], 0: l]

    if k == 3:
        square = image[0: l, 0: l]

    if k == 4:
        square = image[0: l, N[1] - l: N[1]]

    neg_value = np.where(square < 0)
    square[neg_value] = 0

    rms = np.sqrt(np.mean(square**2))

    save_image(square, 'square%d' % k, save_steps)

    return rms #np.mean(square)/rms


def open_all_images(k, save_steps):
    """
    Open all patches of a same portion of the mosaic and
    calculate the noise on these images.

    Return a dictionary gathering all the patches i (i=1,2,3,4),
    a list gathering the noise corresponding to each patch and the
    number of images in the dictionary.

    Args:
        k (int): number of the image (1, 2, 3 or 4)
            depending on what part of the mosaic the
            image corresponds to (top, bottom, right,
            left).

        save_steps (string): True if you want to save
           transitory steps.

    Return:
        image_data (dictionary): each key corresponds to
            one patch.

        noise (list of double): nlist of the noise for each
            patch in the dictionary.

        len(name_all_images) (int): the number of patches in
            the dictionary.
    """
    name_all_images = glob.glob('.\Reduced_patch\patch%d*.fit' % k)

    noise = []
    image_data = {}
    for i in range(0, len(name_all_images) - 1):
        image_file = get_pkg_data_filename(name_all_images[i])
        image_data['image%d' % i] = astro.fits.getdata(image_file, ext=0)
        noise.append(noise_calculation(image_data['image%d' % i], k, save_steps))

    return image_data, noise, len(name_all_images)


def image_combination_choosing_patch(k, save_steps):
    """
    Aligned the patches and combined them trying to optimize
    the noise.

    Return a mean image of the patches chosen in the optimization
    of the noise.

    Args:
         k (int): number of the image (1, 2, 3 or 4)
            depending on what part of the mosaic the
            image corresponds to (top, bottom, right,
            left).

        save_steps (string): True if you want to save
           transitory steps.

    Return:
        image_mid (array-like): mean of the chosen patches
            previously aligned.
    """
    image_data, noise, n = open_all_images(k, save_steps)

    # We align w.r.t the image with the smallest noise.
    index = noise.index(min(noise))
    min_noise = noise[index]

    image_data0 = image_data['image%d' % index]
    image_mid = image_data0

    for i in range(0, n-1):
        if i != index:
            image_aligned = aa.register(image_data['image%d' % i], image_data0)
            new_image_mid = (image_mid + image_aligned)/2

            new_noise = noise_calculation(new_image_mid, k, save_steps)

            # We keep the combination of the image with the new_patch only
            # if the noise of the later improve the global quality:
            if new_noise < min_noise:
                image_mid = new_image_mid
                min_noise = new_noise
            else:
                print('We did not use patch',  i+1, 'for image', k)

    save_image(image_mid, 'patch%d_aligned_and_combined' % k, save_steps)

    return image_mid


def image_combination_using_all_patch(k, save_steps):
    """
    Aligned all patches and combined them.

    Return a mean image of the patches.

    Args:
         k (int): number of the image (1, 2, 3 or 4)
            depending on what part of the mosaic the
            image corresponds to (top, bottom, right,
            left).

        save_steps (string): True if you want to save
           transitory steps.

    Return:
        image_mid (array-like): mean of all patches
            previously aligned.
    """
    image_data, noise, n = open_all_images(k, save_steps)

    image_data0 = image_data['image0']
    image_mid = image_data0

    for i in range(1, n-1):
        image_aligned = aa.register(image_data['image%d' %i], image_data0)
        image_mid = (image_mid + image_aligned)/2

    save_image(image_mid, 'patch%d_aligned_and_combined_using_all_patch' % k, save_steps)

    return image_mid

