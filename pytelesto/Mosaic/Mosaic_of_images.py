"""
MOSAIC_OF_IMAGES : from 4 images that superpose partly, creates a fifth,
bigger, image that is a combination of the previous ones.
Importance of the PLACE of the four images in the mosaic. Name them such
that :
  - image_data1 = top right of the mosaic
  - image_data2 = top left of the mosaic
  - image_data3 = bottom left of the mosaic
  - image_data4 = bottom right of the mosaic
"""

import os
from astropy.io import fits
from pylab import *
from Mosaic_functions import disp_func, adding_frame_to_image, alignment, cut_image
from Mosaic_basic_function import save_image

if len(os.listdir('.\Mosaic')) == 4:
    image_data1 = fits.open('.\Mosaic\image_top_right.fit')[0].data
    image_data2 = fits.open('.\Mosaic\image_top_left.fit')[0].data
    image_data3 = fits.open('.\Mosaic\image_bottom_left.fit')[0].data
    image_data4 = fits.open('.\Mosaic\image_bottom_right.fit')[0].data

    # image_data1 = fits.open('.\Mosaic_all_patch\image_top_right.fit')[0].data
    # image_data2 = fits.open('.\Mosaic_all_patch\image_top_left.fit')[0].data
    # image_data3 = fits.open('.\Mosaic_all_patch\image_bottom_left.fit')[0].data
    # image_data4 = fits.open('.\Mosaic_all_patch\image_bottom_right.fit')[0].data
else:
    raise Exception('Mosaic images are not correctly placed in the Mosaic file. See instructions in README file.')


save_all_steps = input('Do you want to save step by step  the treatment of the images in a file? (Yes/No) ')
if save_all_steps == 'Yes' or save_all_steps == 'yes' or save_all_steps == 'y' or save_all_steps == 'Y':
    save_all_steps = True
elif save_all_steps == 'No' or save_all_steps == 'no' or save_all_steps == 'n' or save_all_steps == 'N':
    save_all_steps = False
else:
    raise Exception('Invalid entry, you should enter: Yes or No')

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

# Suppose all images have the same shape:
N = image_data1.shape

# The boarders of initial images are too noisy,
# thus we cut them a little to avoid noisy bands
# during the alignment of the images. To do that
# we replace a band on each side of the image with
# a band of 1 (= False for the masked array):
cut_boarder = 30
horizontal_boarder = np.ones((N[0], cut_boarder), dtype='f')
vertical_boarder = np.ones((cut_boarder, N[1]), dtype='f')

# Initial images are placed in a frame of width
# "width + add_frame" and height "height + add_frame".
# add_frame controls the number of pixels you want to
# have between the boarders of your image and the ones
# of the frame:
add_frame = 1600
width = int(3/2*N[0])
height = int(3/2*N[1])

# "frame" is a set of 0 on which we will paste the initial
# images (image_datai i=1,2,3,4):
frame = np.zeros((width + add_frame, height + add_frame), dtype='f')

# "frame_mask" is a set of 1. The pixels corresponding to
# the initial image will be later set to 0 using "mask_normal_size":
frame_mask = np.ones((width + add_frame, height + add_frame), dtype='int')
mask_normal_size = np.zeros((N[0], N[1]), dtype='int')


image_top_right = adding_frame_to_image(image_data1, 'top_right', frame, mask_normal_size, frame_mask, cut_boarder,
                                        add_frame, N, width, height, save_all_steps)

image_top_left = adding_frame_to_image(image_data2, 'top_left', frame, mask_normal_size, frame_mask, cut_boarder,
                                       add_frame, N, width, height, save_all_steps)

image_bottom_left = adding_frame_to_image(image_data3, 'bottom_left', frame, mask_normal_size, frame_mask, cut_boarder,
                                          add_frame, N, width, height, save_all_steps)

image_bottom_right = adding_frame_to_image(image_data4, 'bottom_right', frame, mask_normal_size, frame_mask, cut_boarder
                                           , add_frame, N, width, height, save_all_steps)

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

combined_image_of_tops = alignment(image_top_left, image_top_right, '_tops', frame, save_all_steps)

combined_image_of_tops_and_bottom = alignment(image_bottom_left, combined_image_of_tops,
                                              '_tops_and_bottom', frame, save_all_steps)

combined_image_of_tops_and_bottoms = alignment(image_bottom_right, combined_image_of_tops_and_bottom,
                                               '_tops_and_bottoms', frame, save_all_steps)

fig, axes = plt.subplots(2, 2, figsize=(10, 10))
ax = axes[0, 0]
ax.imshow(disp_func(combined_image_of_tops), origin='lower', cmap='flag')
ax = axes[0, 1]
ax.imshow(disp_func(combined_image_of_tops_and_bottom), origin='lower', cmap='flag')
ax = axes[1, 0]
ax.imshow(disp_func(combined_image_of_tops_and_bottoms), origin='lower', cmap='flag')
plt.show()

final_mosaic = cut_image(combined_image_of_tops_and_bottoms, save_all_steps)

fig, axes = plt.subplots()
axes = gca()
axes.imshow(disp_func(final_mosaic), cmap='flag')
starty, endy = axes.get_ylim()
axes.yaxis.set_ticks(np.arange(starty, endy, 100))
axes.xaxis.set_major_locator(MultipleLocator(500))
axes.xaxis.set_minor_locator(MultipleLocator(100))
axes.yaxis.set_major_locator(MultipleLocator(500))
axes.yaxis.set_minor_locator(MultipleLocator(100))
axes.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
axes.grid(which='minor', axis='x', linewidth=0.35, linestyle='-', color='0.75')
axes.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
axes.grid(which='minor', axis='y', linewidth=0.35, linestyle='-', color='0.75')
axes.set_xlabel('pixels')
axes.set_ylabel('pixels')
plt.show()

satisfied = False
while not satisfied:
    height = input('What height do you want for your image? (in number of pixels)')
    height = int(height)
    width = input('What width do you want for your image? (in number of pixels)')
    width = int(width)
    M = final_mosaic.shape

    final_mosaic_crop = final_mosaic[int(M[0]/2.-height/2.):int(M[0]/2+(height/2.)), int(M[1]/2.-width/2.):int(M[1]/2+(width/2.))]

    fig, axes = plt.subplots()
    axes = gca()
    axes.imshow(disp_func(final_mosaic_crop), cmap='flag')
    starty, endy = axes.get_ylim()
    axes.yaxis.set_ticks(np.arange(starty, endy, 100))
    axes.xaxis.set_major_locator(MultipleLocator(500))
    axes.xaxis.set_minor_locator(MultipleLocator(100))
    axes.yaxis.set_major_locator(MultipleLocator(500))
    axes.yaxis.set_minor_locator(MultipleLocator(100))
    axes.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
    axes.grid(which='minor', axis='x', linewidth=0.35, linestyle='-', color='0.75')
    axes.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
    axes.grid(which='minor', axis='y', linewidth=0.35, linestyle='-', color='0.75')
    axes.set_xlabel('pixels')
    axes.set_ylabel('pixels')
    plt.show()

    satisfied = input('Are you satisfied with the shape of this image? (Yes/No)')
    if satisfied == 'Yes' or satisfied == 'yes' or satisfied == 'y' or satisfied == 'Y':
        satisfied = True
        save_image(final_mosaic_crop, '6_final_mosaic_crop', save_all_steps)
    elif satisfied == 'No' or satisfied == 'no' or satisfied == 'n' or satisfied == 'N':
        satisfied = False
    else:
        raise Exception('Invalid entry, tou should enter: Yes or No')