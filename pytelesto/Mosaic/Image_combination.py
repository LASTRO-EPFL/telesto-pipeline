"""
IMAGE_COMBINATION: Combination and alignment of all the patches of a same portion
of the sky, into one single images (mean of the images pixel by pixel).
"""

from Image_combination_functions import image_combination_choosing_patch, image_combination_using_all_patch
from Mosaic_basic_function import save_image
import os

os.makedirs('Mosaic', exist_ok=True)
os.makedirs('Mosaic_all_patch', exist_ok=True)

save_all_steps = True

# Optimizing the noise:
image1 = image_combination_choosing_patch(1, save_all_steps)
image2 = image_combination_choosing_patch(2, save_all_steps)
image3 = image_combination_choosing_patch(3, save_all_steps)
image4 = image_combination_choosing_patch(4, save_all_steps)

save_image(image1, '.\Mosaic\image_top_right', True)
save_image(image2, '.\Mosaic\image_top_left', True)
save_image(image3, '.\Mosaic\image_bottom_left', True)
save_image(image4, '.\Mosaic\image_bottom_right', True)

# Using all patches:
image1_all_patch = image_combination_using_all_patch(1, save_all_steps)
image2_all_patch = image_combination_using_all_patch(2, save_all_steps)
image3_all_patch = image_combination_using_all_patch(3, save_all_steps)
image4_all_patch = image_combination_using_all_patch(4, save_all_steps)

save_image(image1, '.\Mosaic_all_patch\image_top_right', True)
save_image(image2, '.\Mosaic_all_patch\image_top_left', True)
save_image(image3, '.\Mosaic_all_patch\image_bottom_left', True)
save_image(image4, '.\Mosaic_all_patch\image_bottom_right', True)