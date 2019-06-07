# Mosaic of images

Creation of a mosaic of images for Telesto, the 60 cm telescope of the Geneva Observatory &amp; EPFL. You should have four images of the same astronomical object, one approximately covering the top_left part of the object, and the others covering the top_right, bottom_left and bottom_right parts of the object.  


### Language
Python 3, package `pytelesto`.

### How it works
1) Create a file named 'Mosaic' where to put the four images for the future mosaic. Name them 'image_top_right', 'image_top_left', 'image_bottom_left' and 'image_bottom_right' according to their future position in the Mosaic. 

2) Don't forget to install the astroalign package from github (https://github.com/toros-astro/astroalign) before running the code.


NB: -You can choose whether or not you want to save all the steps of the mosaic creation. If Yes  it will save in particular the alignement of the image top right and top left, then the alignment with the image bottom_left etc.

     -At the end, you can crop the final image at the desired size. 

     -You will also be asked to set the boundary limits of the intensity's interval of all your pixels. It allows to redefined contrasts. The pixels with values outside the interval are reset to the closer boundary value.
 
     -If you are not satisfied, several parameters of the Mosaic_of_images.py can changed: map, cut_boarder, add_frame, width, height. 

# Image_combination

Combined four sets of patches into four images. The combination is optimized, minimizing twice the rms value on the sky background. The patches are aligned, and the images resulting are strored in a file called "Mosaic". A second set of four images is also created, without optimization, thus using all patches. The output images are stored in a file "Mosaic_all_patch".

### Language
Python 3, package `pytelesto`.

### How it works
1) In Image_combination_functions.py, in function open_all_images, change the path name so the program finds your patches. 

2) That's all!







