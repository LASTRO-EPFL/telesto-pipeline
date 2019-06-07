from Sorting_images import *
from Reduction import *
from Alignment import *
from Gaussian_filter import *


make_sort()             #Sort all the images in 3 different folders : Images, Dark_frame and Flat_field.
make_reduction()        #Reduce all the raw images ("Subtraction" of Master Dark and Master Flat). Sort the output files by object.
align()                 #Align all the images and the combine them in one single image per filter.
make_gaussian_filter()  #Solve the different resolution of the filters.
