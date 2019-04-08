import numpy as np
import math as m
import astropy
from astropy.io import fits
import os
import shutil

# Folder creation
path = os.path.dirname(os.path.abspath('__Test_shuttle_tri.py__'))
if not os.path.exists(path+'\Reduced'):  #
    os.mkdir(path+'\Reduced')

files=os.listdir(path+'\image')
print(files)
for file in files:
    if "fit" in file:
        Header = fits.getheader('file')
        print(type(Header))

        #Ecrire une ligne pour avoir la déclinaison et l ascenssion

        Darks=os.listdir(path+'\Dark_frame')
        Flats=os.listdir(path+'\Flat_field')