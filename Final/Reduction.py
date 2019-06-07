import numpy as np
import math as m
import astropy
from astropy.io import fits
import os
import shutil


def make_reduction():
    # Folder creation
    path = os.path.dirname(os.path.abspath('__Main_Reducing_process.py__'))
    if not os.path.exists(path+'\Object_1'):
        os.mkdir(path+'\Object_1')

    Darks = os.listdir(path + '\Dark_frame')
    Mean_Dark_tot = 0
    Compteur_D = 0
    Dark_array = []
    for Dark in Darks:
        hdul_dark = fits.open(path + '\Dark_frame\\' + Dark)  # open a FITS file
        hdr_dark = hdul_dark[0].header  # the primary HDU header
        Dark_data = hdul_dark[0].data
        Dark_data = Dark_data / hdr_dark['EXPTIME']
        numrow_dark = len(Dark_data)
        numcol_dark = len(Dark_data[0])
        Dark_data = Dark_data[30:numrow_dark - 30, 30:numcol_dark - 30] #Solving edge problems (Black lines)

        Compteur_D = Compteur_D + 1
        Median_Dark = np.median(Dark_data)
        Dark_array.append(Median_Dark)

        hdul_dark.close()

    if Compteur_D != 0: #check that we don't divide by 0
        Mean_Dark_tot = np.median(Dark_array)

    if Compteur_D == 0:
        print('There is no dark frame available to reduce the image.')

    files = os.listdir(path+'\images')
    Compteur_object = 1
    First_object = 1
    Mean_Dark_tot_before = Mean_Dark_tot
    for file in files:
        if "fit" in file:
            print(file)
            hdul_data = fits.open(path+'\images\\'+file)  # open a FITS file
            hdr_data = hdul_data[0].header  # the primary HDU header
            file_data_temp = fits.getdata(path+'\images\\'+file)

            numrow_data = len(file_data_temp)
            numcol_data = len(file_data_temp[0])
            file_data = file_data_temp[30:numrow_data-30,30:numcol_data-30]

            Mean_Dark_tot = Mean_Dark_tot_before * hdr_data['EXPTIME']

            Flats = os.listdir(path + '\Flat_field')
            Bool_F = 0
            Stack = []
            numpy_stack = []
            for Flat in Flats:
                hdul_flat = fits.open(path + '\Flat_field\\' + Flat)  # open a FITS file
                hdr_flat = hdul_flat[0].header  # the primary HDU header
                Flat_data = hdul_flat[0].data
                numrow_flat = len(Flat_data)
                numcol_flat = len(Flat_data[0])
                Flat_data = Flat_data[30:numrow_flat - 30, 30:numcol_flat - 30]
                if hdr_data['FILTER'] == hdr_flat['FILTER']:
                    if numrow_flat == numrow_data and numcol_flat == numcol_data:  # Check that the Flatfield and the data have the same pixel size
                        Bool_F = 1
                        Flat_data_median = (Flat_data - Mean_Dark_tot) / np.median(Flat_data)
                        Stack.append(Flat_data_median)
                hdul_flat.close()

            if Bool_F == 1:  # check if we actually used a flat field and thus not divide by 0
                numpy_stack = np.dstack(Stack)
                MasterFlat = np.median(numpy_stack, 2)

                file_data_reduced = (file_data - Mean_Dark_tot) / MasterFlat
                print('Reduced successfully')
                hdu = fits.PrimaryHDU(file_data_reduced)
                hdul = fits.HDUList([hdu])

                if First_object == 0:
                    Bool_check = 0
                    for f in range(1, Compteur_object + 1):
                        Checks = os.listdir(path + '\Object_' + str(Compteur_object))
                        for check in Checks:
                            if "fit" in check:
                                hdul_check = fits.open(path + '\Object_' + str(Compteur_object) + '\\' + check)  # open a FITS file
                                hdr_check = hdul_check[0].header  # the primary HDU header
                                if hdr_data['OBJCTRA'][:7] == hdr_check['OBJCTRA'][:7] and hdr_data['OBJCTDEC'][:7] == hdr_check['OBJCTDEC'][:7] and Bool_check == 0: #Precision on the same position on the sky
                                            fits.writeto(path + '\Object_' + str(Compteur_object) + '\\' + file, file_data_reduced, header=hdr_data,
                                                 overwrite=True)
                                            Bool_check = 1
                                            hdul_check.close()

                    if Bool_check == 0:
                        Compteur_object = Compteur_object + 1
                        os.mkdir(path + '\Object_' + str(Compteur_object))
                        fits.writeto(path + '\Object_' + str(Compteur_object) + '\\' + file + 's', file_data_reduced, header=hdr_data,
                                         overwrite=True)
                if First_object == 1:
                    fits.writeto(path + '\Object_1\\' + file + 's', file_data_reduced, header=hdr_data, overwrite=True)
                    First_object = 0

            if Bool_F == 0:
                print('No Flatfield corresponding to ' + file)

            hdul_data.close()
