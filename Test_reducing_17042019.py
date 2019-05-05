import numpy as np
import math as m
import astropy
from astropy.io import fits
import os
import shutil

# Folder creation
path = os.path.dirname(os.path.abspath('__Test_reducing_17042019.py__'))
if not os.path.exists(path+'\Reduced_without_filter'):  #
    os.mkdir(path+'\Reduced_without_filter')
if not os.path.exists(path+'\Reduced_Halpha'):  #
    os.mkdir(path+'\Reduced_Halpha')
if not os.path.exists(path+'\Reduced_R'):  #
    os.mkdir(path+'\Reduced_R')
if not os.path.exists(path+'\Reduced_V'):  #
    os.mkdir(path+'\Reduced_V')
if not os.path.exists(path+'\Reduced_B'):  #
    os.mkdir(path+'\Reduced_B')
if not os.path.exists(path+'\Dark'):  #
    os.mkdir(path+'\Dark')
if not os.path.exists(path+'\Flat'):  #
    os.mkdir(path+'\Flat')


files = os.listdir(path+'\images')

for file in files:

    if "fit" in file:
        print(file)
        hdul_data = fits.open(path+'\images\\'+file)  # open a FITS file
        hdr_data = hdul_data[0].header  # the primary HDU header

        file_data_temp = fits.getdata(path+'\images\\'+file)
        numrow_data = len(file_data_temp)
        numcol_data = len(file_data_temp[0])
        file_data = file_data_temp[30:numrow_data-30,30:numcol_data-30]
        Darks=os.listdir(path+'\Dark_frame')
        Flats=os.listdir(path+'\Flat_field')

        #Calculus of the corresponding master dark
        Mean_Dark_tot = 0
        Compteur_D = 0
        Stack_D = []
        Bool_D = 0
        for Dark in Darks:
            #Dark_data_temp = fits.getdata(path+'\Dark_frame\\'+Dark)
            hdul_dark = fits.open(path + '\Dark_frame\\' + Dark)  # open a FITS file
            hdr_dark = hdul_dark[0].header  # the primary HDU header
            Dark_data = hdul_dark[0].data
            Dark_data = Dark_data * hdr_data['EXPTIME'] / hdr_dark['EXPTIME']

            numrow_dark = len(Dark_data)
            numcol_dark = len(Dark_data[0])

            Dark_data = Dark_data[30:numrow_dark - 30, 30:numcol_dark - 30]

            if numrow_dark == numrow_data and numcol_dark == numcol_data:   #Check if the dark field and the datas has the same pixel size
                Bool_D = 1
                Compteur_D = Compteur_D + 1
                #Stack_D.append(Dark_data)
                #print(Dark_data)
                #print(len(Dark_data[0]))
                #print(np.size(Dark_data))
                #print(np.shape(Dark_data)[0]*np.shape(Dark_data)[1])
                #print(np.sum(Dark_data))
                Mean_Dark = Dark_data.mean() # np.sum(Dark_data)/(len(Dark_data)*len(Dark_data[0]))
                #print(Mean_Dark)
                Mean_Dark_tot = Mean_Dark_tot+Mean_Dark

                hdul_dark.close()

        if Bool_D == 1: #check that we don t divide by 0
            Mean_Dark_tot = Mean_Dark_tot/Compteur_D
            #numpy_stack_D = np.dstack(Stack_D)
            #Mean_Dark_tot = np.median(numpy_stack_D, 2)
            print(Mean_Dark_tot)

            hdu_D = fits.PrimaryHDU(Mean_Dark_tot)
            hdul_D = fits.HDUList([hdu_D])
            hdul_D.writeto(path + '\Dark\\' + file,overwrite=True)

        #Calculus of the MasterFlats
        Sum_F = np.zeros((numrow_data,numcol_data),dtype='i')
        Bool_F = 0
        Stack = []
        for Flat in Flats:
            hdul_flat = fits.open(path+'\Flat_field\\'+Flat)  # open a FITS file
            hdr_flat = hdul_flat[0].header  # the primary HDU header
            #Flat_data_temp = fits.getdata(path+'\Flat_field\\'+Flat)
            Flat_data = hdul_flat[0].data
            #Flat_data = np.array(Flat_data_temp)
            numrow_flat = len(Flat_data)
            numcol_flat = len(Flat_data[0])
            Flat_data = Flat_data[30:numrow_flat - 30, 30:numcol_flat - 30]
            if hdr_data['FILTER'] == hdr_flat['FILTER']:
                print(hdr_flat['FILTER'])
                if numrow_flat == numrow_data and numcol_flat == numcol_data:  # Check that the Flatfield and the data have the same pixel size
                    Bool_F = 1
                    Flat_data_median = (Flat_data - Mean_Dark_tot)/np.median(Flat_data)
                    Stack.append(Flat_data_median)
            hdul_flat.close()

        if Bool_F == 1:  #check if we actually used a flat field and thus not divide by 0
            numpy_stack = np.dstack(Stack)
            MasterFlat = np.median(numpy_stack, 2)
            hdu_F = fits.PrimaryHDU(MasterFlat)
            hdul_F = fits.HDUList([hdu_F])
            hdul_F.writeto(path + '\Flat\\' + file,overwrite=True)

            file_data_reduced = (file_data - Mean_Dark_tot) / MasterFlat
            print('Reduced successfully')
            if 'Halpha' in hdr_data['FILTER']:
                string = file.replace(".fit", "_reduced_Halpha.fits")
                hdu = fits.PrimaryHDU(file_data_reduced)
                hdul = fits.HDUList([hdu])
                hdul.writeto(path + '\Reduced_Halpha\\' + string ,overwrite=True)

            if 'B' in hdr_data['FILTER']:
                file.replace(".fit", "_reduced_B.fits")
                hdu = fits.PrimaryHDU(file_data_reduced)
                hdul = fits.HDUList([hdu])
                hdul.writeto(path + '\Reduced_B\\' + file)

            if 'R' in hdr_data['FILTER']:
                file.replace(".fit", "_reduced_R.fits")
                hdu = fits.PrimaryHDU(file_data_reduced)
                hdul = fits.HDUList([hdu])
                hdul.writeto(path + '\Reduced_R\\' + file)

            if 'V' in hdr_data['FILTER']:
                file.replace(".fit", "_reduced_V.fits")
                hdu = fits.PrimaryHDU(file_data_reduced)
                hdul = fits.HDUList([hdu])
                hdul.writeto(path + '\Reduced_V\\' + file)

            if 'sans filtre' in hdr_data['FILTER']:
                file.replace(".fit", "_reduced_without_filter.fits")
                hdu = fits.PrimaryHDU(file_data_reduced)
                hdul = fits.HDUList([hdu])
                hdul.writeto(path + '\Reduced_without_filter\\' + file)

        if Bool_F == 0:
            print('No Flatfield corresponding to '+file)

        hdul_data.close()

#chaque flat diviser par sa mediane. np.stack, numpy je prends la median selon le 3 eme axe
#Pour les darks, faire comme les flat sans diviser par la médian initiale
#mettre les fit.close pour éviter les memory errors
#overwrite
