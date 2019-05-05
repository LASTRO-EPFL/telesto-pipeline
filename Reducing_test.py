import numpy as np
import math as m
import astropy
from astropy.io import fits
import os
import shutil

# Folder creation
path = os.path.dirname(os.path.abspath('__Reducing_test.py__'))
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


files = os.listdir(path+'\images')
for file in files:

    if "fit" in file:
        hdul_data = fits.open(path+'\images\\'+file)  # open a FITS file
        hdr_data = hdul_data[0].header  # the primary HDU header

        file_data = fits.getdata(path+'\images\\'+file)
        numrow_data = len(file_data)
        numcol_data = len(file_data[0])
        Darks=os.listdir(path+'\Dark_frame')
        Flats=os.listdir(path+'\Flat_field')

        #Calculus of the corresponding master dark
        Mean_Dark_tot = 0
        Compteur_D = 0
        for Dark in Darks:
            Bool_D = 0
            Dark_data = fits.getdata(path+'\Dark_frame\\'+Dark)
            numrow_dark = len(Dark_data)
            numcol_dark = len(Dark_data[0])
            Mean_Dark=0
            if numrow_dark == numrow_data and numcol_dark == numrow_data:   #Check if the dark field and the datas has the same pixel size
                Compteur_D = Compteur_D+1
                Bool_D = 1
                Mean_Dark = np.sum(Dark_data)/(numrow_dark*numcol_dark)
                Mean_Dark_tot = Mean_Dark_tot+Mean_Dark

        if Bool_D == 1: #check that we don t divide by 0
            Mean_Dark_tot = Mean_Dark_tot/Compteur



        #Calculus of the MasterFlats
        Sum_F = np.zeros((numrow_data,numcol_data),dtype='i')
        for Flat in Flats:
            Bool_F = 0
            hdul_flat = fits.open(path+'\Flat_field\\'+Flat)  # open a FITS file
            hdr_flat = hdul_flat[0].header  # the primary HDU header
            Flat_data = fits.getdata(path+'\Flat_field\\'+Flat)
            numrow_flat = len(Flat_data)
            numcol_flat = len(Flat_data[0])
            Compteur_F = 0
            if hdr_data['FILTER'] == hdr_flat['FILTER']:
                if numrow_flat == numrow_data and numcol_flat == numrow_data:  # Check that the Flatfield and the data have the same pixel size
                    Bool_F = 1
                    Compteur_F = Compteur_F + 1
                    Flat_data = Flat_data-Mean_Dark_tot
                    Sum_F = Sum_F + Flat_data
        if Bool_F == 1:
            Mean_F = Sum_F/Compteur_F
            MasterFlat = Mean_F / np.median(Mean_F)


        if Bool_F == 1: #check if we actually used a flat field and thus not divide by 0
            file_data_reduced = (file_data - Mean_Dark_tot) / MasterFlat
        if Bool_F == 0:
            print('No Flatfield corresponding to '+file)

        if 'Halpha' in hdr_data['FILTER']:
            file.replace(".fit", "_reduced_Halpha.fits")
            hdu = fits.PrimaryHDU(file_data_reduced)
            hdul = fits.HDUlist([hdu])
            hdul.writeto(path+'\Reduced_Halpha\\'+file)

        if 'B' in hdr_data['FILTER']:
            file.replace(".fit", "_reduced_B.fits")
            hdu = fits.PrimaryHDU(file_data_reduced)
            hdul = fits.HDUlist([hdu])
            hdul.writeto(path + '\Reduced_B\\' + file)

        if 'R' in hdr_data['FILTER']:
            file.replace(".fit", "_reduced_R.fits")
            hdu = fits.PrimaryHDU(file_data_reduced)
            hdul = fits.HDUlist([hdu])
            hdul.writeto(path + '\Reduced_R\\' + file)

        if 'V' in hdr_data['FILTER']:
            file.replace(".fit", "_reduced_V.fits")
            hdu = fits.PrimaryHDU(file_data_reduced)
            hdul = fits.HDUlist([hdu])
            hdul.writeto(path + '\Reduced_V\\' + file)

        if 'sans filtre' in hdr_data['FILTER']:
            file.replace(".fit", "_reduced_without_filter.fits")
            hdu = fits.PrimaryHDU(file_data_reduced)
            hdul = fits.HDUlist([hdu])
            hdul.writeto(path + '\Reduced_without_filter\\' + file)

