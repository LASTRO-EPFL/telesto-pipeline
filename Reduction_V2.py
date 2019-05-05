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
if not os.path.exists(path+'\Reduced'):  #
    os.mkdir(path+'\Reduced')
if not os.path.exists(path+'\Reduced_V'):  #
    os.mkdir(path+'\Reduced_V')
if not os.path.exists(path+'\Reduced_B'):  #
    os.mkdir(path+'\Reduced_B')
if not os.path.exists(path+'\Dark'):  #
    os.mkdir(path+'\Dark')
if not os.path.exists(path+'\Flat'):  #
    os.mkdir(path+'\Flat')
if not os.path.exists(patj+'\Object_1')
    os.mkdir(path+'\Object_1')

Darks = os.listdir(path + '\Dark_frame')
Mean_Dark_tot = 0
Compteur_D = 0
for Dark in Darks:
    #Dark_data_temp = fits.getdata(path+'\Dark_frame\\'+Dark)
    hdul_dark = fits.open(path + '\Dark_frame\\' + Dark)  # open a FITS file
    hdr_dark = hdul_dark[0].header  # the primary HDU header
    Dark_data = hdul_dark[0].data
    Dark_data = Dark_data / hdr_dark['EXPTIME'] *20 #PROBLEME A REGLER!!!!!!!!!!!!!!!!!!!
    numrow_dark = len(Dark_data)
    numcol_dark = len(Dark_data[0])
    Dark_data = Dark_data[30:numrow_dark - 30, 30:numcol_dark - 30]

    Compteur_D = Compteur_D + 1
    Mean_Dark = Dark_data.mean()
    Mean_Dark_tot = Mean_Dark_tot+Mean_Dark

    hdul_dark.close()

if Compteur_D != 0: #check that we don t divide by 0
    Mean_Dark_tot = Mean_Dark_tot/Compteur_D
    hdu_D = fits.PrimaryHDU(Mean_Dark_tot)
    hdul_D = fits.HDUList([hdu_D])
    hdul_D.writeto(path + '\Dark\\' + 'Dark_frame_tot.fit', overwrite=True)
    print(Mean_Dark_tot)

if Compteur_D == 0:
    print('There is no dark frame available to reduce the image.')


#Calculus of the Masterflat
Flats = os.listdir(path + '\Flat_field')
Bool_Halpha = 0
Bool_B = 0
Bool_V = 0
Bool_R = 0
Bool_without_filer = 0
Stack_Haplha = []
Stack_B = []
Stack_V = []
Stack_R = []
Stack_without_filter = []



for Flat in Flats:
    hdul_flat = fits.open(path+'\Flat_field\\'+Flat)  # open a FITS file
    hdr_flat = hdul_flat[0].header  # the primary HDU header
    Flat_data = hdul_flat[0].data

    numrow_flat = len(Flat_data)
    numcol_flat = len(Flat_data[0])
    Flat_data = Flat_data[30:numrow_flat - 30, 30:numcol_flat - 30]


    if 'Halpha' in hdr_flat['FILTER']:
        print(hdr_flat['FILTER'])
        Bool_Halpha = 1
        Flat_data_median = (Flat_data - Mean_Dark_tot)/np.median(Flat_data)
        Stack_Haplha.append(Flat_data_median)

    if 'B' in hdr_flat['FILTER']:
        Bool_B = 1
        Flat_data_median = (Flat_data - Mean_Dark_tot)/np.median(Flat_data)
        Stack_B.append(Flat_data_median)

    if 'V' in hdr_flat['FILTER']:
        Bool_V = 1
        Flat_data_median = (Flat_data - Mean_Dark_tot)/np.median(Flat_data)
        Stack_V.append(Flat_data_median)

    if 'R' in hdr_flat['FILTER']:
        Bool_R = 1
        Flat_data_median = (Flat_data - Mean_Dark_tot)/np.median(Flat_data)
        Stack_R.append(Flat_data_median)

    if 'sans filtre' in hdr_flat['FILTER']:
        Bool_without_filer = 1
        Flat_data_median = (Flat_data - Mean_Dark_tot) / np.median(Flat_data)
        Stack_without_filter.append(Flat_data_median)


    hdul_flat.close()

if Bool_Halpha == 1:  #check if we actually used a flat field
    numpy_stack = np.dstack(Stack_Haplha)
    MasterFlat_Halpha = np.median(numpy_stack, 2)
    hdu_F = fits.PrimaryHDU(MasterFlat_Halpha)
    hdul_F = fits.HDUList([hdu_F])
    hdul_F.writeto(path + '\Flat\\' + 'MasterFlat_Halpha.fit',overwrite=True)

if Bool_B == 1:  #check if we actually used a flat field
    numpy_stack = np.dstack(Stack_B)
    MasterFlat_B = np.median(numpy_stack, 2)
    hdu_F = fits.PrimaryHDU(MasterFlat_B)
    hdul_F = fits.HDUList([hdu_F])
    hdul_F.writeto(path + '\Flat\\' + 'MasterFlat_B.fit',overwrite=True)

if Bool_V == 1:  #check if we actually used a flat field
    numpy_stack = np.dstack(Stack_V)
    MasterFlat_V = np.median(numpy_stack, 2)
    hdu_F = fits.PrimaryHDU(MasterFlat_V)
    hdul_F = fits.HDUList([hdu_F])
    hdul_F.writeto(path + '\Flat\\' + 'MasterFlat_V.fit',overwrite=True)

if Bool_R == 1:  #check if we actually used a flat field
    numpy_stack = np.dstack(Stack_R)
    MasterFlat_R = np.median(numpy_stack, 2)
    hdu_F = fits.PrimaryHDU(MasterFlat_R)
    hdul_F = fits.HDUList([hdu_F])
    hdul_F.writeto(path + '\Flat\\' + 'MasterFlat_R.fit',overwrite=True)

if Bool_without_filer == 1:  #check if we actually used a flat field
    numpy_stack = np.dstack(Stack_without_filter)
    MasterFlat_without_filter = np.median(numpy_stack, 2)
    hdu_F = fits.PrimaryHDU(MasterFlat_without_filter)
    hdul_F = fits.HDUList([hdu_F])
    hdul_F.writeto(path + '\Flat\\' + 'MasterFlat_without_filter.fit',overwrite=True)


files = os.listdir(path+'\images')
Compteur_object = 1
First_object = 1
for file in files:
    if "fit" in file:
        print(file)
        hdul_data = fits.open(path+'\images\\'+file)  # open a FITS file
        hdr_data = hdul_data[0].header  # the primary HDU header
        file_data_temp = fits.getdata(path+'\images\\'+file)
        numrow_data = len(file_data_temp)
        numcol_data = len(file_data_temp[0])
        file_data = file_data_temp[30:numrow_data-30,30:numcol_data-30]

        Bool_reduced = 0
        if numcol_data == 2086 and numrow_data == 2074: # A REGARDE COMMENT TRAITER LA TAILLE DES IMAGES
            if 'Halpha' in hdr_data['FILTER'] and Bool_Halpha == 1:
                file_data_reduced = (file_data - Mean_Dark_tot) / MasterFlat_Halpha
                Bool_reduced = 1
            if 'B' in hdr_data['FILTER'] and Bool_B == 1:
                file_data_reduced = (file_data - Mean_Dark_tot) / MasterFlat_B
                Bool_reduced = 1
            if 'V' in hdr_data['FILTER'] and Bool_V == 1:
                file_data_reduced = (file_data - Mean_Dark_tot) / MasterFlat_V
                Bool_reduced = 1
            if 'R' in hdr_data['FILTER'] and Bool_R == 1:
                file_data_reduced = (file_data - Mean_Dark_tot) / MasterFlat_R
                Bool_reduced = 1
            if 'sans filtre' in hdr_data['FILTER'] and Bool_without_filer == 1:
                file_data_reduced = (file_data - Mean_Dark_tot) / MasterFlat_without_filter
                Bool_reduced = 1

        if Bool_reduced == 1:
            hdu = fits.PrimaryHDU(file_data_reduced)
            hdul = fits.HDUList([hdu])
            hdul.writeto(path + '\Reduced\\' + file, overwrite=True)

        else:
            print('Not able to reduce the following file : ' + file)

        hdul_data.close()







'''

        if First_object == 1:
            hdu = fits.PrimaryHDU(file_data_reduced)
            hdul = fits.HDUList([hdu])
            hdul.writeto(path + '\Object_1\\' + file)
            First_object = 2

        if First_object == 2:
            for f in range(1, Compteur_object):
                Checks = os.listdir(path + '\Object_' + Compteur_object)
                for check in Checks:
                    if
'''
#chaque flat diviser par sa mediane. np.stack, numpy je prends la median selon le 3 eme axe
#Pour les darks, faire comme les flat sans diviser par la median initiale
#mettre les fit.close pour eviter les memory errors
#overwrite




#Probleme avec le temps d exposition pour le soustraire au flat...
#probleme avec le triage d objet puisque quand j ai reduis j ai plus les infos...
#du coup probleme avec les filtres
