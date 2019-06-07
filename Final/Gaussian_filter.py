import numpy as np
import math as m
import astropy
from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize as optimize
import os
import shutil


def __twoDGauss(tabxy,amplitude,x0,y0,sigma,offset):
    x=tabxy[0]
    y=tabxy[1]
    return offset + amplitude*np.exp(-(((x-x0)**2+(y-y0)**2)/(2*sigma**2)))

def make_gaussian_filter():
    path = os.path.dirname(os.path.abspath('__Main_Reducing_process.py__'))
    Folders = os.listdir(path)
    for Folder in Folders:
        if "Object" in Folder:
            images_to_filter = os.listdir(path + '\\' + Folder + '\Aligned' + '\Combined')
            sigma_filter = []
            filter = []
            for image in images_to_filter:
                Nb_stars = 0
                sigma = []
                if "fit" in image:
                    while Nb_stars != 3: #Number of fit Gaussian per filter
                        print('Please introduce the position ((x,y) on DS9)) of the center of three distinguishable stars (precise center)'
                            'On the following image : ' + path + '\\' + Folder + '\Aligned' + '\Combined' + "\\" + image)
                        y1 = int(round(float(input('x1='))))
                        x1 = int(round(float(input('y1='))))
                        data_image = fits.getdata(path + '\\' + Folder + '\Aligned\Combined' + '\\' + image, ext=0)
                        xpos1 = []
                        ypos1 = []
                        zpos1 = []
                        for i in range(0, 13):
                            for j in range(0, 13):
                                xpos1.append(i + x1 - 7)
                                ypos1.append(j + y1 - 7)
                                zpos1.append(data_image[i + x1 - 7, j + y1 - 7])

                        params1, pcov1 = optimize.curve_fit(__twoDGauss, [xpos1, ypos1], zpos1, p0=[20000, x1, y1, 1, 0])
                        fig1 = plt.figure()
                        ax1 = fig1.add_subplot(111, projection='3d')
                        # ax1.plot_trisurf(haxpos1,haypos1,hazpos1) #For a 3D plot
                        ax1.plot(xpos1, ypos1, zpos1)
                        ax1.plot(xpos1, ypos1,
                                 __twoDGauss([xpos1, ypos1], params1[0], params1[1], params1[2], params1[3], params1[4]))
                        print('Does the fit seems ok? (CLOSE THE IMAGE AND THEN PRESS Y OR N)')
                        plt.show()
                        test = input()
                        if test == 'y' or test == 'Y':
                            Nb_stars = Nb_stars + 1
                            sigma.append(params1[3])
                        if test != 'y' and test != 'Y' and test !='n' and test != 'N':
                            print('Error, the input needs a "Y" or a "N"')

                    sigma_filter.append(np.mean(sigma))

            if not os.path.exists(path + '\\' + Folder + '\Aligned' + '\Combined' + '\Gaussian_filter'):  #
                os.mkdir(path + '\\' + Folder + '\Aligned' + '\Combined' + '\Gaussian_filter')

            print(sigma_filter)
            max_pos = sigma_filter.index(max(sigma_filter))
            print(max_pos)
            c = 0
            for image in images_to_filter:
                if "fit" in image:
                    data_image_gaussian = fits.getdata(path + '\\' + Folder + '\Aligned' + '\\' + '\Combined' + '\\'
                                                       + image, ext=0)
                    hdul_image = fits.open(path + '\\' + Folder + '\Aligned' + '\\' + '\Combined' + '\\'
                                                       + image)  # open a FITS file
                    hdr_image = hdul_image[0].header  # the primary HDU header
                    hdr_image['RA'] = hdr_image['OBJCTRA']
                    hdr_image['DEC'] = hdr_image['OBJCTDEC']
                    if c != max_pos:
                        sigma = abs(max(sigma_filter)-sigma_filter[c])
                        data_gaussian_filter = gaussian_filter(data_image_gaussian, sigma=sigma)
                        fits.writeto(path + '\\' + Folder + '\Aligned\Combined\Gaussian_filter' + '\\' + image,
                                     data_gaussian_filter, header=hdr_image, overwrite=True)
                    if c == max_pos:
                        fits.writeto(path + '\\' + Folder + '\Aligned\Combined\Gaussian_filter' + '\\' + image,
                                    data_image_gaussian, header=hdr_image, overwrite=True)
                    hdul_image.close()

                    c = c + 1
