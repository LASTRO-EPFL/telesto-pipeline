import os
import astroalign as aa
import astropy
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def align():
    path = os.path.dirname(os.path.abspath('__Main_Reducing_process.py__'))
    Folders = os.listdir(path)
    for Folder in Folders:
        if "Object" in Folder:
            images_to_align = os.listdir(path + '\\' + Folder)
            bool = 0
            compt_filter = 1
            list_filter = []

            for image in images_to_align:
                print(image)
                if "fit" in image and bool == 0:
                    hdul_image = fits.open(path + '\\' + Folder + '\\' + image)  # open a FITS file
                    hdr_image = hdul_image[0].header  # the primary HDU header
                    data_image = fits.getdata(path + '\\' + Folder + '\\' + image, ext=0)
                    bool=1
                    ref_image = data_image
                    print('reference image')
                    if not os.path.exists(path + '\\' + Folder + '\Aligned'):  #
                        os.mkdir(path + '\\' + Folder + '\Aligned')
                    fits.writeto(path + '\\' + Folder + '\Aligned' + '\\' + image, ref_image,
                                 header=hdr_image,
                                 overwrite=True)
                    list_filter.append(hdr_image['FILTER'])
                    hdul_image.close()

                if "fit" in image and bool == 1:
                    hdul_image = fits.open(path + '\\' + Folder + '\\' + image)  # open a FITS file
                    hdr_image = hdul_image[0].header  # the primary HDU header
                    data_image = fits.getdata(path + '\\' + Folder + '\\' + image, ext=0)
                    image_aligned = aa.register(data_image, ref_image)
                    fits.writeto(path + '\\' + Folder + '\Aligned' + '\\' + image, image_aligned,
                                 header=hdr_image,
                                 overwrite=True)
                    print('aligned')
                    bool_filter_check = 0
                    for i in range(0,compt_filter):
                        if hdr_image['FILTER'] == list_filter[i]:
                            bool_filter_check = 1

                    if bool_filter_check == 0:
                        compt_filter = compt_filter + 1
                        list_filter.append(hdr_image['FILTER'])

                    hdul_image.close()

            alligned_images = os.listdir(path+'\\'+Folder+'\Aligned')
            nb_filter = len(list_filter)
            if not os.path.exists(path + '\\' + Folder + '\Aligned' + '\Combined'):  #
                os.mkdir(path + '\\' + Folder + '\Aligned' + '\Combined')
            for i in range(0,nb_filter):
                test = 1
                Combined = []
                stack = []
                final = []
                for image in alligned_images:
                    if "fit" in image:
                        hdul_image = fits.open(path + '\\' + Folder + '\Aligned' + '\\' + image)  # open a FITS file
                        hdr_image = hdul_image[0].header  # the primary HDU header
                        data_image = fits.getdata(path + '\\' + Folder + '\Aligned' + '\\' + image, ext=0)
                        if hdr_image['FILTER'] == list_filter[i]:
                            Combined.append(data_image)
                            final_header = hdr_image
                        hdul_image.close()
                stack = np.dstack(Combined)
                final = np.mean(stack,2)
                fits.writeto(path + '\\' + Folder + '\Aligned\Combined' + '\\' + list_filter[i].split(' ')[0] + '.fits', final, header=final_header,
                                 overwrite=True)
