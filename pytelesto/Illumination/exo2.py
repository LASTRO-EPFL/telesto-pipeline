import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
import pysex
import glob
import os

#pixel size
focal_length = 2280. #in mm
CCD_pixel_size = 9. #in microns
px = (CCD_pixel_size/focal_length)*206.3

def _find_pos(fitsFile):
   hdulist = fits.open(fitsFile)
   header = hdulist[0].header
   if "OBJCTRA" in header and "OBJCTDEC" in header:
      R = header['OBJCTRA'].split()
      D = header['OBJCTDEC'].split()

      h,m,s = float(R[0]), float(R[1]), float(R[2])
      deg,arc_m,arc_s = float(D[0]), float(D[1]), float(D[2])

      RA = (15.*(h*3600.))+ (15.*(m*60.)) + (15.*s)
      DEC = (deg*3600) + (arc_m*60.) + (arc_s) 
      RA = RA/px
      DEC = DEC/px
   else:
      raise Exception("Right Ascension and Declination of Image could not be found")
      
   hdulist.close()
   return RA, DEC

#use pysex to extract all stars in every catalogue

fileList = sorted(glob.glob('*.fit')) #or the directory to where the images (to be extracted) are.
print(fileList)
for fitsFile in fileList:
   cat = pysex.run(fitsFile, params=['X_IMAGE', 'Y_IMAGE', 'FLUX_APER', 'FLUXERR_APER'], conf_args={'PHOT_APERTURES':5})

catList = sorted(glob.glob('*.pysexcat'))

star_flux = []
star_error = []
x = []
y = []

xs = 1370.
ys = 1000.

#iRA0 = 832331.6461463886
#iDEC0 = 198302.75133301015

for catalogue in catList:   
 
   i = 0 #counter
   data = np.loadtxt(catalogue, dtype = float, skiprows = 4)
   iRA, iDEC = _find_pos(fileList[i])   #image RA, DEC in pixel
   print("image coord", iRA, iDEC)
   X = np.array(data[:, 0], dtype = 'float')
   Y = np.array(data[:, 1], dtype = 'float')
   FLUX = np.array(data[:, 2], dtype = 'float')
   FLUXERR = np.array(data[:, 3], dtype = 'float')
   m = len(X)

   XX = X #+ iRA - iRA0
   YY = Y #+ iDEC - iDEC0

   dist = []
   for k in range(m):
      dist.append(np.sqrt((xs - XX[k])**2 + (ys - YY[k])**2))
   
   ind = dist.index(min(dist))   
   star_flux.append(FLUX[ind])
   print(X[ind], Y[ind])
   i+=1

plt.plot(star_flux, '.')
plt.xlabel("Image index")
plt.ylabel("Flux (count)")
plt.show()








