import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
import pysex
import glob
import os

FLX_bright = []
X_bright = []
Y_bright = []
X_bright_abs = [] #ensuring the correct star is picked
Y_bright_abs = [] #ensuring the correct star is picked


#find the position of each image from the header
def _find_pos(fitsFile):
   hdulist = fits.open(fitsFile)
   header = hdulist[0].header
   if "RA" in header and "DEC" in header:
      RA = header['RA']
      DEC = header['DEC']
      
   else:
      RA = 0   #some fits file headers don't have these values, 
      DEC = 0  #so if they don't set them as 0 for the time being untill i ask
      
   hdulist.close()
   return RA, DEC

####################################################
#Rename all objects that end with .fit into something Sextractor can open
#The new name is <number of fitsFile>.fit
i = 0
for filename in os.listdir("."):
   if filename.endswith(".fit"):
      i = i + 1
      new_name = str(i) + ".fit"
      os.rename(filename, new_name)
   else:
      continue
  

#Apply Sextractor to obtain the catalogues
fileList = glob.glob('*.fit') #or the directory to where the images are.
print(fileList)

i = 0 #counter
for fitsFile in fileList:

   cat = pysex.run(fitsFile, params=['X_IMAGE', 'Y_IMAGE', 'ALPHAPEAK_J2000', 'DELTAPEAK_J2000', 'FLUX_APER'], conf_args={'PHOT_APERTURES':5})


####################################################
catList = glob.glob('*.pysexcat')
print(catList)

## This is the main part of the code: a loop over all the catalogues extracted ##
i = 0 #define a counter
for catalogue in catList:
   
   data = np.loadtxt(catalogue, dtype = float, skiprows = 5)
   RA, DEC = _find_pos(fileList[i])   
   FLUX = data[:, 4]

   if RA == 0 and DEC == 0: 
      X = data[:, 0]
      Y = data[:, 1]
   else:
      X = data[:, 2]
      Y = data[:, 3]


   #remove any negative flux entries
   FLX_clean = FLUX[:].tolist()
   X_clean = X[:].tolist()
   Y_clean = Y[:].tolist()

   for f in FLUX:
      if f<0:
         ind = np.where(FLUX == f)

         FLX_clean.remove(f)
         X_clean.remove(X[ind])
         Y_clean.remove(Y[ind])         

   #find absolute position of the stars in the sky 
   X_abs = []
   Y_abs = []
   for n in range(len(FLX_clean)):
      X_abs.append(X_clean[n] + RA)
      Y_abs.append(Y_clean[n] + DEC)


   #find and append the brightest star in each image [absolute(rounded to 2 decimals) and regular coordinates]
   #First, round the numbers then append
   max_flx = max(FLX_clean)
   ind = FLX_clean.index(max_flx)
 
   X_bright_abs.append(round(X_abs[ind],2))
   Y_bright_abs.append(round(Y_abs[ind],2))
   
   FLX_bright.append(max_flx)
   X_bright.append(X_clean[ind])
   Y_bright.append(Y_clean[ind])

   i=i+1
      
###################################################################
## make sure the brightest star is the same one chosen in all the images ##

range_Xbright = 1.2*(np.median(X_bright_abs))
range_Ybright = 1.2*(np.median(Y_bright_abs))
for x in X_bright_abs:
   ind = X_bright_abs.index(x)
   y = Y_bright_abs[ind]
   #if the absolute coordinates of the brightest star in this catalogue is within 1.2xMedian of coordinates 
   #in all the catalogues, keep the catalogue, otherwise, it's not the same star ==> remove the values of that
   #star from the array. 
   if abs(x) <= range_Xbright and abs(y) <= range_Ybright:
      continue
   else:
      X_bright.remove(X_bright[ind])
      Y_bright.remove(Y_bright[ind])
      FLX_bright.remove(FLX_bright[ind])



print(FLX_bright)
print(X_bright)
print(Y_bright)
##############################################################
#plot the corresponding axes

fig = plt.figure()
ax = plt.axes(projection = '3d')

ax.scatter3D(X_bright, Y_bright, FLX_bright)
plt.show()

##############################################################
#fit to a polynomial







