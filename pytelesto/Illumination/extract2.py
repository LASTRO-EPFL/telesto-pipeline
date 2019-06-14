import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
import pysex
import glob
import os
import scipy.linalg
import astropy.units as u
from astropy.coordinates import SkyCoord

#calculate the pixel size of the telescope:
#Telesto example:

focal_length = 2280. #in mm
CCD_pixel_size = 9. #in microns
px = (CCD_pixel_size/focal_length)*206.3   #usable pixel size 
print("px",px)

FLX_bright = []
FLXERR_bright = []
X_bright = []
Y_bright = []
X_bright_abs = [] #ensuring the correct star is picked
Y_bright_abs = [] #ensuring the correct star is picked


#find the position of each image from the header
def _find_pos(fitsFile):
   hdulist = fits.open(fitsFile)
   header = hdulist[0].header
   if "OBJCTRA" in header and "OBJCTDEC" in header:
      RA = float(header['OBJCTRA'].replace(" ", ""))
      DEC = float(header['OBJCTDEC'].replace(" ", ""))
 
   else:
      raise Exception("Right Ascension and Declination of Image could not be found")
      #some fits file headers might not have these values, 
      #discard these images if so (do this later)
      
   hdulist.close()
   return RA, DEC

#order = 4   # 1,2,3,4
def fit_func(order):
   if order == 1:
      # best-fit linear plane
      A = np.c_[np.ones(data.shape[0]), data[:,0], data[:,1]]
      C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
    
      return C 

   elif order == 2:
      # best-fit quadratic curve
      A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
      C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
    
      return C

   elif order == 3:

      a = data[:,0]
      b = data[:,1]
      a2 = a**2
      b2 = b**2
   
      A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2, a*b2, a2*b , data[:,:2]**3] 
      C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])   #coefficients
   
      return C

   elif order == 4:

      a = data[:,0]
      b = data[:,1]
      a2 = a**2
      b2 = b**2  
      a3 = a**3
      b3 = b**3
   
      A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2, a*b2, a2*b , data[:,:2]**3, a*b3, a3*b, a2*b2, data[:,:2]**4] 
      C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])   #coefficients
   
      return C


####################################################

#Apply Sextractor to obtain the catalogues
fileList = glob.glob('*.fit') #or the directory to where the images (to be extracted) are.
print(fileList)

for fitsFile in fileList:
   cat = pysex.run(fitsFile, params=['X_IMAGE', 'Y_IMAGE', 'FLUX_APER', 'FLUXERR_APER'], conf_args={'PHOT_APERTURES':5})

#X_IMAGE
#Y_IMAGE
#ALPHAPEAK_J2000
#DELTAPEAK_J2000

####################################################
catList = glob.glob('*.pysexcat')
print(catList)


## This is the main part of the code: a loop over all the catalogues extracted ##
i = 0 #define a counter
for catalogue in catList:
   
   data = np.loadtxt(catalogue, dtype = float, skiprows = 4)
   RA, DEC = _find_pos(fileList[i])   

   X = data[:, 0]
   Y = data[:, 1]
   FLUX = data[:, 2]
   FLUXERR = data[:, 3]

   #remove any negative flux entries
   FLX_clean = FLUX[:].tolist()
   FLXERR_clean = FLUXERR[:].tolist()
   X_clean = X[:].tolist()
   Y_clean = Y[:].tolist()
   

   for f in FLUX:
      if f<0:
         ind = np.where(FLUX == f)

         FLX_clean.remove(f)
         FLXERR_clean.remove(FLX_clean[ind])
         X_clean.remove(X[ind])
         Y_clean.remove(Y[ind])         



   max_flx = max(FLX_clean)
   ind = FLX_clean.index(max_flx)
    
   FLX_bright.append(max_flx)
   FLXERR_bright.append(FLXERR_clean[ind])
   X_bright.append(X_clean[ind])   #Use these to plot in the end
   Y_bright.append(Y_clean[ind])   #

 
   med = np.mean(FLX_bright)
   std = np.std(FLX_bright)
   #print(med, std) 
   f_OK = [] #flux  
   e_OK = [] #error
   x_OK = [] #x
   y_OK = [] #y

   for f in FLX_bright:
      if f > (med - 2*std):
         ind = FLX_bright.index(f)
         f_OK.append(f)  
         e_OK.append(FLXERR_bright[ind])
         x_OK.append(X_bright[ind])
         y_OK.append(Y_bright[ind])    
      else:
         continue

   i=i+1
 
##############################################################
#plot the corresponding axes
'''
fig = plt.figure()
ax = plt.axes(projection = '3d')

ax.scatter3D(x_OK, y_OK, f_OK)
plt.show()
'''
##############################################################
#function to fit 3D data to a polynomial

#data
n = len(x_OK)
data = np.column_stack((x_OK, y_OK, f_OK))
print("data", data)

mn = np.min(data, axis = 0)
mx = np.max(data, axis = 0)
x,y = np.meshgrid(np.linspace(mn[0], mx[0], n), np.linspace(mn[1], mx[1], n))

XX = x.flatten()
YY = y.flatten()


# plot points and fitted surface
order = 1
C = fit_func(order)
print("C", C)

fig = plt.figure()
#plt.contourf(x, y, Z, origin = 'lower')
#plt.imshow(Z, origin = 'lower')
ax = fig.gca(projection='3d')
ax.plot_surface(x_OK, y_OK, Z, alpha = 0.2)
ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
plt.xlabel('X')
plt.ylabel('Y')
ax.set_zlabel('Z')

plt.show()

###########################################################
#check the quality of the fit using X^2 method

x_OK
y_OK
n = len(x_OK)
m = len(C)
f_OK
e_OK

z_OK = []

if order == 1:
   for i in range(n):
      z_OK.append(C[0] + C[1]*x_OK[i] + C[2]*y_OK[i])

elif order == 2:
   for i in range(n):
      z_OK.append(C[0] + C[1]*x_OK[i] + C[2]*y_OK[i] + C[3]*x_OK[i]*y_OK[i] + C[4]*(x_OK[i]**2) + C[5]*(y_OK[i]**2))

elif order == 3:
   for i in range(n):
      z_OK.append(C[0] + C[1]*x_OK[i] + C[2]*y_OK[i] + C[3]*x_OK[i]*y_OK[i] + C[4]*(x_OK[i]**2) + C[5]*(y_OK[i]**2) + C[6]*(x_OK[i]*(y_OK[i]**2)) + C[7]*((x_OK[i]**2)*y_OK[i]) + C[8]*(x_OK[i]**3) + C[9]*(y_OK[i]**3))

elif order == 4:
   for i in range(n):
      z_OK.append(C[0] + C[1]*x_OK[i] + C[2]*y_OK[i] + C[3]*x_OK[i]*y_OK[i] + C[4]*(x_OK[i]**2) + C[5]*(y_OK[i]**2) + C[6]*(x_OK[i]*(y_OK[i]**2)) + C[7]*((x_OK[i]**2)*y_OK[i]) + C[8]*(x_OK[i]**3) + C[9]*(y_OK[i]**3) + C[10]*(x_OK[i]*(y_OK[i]**3)) + C[11]*((x_OK[i]**3)*y_OK[i]) + C[12]*((x_OK[i]**2)*(y_OK[i]**2)) + C[13]*(x_OK[i]**4) + C[14]*(y_OK[i]**4))


value = []

for i in range(n):
   value.append(((f_OK[i] - z_OK[i])**2)/(n*(e_OK[i]**2)))


#print(value)
chi2 = sum(value)/(n-m)
print(chi2)



 














