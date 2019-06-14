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

#from reduced illumination images

C1 = [3.71567046e+05, -9.89670137e-01, -4.36761152e+00]

C2 = [3.69874645e+05, -7.56598474e+00, 6.40014321e+00, -2.22394262e-04, 3.50300281e-03, -5.15523021e-03]

C3 = [3.62476746e+05, 3.08403121e+00, 3.95365567e+01, 4.71614609e-03, -1.63380948e-02, -4.36081134e-02, 4.56980354e-06, -8.09013366e-06, 9.58817141e-06, 1.11045770e-05]

C4 = [3.54782302e+05, -1.26081599e+01, 8.02355097e+01, -8.18679212e-02, 8.76330633e-02, -7.78307819e-02, 4.10779938e-05, 5.61981006e-05, -1.01574313e-04, 2.06315540e-05, -7.14647808e-09, -1.84659550e-08, 6.02417938e-09, 3.37452954e-08, -1.69484871e-10]

#from unreduced illumination images 
'''
C1 = [3.86763870e+05, -1.00704422e+01, -1.42214184e+01] 

C2 = [3.44618481e+05, 2.99232173e+01, 6.73486922e+01, -1.17148599e-03, -1.98171479e-02, -3.74318355e-02]


C3 = [3.55663437e+05, 2.39505267e+01, 4.08958194e+01, 3.95302465e-02, -4.51743514e-02, -2.45510200e-02, -9.64744766e-07,  -2.02177377e-05, 1.62483705e-05, -3.13011578e-06]

C4 = [3.15324661e+05, 8.48910543e+01, 2.16550518e+02, -1.68077457e-01, -1.61506660e-02, -2.49863073e-01, 1.10382380e-04, 9.53246445e-05, -6.26516471e-05, 1.15227683e-04, -2.63544280e-08, -3.00375014e-08, -1.11168390e-08, 2.89992386e-08, -2.12920170e-08]
'''


def get_XY(image):
   hdulist = fits.open(image)
   img_data = hdulist[0].data
   n1 = len(img_data[:,0])
   n2 = len(img_data[0,:])
   a1 = np.arange(n1)
   a2 = np.arange(n2)
   mesh = np.dstack(np.meshgrid(a1,a2)).reshape(-1,2).astype(float)
   x = mesh[:,0]
   y = mesh[:,1]
   return x,y,img_data

x,y,data = get_XY("reduced_00009646.M42.fit") #replace with the image you want to correct

n1, n2 = data.shape
print(n1, n2)
print(x.shape)
print(y.shape)

order = 2

if order == 1:
   Z = np.dot(np.c_[np.ones(x.shape), x, y], C1)
elif order == 2:
   Z = np.dot(np.c_[np.ones(x.shape), x, y, x*y, x**2, y**2], C2)
elif order == 3:
   Z = np.dot(np.c_[np.ones(x.shape), x, y, x*y, x**2, y**2, x*(y**2), (x**2)*y, x**3, y**3], C3)
elif order == 4:
   Z = np.dot(np.c_[np.ones(x.shape),x,y,x*y,x**2,y**2,x*(y**2),(x**2)*y,x**3,y**3,x*(y**3),(x**3)*y,(x**2)*(y**2),x**4,y**4], C4)


Z = Z.reshape((n2,n1))
Z /= Z.mean()
#plt.imshow(Z, origin = 'lower')


data = (data.T)/Z
data = data.T
hdu_corrected = fits.PrimaryHDU(data)
hdu_corrected.writeto('corrected.fit')

#fig = plt.figure()
#plt.contourf(x, y, Z, origin = 'lower')
#plt.imshow(data, origin = 'lower')
#plt.imshow(Z, origin = 'lower')
#ax = fig.gca(projection='3d')
#ax.plot_surface(X, Y, Z, alpha = 0.2)

#plt.xlabel('X')
#plt.ylabel('Y')
#ax.set_zlabel('Z')

plt.show()














