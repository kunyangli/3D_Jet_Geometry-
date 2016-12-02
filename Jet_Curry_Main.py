__author__="Katie Kosak"
__email__="katie.kosak@gmail.com"
__credits__=["KunYang Li", "Dr. Eric Perlman", "Dr.Sayali Avachat"]
__status__="production"
from scipy import array
import pyfits
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import *
import Jet_Curry as jet
## Open the image from fits file 
## Obtain image parameters
np.seterr(all='ignore')
s=[]
eta=[]
theta=0.261799388 # Line of Sight (radians)
core_points=np.array([20.5,16]) # Position of Nucleus on Image
nameofknots='KnotD_Radio'
filename=nameofknots+'.fits'
Downstream_Bounds=np.array([80,10])
Upstream_Bounds=np.array([20.5,16])
number_of_points=Downstream_Bounds[0]-Upstream_Bounds[0]

nameofknots='KnotD_Radio'
filename=nameofknots+'.fits'

##  Obtain Information from the Image/Show the Image
file1=pyfits.getdata(filename)
pixel_min=np.nanmin(file1)
pixel_maxima=np.nanmax(file1)
file1=jet.imagesqrt(file1,pixel_min,pixel_maxima) # Square Root Scaling for fits image
plt.imshow(file1)

## Go column by column to calculate the max flux for each column
x,y,x_smooth,y_smooth,intensity_max=jet.Find_MaxFlux(file1,Upstream_Bounds,number_of_points)

## Plot the Max Flux over the Image
plt.contour(file1,10, cmap='gray')
#plt.scatter(x_smooth,y_smooth,c='b')
plt.scatter(x,y,c='b')
plt.title('M87:Knot D Outline of Jet Stream')
ax = plt.gca()
ax.invert_yaxis()
plt.show()

## Calculate the s and eta values
## s,eta, x_smooth,y_smooth values will 
## be stored in parameters.txt file
s,eta=jet.Calculate_s_and_eta(x_smooth,y_smooth,core_points)
### Run the First MCMC Trial in Parallel
jet.MCMC1_Parallel(s,eta,theta)
jet.MCMC2_Parallel(s,eta,theta)
jet.Annealing1_Parallel(s,eta,theta)
jet.Annealing2_Parallel(s,eta,theta) 
x_coordinates,y_coordinates,z_coordinates=jet.Convert_Results_Cartesian(s,eta,theta)

# Plot the Results on Image   
#plt.scatter(x_coordinates,y_coordinates,c='y')
#plt.scatter(x_smooth,y_smooth,c='r')
#ax = plt.gca()
#ax.invert_yaxis()
#plt.show()