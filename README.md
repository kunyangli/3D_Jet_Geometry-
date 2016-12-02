# 3D_Jet_Geometry-
Python code that models the 3D geometry of AGN jets

Make sure the files are in the same directory. 

Jet_Curry has the functions for performing the calculation from 2D to 3D. 
There is no need to change this code to perform calculations. It is best not to 
mess around with this code. 

Jet_Curry_Main is the code that the user will use to define the fits file, define the 
core location, etc. 

To call in the tools used from Jet_Curry, the following command is used:

import Jet_Curry as jet

To use the functions necessary for the calculations, use the following:

jet.Find_MaxFlux(file1,Upstream_Bounds,number_of_points)

where jet tells the code to pull the function from the code Jet_Curry. Find_MaxFlux is the function in Jet_Curry being called. 

For the functions in Jet_Curry, all arguments must be included as the example. 

The following must be defined:

1) The core location is the point where the stream originates in the code. 

2) Upstream is the core points. However, it still needs to be defined just in case the user didn't want the core in the geometry calculations.

3) Downstream is the point where the code should stop calculating the 3D Geometry. 

4) Theta. Line of sight. This value must be entered in radians. 
