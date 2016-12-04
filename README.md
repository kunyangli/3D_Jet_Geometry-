# 3D_Jet_Geometry-
Last Edited: Dec 2016
Questions: Contact Katie Kosak, katie.kosak@gmail.com

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

To execute the code, execute the Jet_Curry_Main.py code as a normal Python code. If the Jet_Curry.pyc file is in the same directory as Jet_Curry_Main, the functions will be automatically imported. 

where jet tells the code to pull the function from the code Jet_Curry. Find_MaxFlux is the function in Jet_Curry being called. 

For the functions in Jet_Curry, all arguments must be included as the example. 

The following must be defined:

1) The core location is the point where the stream originates in the code. Not necessarily the origin. 

2) Upstream is the core points (by default). However, it still needs to be defined just in case the user didn't want the core in the geometry calculations.

3) Downstream is the point where the code should stop. 

4) Theta, angle of line of sight. This value must be entered in radians. 

5) Fits Image must be a standard NASA FITS file format to be read by astropy. 

----- Visual Python-----

Visual Python can be used to view the 3D structure in a graphics IDE.

The installation of VPython can be found in detail here: vpython.org

The code should be run from the VPython IDE and will automatically read in the x,y,z from the Cartesian_Coordinates file.

With the line of sight, visulatizations, the code must be manually changed by the user. If you have questions about the use of our Vpython code, please ask.
