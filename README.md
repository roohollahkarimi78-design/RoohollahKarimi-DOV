# User guide
The folder of "Conventional_Exact_xi_eta" contain all the MATLAB functions required to compute the components of the Helmert deflection of the vertical (DOV) 
from global geopotential models (GGMs) using both the conventional method and the exact (proposed) method.

The main functions are "Conventional_xi_eta.m" and "Exact_xi_eta.m" that compute the north-south, xi, and east-west, eta, components of the Helmert DOV 
based on the geodetic coordinates of the computation point.

Descriptions for each function, along with all required parameters, including inputs and outputs, are provided within the function itself. 

The computations are based on the EIGEN-6C4 GGM and the GRS80 reference ellipsoid. However, the GGM and the reference ellipsoid can be changed.
The files of "Cnm_EIGEN_6C4_zero_tide.mat" and "Snm_EIGEN_6C4_zero_tide.mat" contain the Cnm and Snm coefficients for the EIGEN-6C4 model in the zero-tide system.

The dataset files used in the study regions, named "IRN.mat", "USA.mat", "ALP.mat", and "AUS.mat", are located in the "Conventional_Exact_xi_eta" folder. 

The format of the dataset files is as follows:

Column 1: Point number

Column 2: Geodetic latitude (in degrees)

Column 3: Geodetic longitude (in degrees)

Column 4: Geodetic height (in meters)

The programs of "Prog_Conventional_Exact_xi_eta_IRN.m", "Prog_Conventional_Exact_xi_eta_USA.m", "Prog_Conventional_Exact_xi_eta_ALP.m", and "Prog_Conventional_Exact_xi_eta_AUS.m" 
compute the components of the Helmert DOV for the study regions. The outputs are saved in the files named "Conventional_Exact_xi_eta_IRN.mat", "Conventional_Exact_xi_eta_USA.mat", 
"Conventional_Exact_xi_eta_ALP.mat", and "Conventional_Exact_xi_eta_AUS.mat".

The format of the output files is as follows:

Column 1: Point number

Column 2: Geodetic latitude (in degrees)

Column 3: Geodetic longitude (in degrees)

Column 4: Geodetic height (in meters)

Column 5: xi from the conventional method (in arcseconds)

Column 6: eta from the conventional method (in arcseconds)

Column 7: xi from the exact method (in arcseconds)

Column 8: eta from the exact method (in arcseconds)

Important note:
Due to computational limitations in the MATLAB software regarding the calculation of Legendre functions at high degrees for specific geodetic latitudes, 
degrees exceeding 2065 should not be used within the geodetic latitude range of 57° 30′ to 60° 15′.

If you have any questions regarding functions, programs, or data, please feel free to contact the corresponding author.
E-mail: karimi@tafreshu.ac.ir, roohollah.karimi.78@gmail.com
 
