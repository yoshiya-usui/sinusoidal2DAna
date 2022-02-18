# Function
This program to calculate analytical appanret resistivity, phasae, and vertical magnetic transfer functions for a sinusoidal interface.
The algorithm is based on the following studies.

*K. Schwalenberg, R. N. Edwards, The effect of seafloor topography on magnetotelluric fields: an analytical formulation confirmed with numerical results, Geophysical Journal International, Volume 159, Issue 2, November 2004, Pages 607–621, https://doi.org/10.1111/j.1365-246X.2004.02280.x*

*Yoshiya Usui, Takafumi Kasaya, Yasuo Ogawa, Hisanori Iwamoto, Marine magnetotelluric inversion with an unstructured tetrahedral mesh, Geophysical Journal International, Volume 214, Issue 2, August 2018, Pages 952–974, https://doi.org/10.1093/gji/ggy171*

Originally, Schwalenberg and Edwards (2004) proposed the formulation, and later, Usui et al. (2018) modified the formulation for the TM-mode.

# How to use
1) Conductivities, amplitude, and wavelength are specified by sigma1, sigma2, delta, and lambda in "cal.m" and those in calcBetaGammaTE and calcBetaGammaTM of "sinusoidal2DAna.cpp".
2) Calculate modified Bessel function by a MATLAB program "cal.m"
3) Calculate analytical response functions b a C++ program "sinusoidal2DAna.cpp". In executing the program, you need to specify the file created by "cal.m", the file name in which location of observation points are written, and the period(sec) by its arguments as follows.

   *sinusoidal2DAna  "File name of modified Bessel function" "File name of observation points" "Period(sec)"*

# File format
  **File for observation points** 
  Number of observation points
  Coordinate (m) of the 1st point
  Coordinate (m) of the 2nd point
  Coordinate (m) of the 3rd point
  ...
  
 
