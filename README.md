# Function
This program to calculate analytical appanret resistivity, phasae, and vertical magnetic transfer functions for a sinusoidal interface.
The algorithm is based on the following studies.

*K. Schwalenberg, and R.N. Edwards. The effect of seafloor topography on magnetotelluric fields: an analytical formulation confirmed with
numerical results, Geophys. J. Int., 159(2), 607–621, 2004.*

*Y. Usui, T. Kasaya, Y. Ogawa and H. Iwamoto, Marine magnetotelluric inversion with an unstructured tetrahedral mesh, Geophys. J. Int., 214(2): 952-974, https://doi.org/10.1093/gji/ggy171, 2018.*

Originally, Schwalenberg and Edwards (2004) proposed the formulation, and later, Usui et al. (2018) modified the formulation for the TM-mode.

# How to use
1) First, modified Bessel function is calculated with MATLAB program "cal.m"
2) Next, analytical response functions are calculated with C++ program "sinusoidal2DAna.cpp". In executing the program, you need to specify the file created by "cal.m", the file name in which location of observation points are written, and the period(sec) by its arguments as follows.
   *sinusoidal2DAna  "File name of modified Bessel function" "File name of observation points" "Period(sec)"*

# File format
Under construction...
