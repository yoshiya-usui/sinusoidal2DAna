# Function
This program calculates analytical MT apparent resistivity, phasae, and vertical magnetic transfer functions for a 2-D sinusoidal interface.
The algorithm is based on the following studies.

*K. Schwalenberg, R. N. Edwards, The effect of seafloor topography on magnetotelluric fields: an analytical formulation confirmed with numerical results, Geophysical Journal International, Volume 159, Issue 2, November 2004, Pages 607–621, https://doi.org/10.1111/j.1365-246X.2004.02280.x*

*Yoshiya Usui, Takafumi Kasaya, Yasuo Ogawa, Hisanori Iwamoto, Marine magnetotelluric inversion with an unstructured tetrahedral mesh, Geophysical Journal International, Volume 214, Issue 2, August 2018, Pages 952–974, https://doi.org/10.1093/gji/ggy171*

Originally, Schwalenberg and Edwards (2004) proposed the formulation, and later, Usui et al. (2018) modified the formulation for the TM-mode.

# How to use
1) Conductivities, amplitude, and wavelength are specified by sigma1, sigma2, delta, and lambda in "cal.m" and those in calcBetaGammaTE and calcBetaGammaTM of "sinusoidal2DAna.cpp".
2) Calculate modified Bessel function by a MATLAB program "cal.m". Periods (sec) should be specified by variable T in "cal.m".
3) Calculate analytical response functions by a C++ program "sinusoidal2DAna.cpp". In executing the program, you need to specify the file created by "cal.m", the file name in which location of observation points are written, and the period(sec) by its arguments as follows.

   *sinusoidal2DAna  "File name of modified Bessel function" "File name of observation points" "Period(sec)"*

# File format
  **Input file about observation points (arbitrary name)**

  *[Number of observation points]*

  *[X-coordinate(m) of the 1st point]*

  *[X-coordinate(m) of the 2nd point]*

  *[X-coordinate(m) of the 3rd point]*

  ...
  
  **Output file for the TE-mode (TE.txt)**
  
  *[X-coordinate(m) of the 1st point] [Z-coordinate(m) of the 1st point] [apparent resistivity(Ohm-m) of the upper region] [phase(deg.) of the upper region] [apparent resistivity(Ohm-m) of the lower region] [phase(deg.) of the lower region] [Re(Tzx) of the upper region] [Im(Tzx) of the upper region] [Re(Tzx) of the lower region] [Im(Tzx) of the lower region]*
  
  *[X-coordinate(m) of the 2nd point] [Z-coordinate(m) of the 2nd point] [apparent resistivity(Ohm-m) of the upper region] [phase(deg.) of the upper region] [apparent resistivity(Ohm-m) of the lower region] [phase(deg.) of the lower region] [Re(Tzx) of the upper region] [Im(Tzx) of the upper region] [Re(Tzx) of the lower region] [Im(Tzx) of the lower region]*
  
  *[X-coordinate(m) of the 3rd point] [Z-coordinate(m) of the 3rd point] [apparent resistivity(Ohm-m) of the upper region] [phase(deg.) of the upper region] [apparent resistivity(Ohm-m) of the lower region] [phase(deg.) of the lower region] [Re(Tzx) of the upper region] [Im(Tzx) of the upper region] [Re(Tzx) of the lower region] [Im(Tzx) of the lower region]*
  
  ...
    
  **Output file for the TM-mode (TM.txt)**
  
  *[X-coordinate(m) of the 1st point] [Z-coordinate(m) of the 1st point] [apparent resistivity(Ohm-m) of the upper region] [phase(deg.) of the upper region] [apparent resistivity(Ohm-m) of the lower region] [phase(deg.) of the lower region]*
  
  *[X-coordinate(m) of the 2nd point] [Z-coordinate(m) of the 2nd point] [apparent resistivity(Ohm-m) of the upper region] [phase(deg.) of the upper region] [apparent resistivity(Ohm-m) of the lower region] [phase(deg.) of the lower region]*
  
  *[X-coordinate(m) of the 3rd point] [Z-coordinate(m) of the 3rd point] [apparent resistivity(Ohm-m) of the upper region] [phase(deg.) of the upper region] [apparent resistivity(Ohm-m) of the lower region] [phase(deg.) of the lower region]*
 
  ...
  
