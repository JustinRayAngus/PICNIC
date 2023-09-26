#ifndef _PicnicConstants_H_
#define _PicnicConstants_H_

#include <cmath>

#include "NamespaceHeader.H"

/**
   A namespace providing definitions of commomly used physical constants
   and conversion factors.
 */
namespace Constants {
   
   // Mathematical constants
   const double PI    = M_PI;
   const double HALFPI = M_PI / 2.0;
   const double TWOPI  = 2.0 * M_PI;
   const double FOURPI = 4.0 * M_PI;
   const double E     = exp( 1.0 );

   // Fundamental Physical Constants
   const double CVAC   = 2.99792458e+08;  // speed of light [m/s]
   const double MU0    = FOURPI*1.0e-7;   // vacuum permeability [H/m]
   const double EP0    = 1.0/CVAC/CVAC/MU0;  // vacuum permittivity [C^2/N-m^2]
   const double ME     = 9.10938370e-31;  // electron mass [kg]
   const double MN     = 1.67492750e-27;  // neutron mass [kg]
   const double MP     = 1.67262192e-27;  // proton mass [kg]
   const double AMU    = 1.66053901e-27;  // atomic mass unit [kg]
   const double QE     = 1.60217663e-19;  // elementary charge [C]
   const double KB     = 1.38064852e-23;  // Boltzmann constant [J/K]
   const double H      = 6.62607015e-34;  // Planck's constant [J-s]
   
   // Derived fundamental constants
   const double HBAR   = H/TWOPI;                            // reduced Planck's constant [J-s]
   const double ABOHR  = 2.0*TWOPI*EP0*HBAR*HBAR/(ME*QE*QE); // Bohr radius [m] 

   // Conversion Factors
   const double JOULE_PER_EV      = QE;
   const double EV_PER_JOULE      = 1.0 / JOULE_PER_EV;
   const double K_PER_EV          = QE/KB;
   const double EV_PER_K          = KB/QE;
   const double CM_PER_M          = 1.0e+02;
   const double CM2_PER_M2        = 1.0e+04;
   const double CM3_PER_M3        = 1.0e+06;
   const double M_PER_CM          = 1.0/CM_PER_M;
   const double M2_PER_CM2        = 1.0/CM2_PER_M2;
   const double M3_PER_CM3        = 1.0/CM3_PER_M3;
   const double SECOND_PER_HOUR   = 3.6e+03;
   const double SECOND_PER_MINUTE = 6.0e+01;

}

#include "NamespaceFooter.H"

#endif
