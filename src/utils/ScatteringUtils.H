#ifndef _SCATTERINGUTILS_H_
#define _SCATTERINGUTILS_H_

#include "RealVect.H"

#include <cmath>
#include "math.h"

#include "NamespaceHeader.H"

namespace ScatteringUtils {
  
   void computeDeltaU( std::array<Real,3>&  a_deltaU,
                 const std::array<Real,3>&  a_vp1,
                 const std::array<Real,3>&  a_vp2,
                 const Real&                a_theta );
   
   void computeDeltaU( std::array<Real,3>&  a_deltaU,
                 const Real&                ux, 
                 const Real&                uy, 
                 const Real&                uz, 
                 const Real&                a_theta );
   
}

#include "NamespaceFooter.H"
#endif

