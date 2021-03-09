#ifndef _SCATTERINGUTILS_H_
#define _SCATTERINGUTILS_H_

#include "RealVect.H"

#include <cmath>
#include "math.h"

#include "NamespaceHeader.H"

namespace ScatteringUtils {
  
   void computeDeltaU( std::array<Real,3>&  a_deltaU,
                 const Real  a_vx1,
                 const Real  a_vy1,
                 const Real  a_vz1,
                 const Real  a_vx2,
                 const Real  a_vy2,
                 const Real  a_vz2,
                 const Real  a_th );
   
}

#include "NamespaceFooter.H"
#endif

