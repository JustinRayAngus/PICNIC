#ifndef _MATHUTILS_H_
#define _MATHUTILS_H_

//#include "Dimensions.H"
//#include "BoundaryBoxLayout.H.multidim"

#include "Box.H"
#include "LevelData.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"

#include "math.h"

#include "NamespaceHeader.H"

namespace MathUtils {
  
   double errorfun( const double&  a_x );
   
   double errorinv( const double&  a_x );

   double rand();
 
   int randInt(const int& A, const int& B);

}

#include "NamespaceFooter.H"
#endif

