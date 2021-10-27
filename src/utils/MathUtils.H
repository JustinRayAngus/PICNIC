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
//#include "JustinsParticlePtr.H"

#include "math.h"
#include <random>

#include "NamespaceHeader.H"
   
extern std::mt19937 global_rand_gen;
   
// This (using extern) may not be the best way to do this. 
// Does this mean that each file that includes the MathUtils.H 
// file will have its own local static random number engine?

namespace MathUtils {

   double errorfun( const double&  a_x );
   
   double errorinv( const double&  a_x );
   
   double gammainc( const Real  a_x, 
                    const Real  a_a );

   double rand();
   
   int randInt(const int& A, const int& B);
   
   double randn();

   void seedRNG( const int&  a_seed );

}

#include "NamespaceFooter.H"
#endif

