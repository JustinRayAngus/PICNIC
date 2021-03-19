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

#include "NamespaceHeader.H"
   
extern std::mt19937 global_rand_gen;
   
// this is not best way to do this. This means
// that each file that has the MathUtils.H header file
// will have its own local static random number engine.
// But I can't figure out how to do it better yet.
//static std::random_device global_rd;
//static std::mt19937 global_rand_gen(global_rd());

namespace MathUtils {

   double errorfun( const double&  a_x );
   
   double errorinv( const double&  a_x );
   
   double gammainc( const Real  a_x, 
                    const Real  a_a );

   double rand();
   
   int randInt(const int& A, const int& B);
   
   double randn();

}

#include "NamespaceFooter.H"
#endif

