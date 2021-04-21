#ifndef _PARTICLEUTILS_H_
#define _PARTICLEUTILS_H_

//#include "Dimensions.H"
//#include "BoundaryBoxLayout.H.multidim"

#include "Box.H"
#include "LevelData.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"

#include "ParticleData.H"
#include "Particle.H"

#include "NamespaceHeader.H"

namespace ParticleUtils {
  
   void moment( LevelData<FArrayBox>&    a_momentum,
          const ParticleData<Particle>&  a_Pdata );
  
   void moment( FArrayBox&  a_momentum,
          const Particle&   a_Pdata );
   
   void borisPusher( std::array<Real,3>&  a_vp,
               const std::array<Real,3>&  a_vpold,
               const std::array<Real,3>&  a_Ep,
               const std::array<Real,3>&  a_Bp,
               const Real&                a_fnorm_const,
               const Real&                a_cnormDt );
      
}

#include "NamespaceFooter.H"
#endif

