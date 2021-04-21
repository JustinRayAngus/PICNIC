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
  
   void borisPusher( std::array<Real,3>&  a_vp,
               const std::array<Real,3>&  a_vpold,
               const std::array<Real,3>&  a_Ep,
               const std::array<Real,3>&  a_Bp,
               const Real&                a_fnorm_const,
               const Real&                a_cnormDt );

   void computeDeltaU( std::array<Real,3>&  a_deltaU, 
                 const std::array<Real,3>&  a_vp1,
                 const std::array<Real,3>&  a_vp2,
                 const Real&                a_theta,
                 const Real&                a_phi );

   void computeDeltaU( std::array<Real,3>&  a_deltaU,
                 const Real&                ux, 
                 const Real&                uy, 
                 const Real&                uz, 
                 const Real&                a_theta,
                 const Real&                a_phi );

   void computeDeltaU( std::array<Real,3>&  a_deltaU,
                 const Real&                ux, 
                 const Real&                uy, 
                 const Real&                uz, 
                 const Real&                costh,
                 const Real&                sinth,
                 const Real&                cosphi,
                 const Real&                sinphi );
      
}

#include "NamespaceFooter.H"
#endif

