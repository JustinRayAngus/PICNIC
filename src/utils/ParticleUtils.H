#ifndef _PARTICLEUTILS_H_
#define _PARTICLEUTILS_H_

#include "RealVect.H"

#include "ParticleData.H"
#include "Particle.H"

#include "ParticleUtilsF_F.H"

#include "NamespaceHeader.H"

namespace ParticleUtils {

   inline      
   void borisPusher( std::array<Real,3>&  a_vp,
                              const std::array<Real,3>&  a_vpold,
                              const std::array<Real,3>&  a_Ep,
                              const std::array<Real,3>&  a_Bp,
                              const Real&                a_fnorm_const,
                              const Real&                a_cnormDt ) 
   {
      CH_TIME("ParticleUtils::borisPusher()");

      FORT_BORIS_PUSHER( CHF_REAL(a_vp[0]),
                         CHF_REAL(a_vp[1]),
                         CHF_REAL(a_vp[2]),
                         CHF_CONST_REAL(a_vpold[0]),
                         CHF_CONST_REAL(a_vpold[1]),
                         CHF_CONST_REAL(a_vpold[2]),
                         CHF_CONST_REAL(a_Ep[0]),
                         CHF_CONST_REAL(a_Ep[1]),
                         CHF_CONST_REAL(a_Ep[2]),
                         CHF_CONST_REAL(a_Bp[0]),
                         CHF_CONST_REAL(a_Bp[1]),
                         CHF_CONST_REAL(a_Bp[2]),
                         CHF_CONST_REAL(a_fnorm_const),
                         CHF_CONST_REAL(a_cnormDt) );

   }

   inline
   void computeDeltaU( std::array<Real,3>&  a_deltaU,
                 const Real&                ux, 
                 const Real&                uy, 
                 const Real&                uz, 
                 const Real&                costh,
                 const Real&                sinth,
                 const Real&                cosphi,
                 const Real&                sinphi )
   {
      CH_TIME("ParticleUtils::computeDeltaU() v3");

      // define u and uperp
      Real u = sqrt(ux*ux + uy*uy + uz*uz);
      Real uperp = sqrt(ux*ux + uy*uy);

      // define deltaU
      a_deltaU[0] = ux*uz/uperp*sinth*cosphi - uy*u/uperp*sinth*sinphi - ux*(1.-costh);   
      a_deltaU[1] = uy*uz/uperp*sinth*cosphi + ux*u/uperp*sinth*sinphi - uy*(1.-costh);   
      a_deltaU[2] = -uperp*sinth*cosphi - uz*(1.-costh);   

      /*
      FORT_COMPUTE_DELTAU( CHF_REAL(a_deltaU[0]),
                           CHF_REAL(a_deltaU[1]),
                           CHF_REAL(a_deltaU[2]),
                           CHF_CONST_REAL(ux),
                           CHF_CONST_REAL(uy),
                           CHF_CONST_REAL(uz),
                           CHF_CONST_REAL(costh),
                           CHF_CONST_REAL(sinth),
                           CHF_CONST_REAL(cosphi),
                           CHF_CONST_REAL(sinphi) );
      */

   }
  
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

}

#include "NamespaceFooter.H"
#endif

