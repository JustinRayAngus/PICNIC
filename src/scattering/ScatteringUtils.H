#ifndef _SCATTERINGUTILS_H_
#define _SCATTERINGUTILS_H_

#include "RealVect.H"

#include "NamespaceHeader.H"

namespace ScatteringUtils {

   inline
   Real getScatteringCos( const Real  a_R,
                          const Real  a_xi )
   {
      Real costh = 1.0 - 2.0*a_R*(1.0 - a_xi)/(1.0 + a_xi*(1.0 - 2.0*a_R));
      return costh;
   }
   
   inline
   void collapseThreeToTwo( std::array<Real,3>&  vp2,
                            Real&                wp2, 
                            std::array<Real,3>&  vp3,
                            Real&                wp3, 
                      const std::array<Real,3>&  vp2p,
                      const Real                 wp2p ) 
   {

      // transform particles 2, 2p and 3 into two equally-weighted particles such
      // that momentum, and direction-dependent energy are conserved
      CH_assert(wp2>wp2p);

      Real wp23 = 0.5*(wp2 + wp3);
      Real c23, d23, arg23;
      for (int dir=0; dir<3; dir++) {
         c23 = ( wp2p*vp2p[dir] + (wp2-wp2p)*vp2[dir] + wp3*vp3[dir] )/wp23;
         d23 = ( wp2p*vp2p[dir]*vp2p[dir] + (wp2-wp2p)*vp2[dir]*vp2[dir] 
             +   wp3*vp3[dir]*vp3[dir] )/wp23;
         //
         arg23 = 2.0*d23 - c23*c23;
         CH_assert(arg23>=0.0);
         vp2[dir] = 0.5*(c23 + sqrt(arg23));
         vp3[dir] = 0.5*(c23 - sqrt(arg23));
      }
      wp2 = wp23;
      wp3 = wp23;

   }

   inline
   void computeDeltaU( std::array<Real,3>&  a_deltaU,
                 const Real                 ux, 
                 const Real                 uy, 
                 const Real                 uz, 
                 const Real                 costh,
                 const Real                 sinth,
                 const Real                 cosphi,
                 const Real                 sinphi )
   {

      // define u and uperp
      Real u = sqrt(ux*ux + uy*uy + uz*uz);
      Real uperp = sqrt(ux*ux + uy*uy);

      // define deltaU
      if(uperp==0.0) {
         a_deltaU[0] = u*sinth*cosphi;   
         a_deltaU[1] = u*sinth*sinphi;   
         a_deltaU[2] = u*costh - u;   
      }
      else {
         a_deltaU[0] = ux*uz/uperp*sinth*cosphi - uy*u/uperp*sinth*sinphi - ux*(1.-costh);   
         a_deltaU[1] = uy*uz/uperp*sinth*cosphi + ux*u/uperp*sinth*sinphi - uy*(1.-costh);   
         a_deltaU[2] = -uperp*sinth*cosphi - uz*(1.-costh);   
      }

   }
   
   void computeDeltaU( std::array<Real,3>&  a_deltaU, 
                 const std::array<Real,3>&  a_vp1,
                 const std::array<Real,3>&  a_vp2,
                 const Real                 a_costh,
                 const Real                 a_sinth,
                 const Real                 a_phi );
   
  void computeDeltaU( std::array<Real,3>&  a_deltaU,
                const Real                 ux, 
                const Real                 uy, 
                const Real                 uz, 
                const Real                 a_costh,
                const Real                 a_sinth,
                const Real                 a_phi );
   
   void computeDeltaU( std::array<Real,3>&  a_deltaU, 
                 const std::array<Real,3>&  a_vp1,
                 const std::array<Real,3>&  a_vp2,
                 const Real                 a_potential,
                 const Real                 a_mu,
                 const Real                 a_costh,
                 const Real                 a_sinth,
                 const Real                 a_phi );

   Real linearInterp( int&                a_index,
                const std::vector<Real>&  a_X,
                const std::vector<Real>&  a_Y,
                const Real                a_X0 );
   
   Real linearInterp( const std::vector<Real>&  a_X,
                      const std::vector<Real>&  a_Y,
                      const Real                a_X0,
                      const int                 a_index );
   
   Real semilogInterp( const std::vector<Real>&  a_X,
                       const std::vector<Real>&  a_Y,
                       const Real                a_X0,
                       const int                 a_index );  
   
   Real loglogInterp( const std::vector<Real>&  a_X,
                      const std::vector<Real>&  a_Y,
                      const Real                a_X0,
                      const int                 a_index );  
   
}

#include "NamespaceFooter.H"
#endif

