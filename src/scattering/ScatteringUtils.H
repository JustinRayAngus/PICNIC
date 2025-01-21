#ifndef _SCATTERINGUTILS_H_
#define _SCATTERINGUTILS_H_

#include "RealVect.H"
#include <array>

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
   void rotateVelocity( std::array<Real,3>&  a_u,
                  const Real                 costh,
                  const Real                 sinth,
                  const Real                 cosphi,
                  const Real                 sinphi )
   {

      // define u and uperp
      Real ux = a_u[0];
      Real uy = a_u[1];
      Real uz = a_u[2];
      Real u = sqrt(ux*ux + uy*uy + uz*uz);
      Real uperp = sqrt(ux*ux + uy*uy);

      // rotate
      if(uperp==0.0) {
         a_u[0] = u*sinth*cosphi;   
         a_u[1] = u*sinth*sinphi;   
         a_u[2] = u*costh;   
      }
      else {
         a_u[0] = ux*uz/uperp*sinth*cosphi - uy*u/uperp*sinth*sinphi + ux*costh;   
         a_u[1] = uy*uz/uperp*sinth*cosphi + ux*u/uperp*sinth*sinphi + uy*costh;   
         a_u[2] = -uperp*sinth*cosphi + uz*costh;   
      }

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

   inline
   void modEnergyPairwise( std::array<Real,3>&  a_betap1,
                           std::array<Real,3>&  a_betap2,
                     const long double          a_wpmp1,
                     const long double          a_wpmp2,
                          const double          a_Erel_frac,
                                double&         a_Erel_cumm,
                           long double&         a_deltaE )
   {
      // method to add/subtract a small percent of energy from a
      // pair of particles such that momentum is not disturbed
      //
      // This is done by treating it as an inelastic scattering event
      // with potential U = a_deltaE and zero scattering angle
      // a_deltaE > 0 ==> need to take energy away
      // a_deltaE < 0 ==> need to add energy

      int sign = 1;
      if (a_deltaE<0.0) { sign = -1; }

      const Real ux = a_betap1[0] - a_betap2[0];
      const Real uy = a_betap1[1] - a_betap2[1];
      const Real uz = a_betap1[2] - a_betap2[2];
      const long double usq = ux*ux + uy*uy + uz*uz;
#ifdef RELATIVISTIC_PARTICLES
      long double gbsq1 = 0.0;
      for (int n=0; n<3; n++) { gbsq1 += a_betap1[n]*a_betap1[n]; }
      long double gbsq2 = 0.0;
      for (int n=0; n<3; n++) { gbsq2 += a_betap2[n]*a_betap2[n]; }
      const long double gamma1 = std::sqrt(1.0 + gbsq1);
      const long double gamma2 = std::sqrt(1.0 + gbsq2);
      const long double E1 = a_wpmp1*gamma1;
      const long double E2 = a_wpmp2*gamma2;
      const long double Etot = E1 + E2;
      const long double pxtot = a_wpmp1*a_betap1[0] + a_wpmp2*a_betap2[0];
      const long double pytot = a_wpmp1*a_betap1[1] + a_wpmp2*a_betap2[1];
      const long double pztot = a_wpmp1*a_betap1[2] + a_wpmp2*a_betap2[2];
      const long double Ecm = std::sqrt(Etot*Etot - pxtot*pxtot - pytot*pytot - pztot*pztot);
      const long double Erel = Ecm - a_wpmp1 - a_wpmp2;
#else
      const long double muR = a_wpmp1*a_wpmp2/(a_wpmp1 + a_wpmp2);
      const long double Erel = muR/2.0*usq;
#endif
      if (Erel<=0.0) { return; }

      long double deltaE = sign*a_Erel_frac*Erel;
      if (std::abs(deltaE) > std::abs(a_deltaE)) {
         deltaE = a_deltaE;
         a_deltaE = 0.0;
      }
      else { a_deltaE -= deltaE; }
      a_Erel_cumm += Erel - deltaE;

#ifdef RELATIVISTIC_PARTICLES
      const long double A = Etot - deltaE;
      const long double D = A*A + E2*E2 - E1*E1;
      const long double p2dotu = a_wpmp2*(a_betap2[0]*ux + a_betap2[1]*uy + a_betap2[2]*uz);
      const long double ptdotu = pxtot*ux + pytot*uy + pztot*uz;

      // compute coefficients for quadratic equation for alpha
      const long double a = A*A*usq - ptdotu*ptdotu;
      const long double b = D*ptdotu - 2*A*A*p2dotu;
      const long double c = A*A*E2*E2 - D*D/4.0;

      const long double root = b*b - 4.0*a*c;
      if (root < 0.0 || a == 0.0) { return; }
      const long double alpha = (-b + std::sqrt(root))/(2.0*a);

      // update particle velocities
      const long double ratio1 = alpha/a_wpmp1;
      const long double ratio2 = alpha/a_wpmp2;
      a_betap1[0] += ratio1*ux;
      a_betap1[1] += ratio1*uy;
      a_betap1[2] += ratio1*uz;
      a_betap2[0] -= ratio2*ux;
      a_betap2[1] -= ratio2*uy;
      a_betap2[2] -= ratio2*uz;
#else
      const long double uprime_over_u = std::sqrt(1.0 - deltaE/Erel);
      std::array<Real,3> deltaU;
      deltaU[0] = uprime_over_u*ux - ux;
      deltaU[1] = uprime_over_u*uy - uy;
      deltaU[2] = uprime_over_u*uz - uz;

      for (int n=0; n<3; n++) {
         a_betap1[n] += muR/a_wpmp1*deltaU[n];
         a_betap2[n] -= muR/a_wpmp2*deltaU[n];
      }
#endif

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

   void LorentzTransform( Real&                gammapst,
                          std::array<Real,3>&  upst,
                    const Real                 gammap,
                    const std::array<Real,3>&  up,
                    const Real                 gammacm,
                    const std::array<Real,3>&  ucm,
                    const bool                 reverse = false );

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

