#ifndef _PicSpeciesUtils_H_
#define _PicSpeciesUtils_H_

#include "Box.H"
#include "ListBox.H"
#include "JustinsParticle.H"
#include "DomainGrid.H"

#include "NamespaceHeader.H"
   
namespace PicSpeciesUtils {

   void applyForces( List<JustinsParticle>&  a_pList,
               const Real                    a_fnorm,
               const Real                    a_cnormDt, 
               const bool                    a_byHalfDt,
               const bool                    a_anticyclic );
   
   void applyForces_CYL_CYL( List<JustinsParticle>&  a_pList,
                       const Real                    a_fnorm,
                       const Real                    a_cnormDt, 
                       const bool                    a_byHalfDt,
                       const bool                    a_anticyclic );
   
   void applyForces_SPH_SPH( List<JustinsParticle>&  a_pList,
                       const Real                    a_fnorm,
                       const Real                    a_cnormDt, 
                       const bool                    a_byHalfDt );
   
   void applyForces_CYL_HYB( List<JustinsParticle>&  a_pList,
                       const Real                    a_fnorm,
                       const Real                    a_cnormDt, 
                       const bool                    a_anticyclic );
   
   void applyForces_SPH_HYB( List<JustinsParticle>&  a_pList,
                       const Real                    a_fnorm,
                       const Real                    a_cnormDt );

#ifdef RELATIVISTIC_PARTICLES
   inline Real getImplicitGamma( const std::array<Real,3>&  upold,
                                 const std::array<Real,3>&  upnew )
   {
      // return gamma factor used to convert proper velocity to velocity for
      // implicit advance methods
      Real gbsq_old = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
      Real gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
      Real gammap = 0.5*(std::sqrt(1.0 + gbsq_old) + std::sqrt(1.0 + gbsq_new));
      return gammap;
   }
#endif

}

#include "NamespaceFooter.H"
#endif

