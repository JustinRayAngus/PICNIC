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

}

#include "NamespaceFooter.H"
#endif

