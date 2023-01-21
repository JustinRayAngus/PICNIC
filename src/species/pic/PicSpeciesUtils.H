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
   
   void applyForcesNew( List<JustinsParticle>&  a_pList,
               const Real                    a_fnorm,
               const Real                    a_cnormDt, 
               const bool                    a_byWholeDt,
               const bool                    a_anticyclic );
   
   void applyForcesAxisymm( List<JustinsParticle>&  a_pList,
                      const Real                    a_fnorm,
                      const Real                    a_cnormDt, 
                      const bool                    a_byHalfDt,
                      const bool                    a_anticyclic,
                      const int                     a_iter_max );

   void stepNormTransfer( List<JustinsParticle>&  a_pList,
                          List<JustinsParticle>&  a_temp_pList,
                    const RealVect&               a_dX,
                    const Real                    a_cnormHalfDt, 
                    const Real                    a_rtol, 
                    const bool                    a_reverse );

   void stepNormTransferInflow( List<JustinsParticle>&  a_pList,
                                List<JustinsParticle>&  a_temp_pList,
                          const RealVect&               a_dX,
                          const Real                    a_Xbdry, 
                          const int                     a_bdry_dir, 
                          const Real                    a_cnormDt, 
                          const Real                    a_rtol, 
                          const bool                    a_reverse );

}

#include "NamespaceFooter.H"
#endif

