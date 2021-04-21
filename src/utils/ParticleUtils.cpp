#include "ParticleUtils.H"
#include "ParticleUtilsF_F.H"

#include "NamespaceHeader.H"
      
void ParticleUtils::moment( LevelData<FArrayBox>&    a_moment,
                      const ParticleData<Particle>&  a_Pdata )
{
   /*
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FluxBox& this_face_phi( a_face_phi[dit] );
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FluxBox& this_norm_vel( a_norm_vel[dit] );
      const Box& this_dbl_box( grids[dit] ); // this box has no ghost cells
      upWindToFaces( this_face_phi, this_cell_phi, this_norm_vel, this_dbl_box, a_method );
       
   } 
   */

}

void ParticleUtils::moment( FArrayBox&  a_moment,
                      const Particle&   a_Pdata )
{

}
   
void ParticleUtils::borisPusher( std::array<Real,3>&  a_vp,
                           const std::array<Real,3>&  a_vpold,
                           const std::array<Real,3>&  a_Ep,
                           const std::array<Real,3>&  a_Bp,
                           const Real&                a_fnorm_const,
                           const Real&                a_cnormDt ) 
{
   CH_TIME("ParticleUtils::borisPusher()");

   const Real zeroValue = 0.0;
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

#include "NamespaceFooter.H"


