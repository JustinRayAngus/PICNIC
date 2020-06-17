#include "ParticleUtils.H"
//#include "ParticleUtilsF_F.H"

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


