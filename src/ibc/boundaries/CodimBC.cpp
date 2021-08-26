#include "CodimBC.H"
#include "CodimBCF_F.H"

#include "BoundaryLookupTable.H"

#include "NamespaceHeader.H"

inline
Box CodimBC::getGhostBox( const int a_codim,
                          const int& a_current_index,
                          const Box& a_box,
                          const Box& a_domain_box,
                          const Vector<int>& a_dir,
                          const Vector<Side::LoHiSide>& a_side,
                          const IntVect& a_ghosts )
{
   const int dir( a_dir[a_current_index] );
   const Side::LoHiSide side( a_side[a_current_index] );
   const int num_ghosts( a_ghosts[dir] );
   const Box gBox( adjCellBox( a_box, dir, side, num_ghosts ) );
   if (a_current_index+1<a_codim) {
      return getGhostBox( a_codim,
                          a_current_index+1,
                          gBox,
                          a_domain_box,
                          a_dir,
                          a_side,
                          a_ghosts );
   }
   Box gDomainBox( a_domain_box );
   gDomainBox.grow( a_ghosts );
   return (gBox & gDomainBox);
}


inline
void CodimBC::fillCodimGhostCells(
   FArrayBox&                     a_this_soln,
   const Box&                     a_boundary_box,
   const Vector<int>&             a_dirs,
   const Vector<Side::LoHiSide>&  a_sides,
   const int&                     a_codim )
{
   Vector<int> sides( a_sides.size() );
   for (int i(0); i<sides.size(); i++) {
      sides[i] = static_cast<int>(a_sides[i]);
   }
   FORT_FILL_CODIM_GHOST_CELLS( CHF_FRA(a_this_soln),
                                CHF_BOX(a_boundary_box),
                                CHF_CONST_I1D(&(a_dirs[0]),a_dirs.size()),
                                CHF_CONST_I1D(&(sides[0]),sides.size()),
                                CHF_CONST_INT(a_codim) );
}

void CodimBC::setCodimCornerValues( LevelData<FArrayBox>&  a_u,
                              const DomainGrid&            a_mesh )
{
   CH_TIME("CodimBC::setCodimCornerValues()");
   
   const BoundaryLookupTable& boundaries( BoundaryLookupTable::getLookupTable() );
               
   const ProblemDomain& prob_domain = a_mesh.getDomain();
   const Box& domain_box = prob_domain.domainBox();

   const DisjointBoxLayout& grids( a_u.disjointBoxLayout() );
   for (DataIterator dit( a_u.dataIterator() ); dit.ok(); ++dit) {
      for (int codim(2); codim<=SpaceDim; codim++) {
         for (int m(0); m<boundaries.numberOfBoundaryCases(codim); m++) {

            const Vector<int>& dirs( boundaries.getDirections( m, codim ) );
            const Vector<Side::LoHiSide>& sides( boundaries.getSides( m, codim ) );
            const Box& box( grids[dit] );

            if(SpaceDim==2) { // hack to not do this if either dim is periodic in 2D (what about 3D?)
               if(prob_domain.isPeriodic(0) || prob_domain.isPeriodic(1)) break; 
            }

            if (isPhysicalCorner( domain_box, box, dirs, sides, codim )) {
               int CURRENT_INDEX(0);

               Box boundary_box( CodimBC::getGhostBox( codim,
                                                       CURRENT_INDEX,
                                                       box,
                                                       domain_box,
                                                       dirs,
                                                       sides,
                                                       a_u.ghostVect() ) );
               if (!boundary_box.isEmpty()) {
                  fillCodimGhostCells( a_u[dit], boundary_box, dirs, sides, codim );
               }
            }

         }
      }
   }

}

#include "NamespaceFooter.H"


