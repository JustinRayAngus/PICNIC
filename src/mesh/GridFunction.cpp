#include "GridFunction.H"
#include "DomainGrid.H"

#include "NamespaceHeader.H"

GridFunction::GridFunction(const int& a_verbosity )

   : m_verbosity(a_verbosity)

{
}


void GridFunction::assign( LevelData<FArrayBox>&  a_data,
                     const DomainGrid&            a_mesh,
                     const Real&                  a_time ) const
{
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const LevelData<FArrayBox>& real_coords =  a_mesh.getXcc();
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      setPointwise( a_data[dit], a_mesh, real_coords[dit] );
   }
   a_data.exchange();
}
      
void GridFunction::assign( FArrayBox&   a_data,
                     const DomainGrid&  a_mesh,
                     const FArrayBox&   a_real_coords,
                     const Real&        a_time ) const
{
   setPointwise( a_data, a_mesh, a_real_coords );
}

      /// Print object parameters.


#include "NamespaceFooter.H"
