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
      setPointwise( a_data[dit], real_coords[dit] );
   }
   //a_data.exchange();
}

void GridFunction::assign( LevelData<FluxBox>&  a_data,
                     const int                  a_dir,
                     const DomainGrid&          a_mesh,
                     const Real&                a_time ) const
{
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const LevelData<FluxBox>& real_coords =  a_mesh.getXfc();
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      FArrayBox& this_data_dir = a_data[dit][a_dir];
      setPointwise( this_data_dir, real_coords[dit][a_dir] );
   }
   //a_data.exchange();

}

void GridFunction::assign( LevelData<EdgeDataBox>&  a_data,
                     const int                      a_dir,
                     const DomainGrid&              a_mesh,
                     const Real&                    a_time ) const
{
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const LevelData<EdgeDataBox>& real_coords =  a_mesh.getXec();
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      FArrayBox& this_data_dir = a_data[dit][a_dir];
      setPointwise( this_data_dir, real_coords[dit][a_dir] );
   }
   //a_data.exchange();

}

void GridFunction::assign( LevelData<NodeFArrayBox>&  a_data,
                     const DomainGrid&                a_mesh,
                     const Real&                      a_time ) const
{
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const LevelData<NodeFArrayBox>& real_coords =  a_mesh.getXnc();
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      FArrayBox& this_data = a_data[dit].getFab();
      const FArrayBox& this_real_coords = real_coords[dit].getFab();
      setPointwise( this_data, this_real_coords );
   }
   //a_data.exchange();

}
      

#include "NamespaceFooter.H"
