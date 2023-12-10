#ifndef _SPACEUTILS_H_
#define _SPACEUTILS_H_

//#include "Dimensions.H"
//#include "BoundaryBoxLayout.H.multidim"

#include <math.h>
#include <vector>

#include "REAL.H"
#include "Box.H"
#include "LevelData.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"
#include "DomainGrid.H"

#include "NamespaceHeader.H"

namespace SpaceUtils {

  void applyBinomialFilter( LevelData<EdgeDataBox>&  a_var );
  
  void applyBinomialFilter( LevelData<NodeFArrayBox>&  a_var );
  
  void applyBinomialFilter( FArrayBox&  a_Q,
                      const Box&        a_grid_box );
      
  void computeLaxSplitting( LevelData<FArrayBox>&,
                            LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>& );
  
  void computeLaxSplitting( LevelData<FArrayBox>&,
                            LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const int );
  
  void computeLaxSplitting( FArrayBox&,
                            FArrayBox&,
                      const FArrayBox&,
                      const FArrayBox&,
                      const FArrayBox& );
  
  void interpToFaces( LevelData<FluxBox>&,
                const LevelData<FArrayBox>&,
                const LevelData<FArrayBox>&,
                const LevelData<FluxBox>&,
                const std::string& );
  
  void interpToFaces( FluxBox&,
                const FArrayBox&,
                const FArrayBox&,
                const FluxBox&,
                const Box&,
                const std::string& );

  void upWindToFaces( LevelData<FluxBox>&,
                const LevelData<FArrayBox>&,
                const LevelData<FluxBox>&,
                const std::string& );
  
  void upWindToFaces( FluxBox&,
                const FArrayBox&,
                const FluxBox&,
                const Box&,
                const std::string& );
  
  
  void interpCellToEdges( LevelData<EdgeDataBox>&,
                    const LevelData<FArrayBox>&,
                    const LevelData<EdgeDataBox>&,
                    const std::string& );
  
  void interpEdgesToCell( LevelData<FArrayBox>&,
                    const LevelData<EdgeDataBox>&,
                    const std::string& );
  
  void interpEdgesToEdges( LevelData<EdgeDataBox>&,
                     const LevelData<EdgeDataBox>&,
                     const std::string& );

  void interpNodesToEdges( LevelData<EdgeDataBox>&,
                     const LevelData<NodeFArrayBox>&,
                     const std::string& );

  void interpNodesToCells( LevelData<FArrayBox>&,
                     const LevelData<NodeFArrayBox>&,
                     const std::string& );

  void interpCellsToNodes( LevelData<NodeFArrayBox>&,
                     const LevelData<FArrayBox>&,
                     const std::string& );

  void PerpEdgeGradientAtCells( FArrayBox&,
                          const EdgeDataBox&,
                          const RealVect&,
                          const Box&,
                          const std::string& );

  void interpFacesToCell( LevelData<FArrayBox>&,
                    const LevelData<FluxBox>&,
                    const std::string& );
  
  void interpToFacesWENO( LevelData<FluxBox>&,
                    const LevelData<FArrayBox>&,
                    const LevelData<FluxBox>&,
                    const LevelData<FArrayBox>&,
                    const std::string& );
  
  void interpToFacesWENO( FluxBox&,
                    const FArrayBox&,
                    const FluxBox&,
                    const FArrayBox&,
                    const Box&,
                    const std::string& );

  void faceInterpolate( const int,
                        const Box&,
                        const int,
                        const FArrayBox&,
                        FArrayBox& );

  void faceInterpolate( const int,
                        const Box&,
                        const Box&,
                        const int,
                        const FArrayBox&,
                        FArrayBox& );

  void interpolateStag( FArrayBox&  a_dst, 
                  const Box&        a_dst_box,
                  const int         a_dst_comp,
                  const FArrayBox&  a_src,
                  const int         a_src_comp,
                  const int         a_dir );
  
  void simpleStagGradComp( FArrayBox&  a_dst,
                     const Box&        a_grid_box,
                     const int         a_dst_comp,
                     const FArrayBox&  a_src,
                     const int         a_src_comp,
                     const Real        a_dX,
                     const int         a_dir_dX,
                     const int         a_additive ); 

  void cellCenteredGradientComponent( const Box&,
                                      const int,
                                      const FArrayBox&,
                                      const RealVect&,
                                      const int,
                                      FArrayBox& );

  void faceCenteredGradientComponent( const Box&,
                                      const int,
                                      const FArrayBox&,
                                      const RealVect&,
                                      const int,
                                      FArrayBox& );

  void extrapBoundaryGhostsForCC(  FArrayBox&, 
                                   const Box&, 
                                   const Box&, 
                                   const int, 
                                   const int, 
                                   const int);
  
  void extrapBoundaryGhostsForFC(  FArrayBox&, 
                                   const Box&, 
                                   const Box&, 
                                   const int, 
                                   const int, 
                                   const int);

  void secondOrderTransExtrapAtDomainBdry(  FArrayBox&,
                                            const int,
                                            const Box&,
                                            const ProblemDomain&,
                                            const int );

  void fillGhostCellsSimple(  FArrayBox&, const Box&, const int);
  
  void copyAndFillGhostCellsSimple( LevelData<FArrayBox>&, 
                                    const LevelData<FArrayBox>& );
  
  void inspectFArrayBox( const FArrayBox&  a_F0,
                         const Box&        a_box,
                         const int         a_comp );

  void inspectFArrayBox( const LevelData<FArrayBox>&,
                         const int,
                         const int                   );
  
  void inspectFluxBox( const LevelData<FluxBox>&,
                       const int                 );
  
  void inspectEdgeDataBox( const LevelData<EdgeDataBox>&,
                           const int                     );
  
  void inspectNodeFArrayBox( const LevelData<NodeFArrayBox>&,
                             const int                     );

  void exchangeFluxBox( LevelData<FluxBox>& );
  
  void exchangeEdgeDataBox( LevelData<EdgeDataBox>& );
  
  void exchangeNodeFArrayBox( LevelData<NodeFArrayBox>& );

  void checkForNAN( const LevelData<FArrayBox>&  a_var,
   	            const std::string&           a_string );
  
  void checkForNAN( const LevelData<FluxBox>&  a_var,
   	            const std::string&         a_string );

  void checkForNAN( const LevelData<EdgeDataBox>&  a_var,
   	            const std::string&             a_string );
  
  void checkForNAN( const LevelData<NodeFArrayBox>&  a_var,
   	            const std::string&               a_string );

  void checkForNAN( const FArrayBox&    a_var,
   	            const std::string&  a_string );

  void copyEdgeDataBox( EdgeDataBox&,
                  const EdgeDataBox& );

  void copy( FArrayBox&,
	     const FArrayBox& );

  void copy( FArrayBox&,
       const FArrayBox&,
       const Box&        );

  void copyNodeToCell( FArrayBox&,
                 const FArrayBox&,
                 const Box&,
                 const int );
  
  void setVal( FArrayBox&  a_dst,
	 const Real        a_val,
	 const int         a_comp );


  void localVectorNorm( FArrayBox&,
                  const FArrayBox& );

  inline void zero(  LevelData<FArrayBox>& a_u ) 
  {
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      a_u[dit].setVal(0.0);
    }
  }

  inline void zero(  LevelData<FArrayBox>& a_u, const int a_comp ) 
  {
    CH_assert((a_comp >= 0) && (a_comp < a_u.nComp()));
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      a_u[dit].setVal(0.0, a_comp);
    }
  }

  inline void setVal(  LevelData<FArrayBox>& a_u, const int a_comp, const Real a_val ) 
  {
    CH_assert((a_comp >= 0) && (a_comp < a_u.nComp()));
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      a_u[dit].setVal(a_val, a_comp);
    }
  }

  inline void zero(  LevelData<FluxBox>& a_u ) 
  {
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      for (auto dir(0); dir<SpaceDim; dir++) {
        a_u[dit][dir].setVal(0.0);
      }
    }
  }

  inline void zero(  LevelData<FluxBox>& a_u, const int a_comp ) 
  {
    CH_assert((a_comp >= 0) && (a_comp < a_u.nComp()));
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      for (auto dir(0); dir<SpaceDim; dir++) {
        a_u[dit][dir].setVal(0.0, a_comp);
      }
    }
  }

  inline void setVal(  LevelData<FluxBox>& a_u, const int a_comp, const Real a_val ) 
  {
    CH_assert((a_comp >= 0) && (a_comp < a_u.nComp()));
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      for (auto dir(0); dir<SpaceDim; dir++) {
        a_u[dit][dir].setVal(a_val, a_comp);
      }
    }
  }

  inline void zero(  LevelData<EdgeDataBox>& a_u ) 
  {
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      for (auto dir(0); dir<SpaceDim; dir++) {
        a_u[dit][dir].setVal(0.0);
      }
    }
  }

  inline void zero(  LevelData<EdgeDataBox>& a_u, const int a_comp ) 
  {
    CH_assert((a_comp >= 0) && (a_comp < a_u.nComp()));
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      for (auto dir(0); dir<SpaceDim; dir++) {
        a_u[dit][dir].setVal(0.0, a_comp);
      }
    }
  }

  inline void setVal(  LevelData<EdgeDataBox>& a_u, const int a_comp, const Real a_val ) 
  {
    CH_assert((a_comp >= 0) && (a_comp < a_u.nComp()));
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      for (auto dir(0); dir<SpaceDim; dir++) {
        a_u[dit][dir].setVal(a_val, a_comp);
      }
    }
  }

  inline void zero(  LevelData<NodeFArrayBox>& a_u ) 
  {
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      a_u[dit].setVal(0.0);
    }
  }

  inline void zero(  LevelData<NodeFArrayBox>& a_u, const int a_comp ) 
  {
    CH_assert((a_comp >= 0) && (a_comp < a_u.nComp()));
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      a_u[dit].setVal(0.0, surroundingNodes(a_u[dit].box()), a_comp, 1);
    }
  }

  inline void setVal(  LevelData<NodeFArrayBox>& a_u, const int a_comp, const Real a_val ) 
  {
    CH_assert((a_comp >= 0) && (a_comp < a_u.nComp()));
    for (auto dit(a_u.dataIterator()); dit.ok(); ++dit) {
      a_u[dit].setVal(a_val, surroundingNodes(a_u[dit].box()), a_comp, 1);
    }
  }

  inline
  void copyLevelData( const LevelData<FArrayBox>& a_src,
                      LevelData<FArrayBox>&       a_dst )
  {
    for (DataIterator dit(a_dst.dataIterator()); dit.ok(); ++dit) {
      a_dst[dit].copy( a_src[dit] );
    }
  }

  inline
  void copyLevelData( const LevelData<FluxBox>& a_src,
                      LevelData<FluxBox>&       a_dst )
  {
    for (DataIterator dit(a_dst.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        a_dst[dit][dir].copy( a_src[dit][dir] );
      }
    }
  }

  inline
  void copyLevelData( const LevelData<EdgeDataBox>& a_src,
                      LevelData<EdgeDataBox>&       a_dst )
  {
    for (DataIterator dit(a_dst.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        a_dst[dit][dir].copy( a_src[dit][dir] );
      }
    }
  }

  inline
  void copyLevelData( const LevelData<NodeFArrayBox>& a_src,
                      LevelData<NodeFArrayBox>&       a_dst )
  {
    for (DataIterator dit(a_dst.dataIterator()); dit.ok(); ++dit) {
      a_dst[dit].copy( a_src[dit] );
    }
  }

  inline
  void addToLevelData(  const LevelData<FArrayBox>& a_src,
                        LevelData<FArrayBox>&       a_dst,
                        const Real                  a_a = 1 )
  {
    for (DataIterator dit(a_dst.dataIterator()); dit.ok(); ++dit) {
      a_dst[dit].plus( a_src[dit], a_a );
    }
  }

  inline
  void addToLevelData(  const LevelData<FluxBox>& a_src,
                        LevelData<FluxBox>&       a_dst,
                        const Real                a_a = 1 )
  {
    for (DataIterator dit(a_dst.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        a_dst[dit][dir].plus( a_src[dit][dir], a_a );
      }
    }
  }

  inline
  void addToLevelData(  const LevelData<EdgeDataBox>& a_src,
                        LevelData<EdgeDataBox>&       a_dst,
                        const Real                    a_a = 1 )
  {
    for (DataIterator dit(a_dst.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        a_dst[dit][dir].plus( a_src[dit][dir], a_a );
      }
    }
  }

  inline
  void addToLevelData(  const LevelData<NodeFArrayBox>& a_src,
                        LevelData<NodeFArrayBox>&       a_dst,
                        const Real                      a_a = 1 )
  {
    for (DataIterator dit(a_dst.dataIterator()); dit.ok(); ++dit) {
      a_dst[dit].plus( a_src[dit], a_a );
    }
  }

  inline
  Real sumOfSquares(  const FArrayBox&  a_fab )
  {
    Real retval = 0;
    for (int comp = 0; comp < a_fab.nComp(); comp++) {    
      retval += a_fab.sumPow( a_fab.box(), 2, comp );
    }
    return retval;
  }

  inline
  Real sumOfSquares(  const NodeFArrayBox&  a_fab )
  {
    Real retval = 0;
    for (int comp = 0; comp < a_fab.nComp(); comp++) {    
      retval += a_fab.sumPow( surroundingNodes(a_fab.box()), 2, comp );
    }
    return retval;
  }

  inline
  void sumOfSquares(  const FArrayBox&    a_fab,
                      std::vector<Real>&  a_vec )
  {
    a_vec.resize( a_fab.nComp() );
    for (int comp = 0; comp < a_fab.nComp(); comp++) {    
      a_vec[comp] = a_fab.sumPow( a_fab.box(), 2, comp );
    }
    return;
  }

  inline
  void sumOfSquares(  const NodeFArrayBox&  a_fab,
                      std::vector<Real>&    a_vec )
  {
    a_vec.resize( a_fab.nComp() );
    for (int comp = 0; comp < a_fab.nComp(); comp++) {    
      a_vec[comp] = a_fab.sumPow( surroundingNodes(a_fab.box()), 2, comp );
    }
    return;
  }

  inline
  Real sumOfSquares( const LevelData<FArrayBox>& a_var )
  {
    Real sum_of_squares = 0;
    for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
      sum_of_squares += sumOfSquares( a_var[dit] );
    }
    Real sum_of_squares_glbl = 0;
#ifdef CH_MPI
    MPI_Allreduce(  &sum_of_squares,
                    &sum_of_squares_glbl,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);
#else
    sum_of_squares_glbl = sum_of_squares;
#endif
    return sum_of_squares_glbl;
  }

  inline
  void sumOfSquares(  const LevelData<FArrayBox>& a_var,
                      std::vector<Real>&          a_vec )
  {
    std::vector<Real> vec_lcl(a_var.nComp(), 0);
    for (DataIterator dit( a_var.dataIterator() ); dit.ok(); ++dit) {
      std::vector<Real> tmp( vec_lcl.size() );
      sumOfSquares( a_var[dit], tmp );
      for (int n = 0; n < a_var.nComp(); n++) {
        a_vec[n] += tmp[n];
      }
    }
    a_vec.resize( a_var.nComp() );
#ifdef CH_MPI
    for (int n = 0; n < a_var.nComp(); n++) {
      a_vec[n] = 0;
    }
    MPI_Allreduce(  vec_lcl.data(),
                    a_vec.data(),
                    a_var.nComp(),
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);
#else
    for (int n = 0; n < a_var.nComp(); n++) {
      a_vec[n] = vec_lcl[n];
    }
#endif
    return;
  }

  inline
  Real sumOfSquares( const LevelData<FluxBox>& a_var )
  {
    Real sum_of_squares = 0;
    for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        sum_of_squares += sumOfSquares( a_var[dit][dir] );
      }
    }
    Real sum_of_squares_glbl = 0;
#ifdef CH_MPI
    MPI_Allreduce(  &sum_of_squares,
                    &sum_of_squares_glbl,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);
#else
    sum_of_squares_glbl = sum_of_squares;
#endif
    return sum_of_squares_glbl;
  }

  inline
  void sumOfSquares(  const LevelData<FluxBox>& a_var,
                      std::vector<Real>&        a_vec )
  {
    std::vector<Real> vec_lcl(a_var.nComp(), 0);
    for (DataIterator dit( a_var.dataIterator() ); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        std::vector<Real> tmp( vec_lcl.size() );
        sumOfSquares( a_var[dit][dir], tmp );
        for (int n = 0; n < a_var.nComp(); n++) {
          a_vec[n] += tmp[n];
        }
      }
    }
    a_vec.resize( a_var.nComp() );
#ifdef CH_MPI
    for (int n = 0; n < a_var.nComp(); n++) {
      a_vec[n] = 0;
    }
    MPI_Allreduce(  vec_lcl.data(),
                    a_vec.data(),
                    a_var.nComp(),
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);
#else
    for (int n = 0; n < a_var.nComp(); n++) {
      a_vec[n] = vec_lcl[n];
    }
#endif
    return;
  }

  inline
  Real sumOfSquares( const LevelData<EdgeDataBox>& a_var )
  {
    Real sum_of_squares = 0;
    for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        sum_of_squares += sumOfSquares( a_var[dit][dir] );
      }
    }
    Real sum_of_squares_glbl = 0;
#ifdef CH_MPI
    MPI_Allreduce(  &sum_of_squares,
                    &sum_of_squares_glbl,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);
#else
    sum_of_squares_glbl = sum_of_squares;
#endif
    return sum_of_squares_glbl;
  }

  inline
  void sumOfSquares(  const LevelData<EdgeDataBox>& a_var,
                      std::vector<Real>&            a_vec )
  {
    std::vector<Real> vec_lcl(a_var.nComp(), 0);
    for (DataIterator dit( a_var.dataIterator() ); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        std::vector<Real> tmp( vec_lcl.size() );
        sumOfSquares( a_var[dit][dir], tmp );
        for (int n = 0; n < a_var.nComp(); n++) {
          a_vec[n] += tmp[n];
        }
      }
    }
    a_vec.resize( a_var.nComp() );
#ifdef CH_MPI
    for (int n = 0; n < a_var.nComp(); n++) {
      a_vec[n] = 0;
    }
    MPI_Allreduce(  vec_lcl.data(),
                    a_vec.data(),
                    a_var.nComp(),
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);
#else
    for (int n = 0; n < a_var.nComp(); n++) {
      a_vec[n] = vec_lcl[n];
    }
#endif
    return;
  }

  inline
  Real sumOfSquares( const LevelData<NodeFArrayBox>& a_var )
  {
    Real sum_of_squares = 0;
    for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
      sum_of_squares += sumOfSquares( a_var[dit] );
    }
    Real sum_of_squares_glbl = 0;
#ifdef CH_MPI
    MPI_Allreduce(  &sum_of_squares,
                    &sum_of_squares_glbl,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);
#else
    sum_of_squares_glbl = sum_of_squares;
#endif
    return sum_of_squares_glbl;
  }

  inline
  void sumOfSquares(  const LevelData<NodeFArrayBox>& a_var,
                      std::vector<Real>&              a_vec )
  {
    std::vector<Real> vec_lcl(a_var.nComp(), 0);
    for (DataIterator dit( a_var.dataIterator() ); dit.ok(); ++dit) {
      std::vector<Real> tmp( vec_lcl.size() );
      sumOfSquares( a_var[dit], tmp );
      for (int n = 0; n < a_var.nComp(); n++) {
        a_vec[n] += tmp[n];
      }
    }
    a_vec.resize( a_var.nComp() );
#ifdef CH_MPI
    for (int n = 0; n < a_var.nComp(); n++) {
      a_vec[n] = 0;
    }
    MPI_Allreduce(  vec_lcl.data(),
                    a_vec.data(),
                    a_var.nComp(),
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);
#else
    for (int n = 0; n < a_var.nComp(); n++) {
      a_vec[n] = vec_lcl[n];
    }
#endif
    return;
  }

  inline
  Real computeNorm( const FArrayBox& a_fab )
  {
    return sqrt(sumOfSquares(a_fab));
  }

  inline
  void computeNorm( const FArrayBox&    a_fab,
                    std::vector<Real>&  a_vec )
  {
    sumOfSquares( a_fab, a_vec );
    for (int comp = 0; comp < a_fab.nComp(); comp++) {    
      a_vec[comp] = sqrt( a_vec[comp] );
    }
    return;
  }

  inline
  Real computeNorm( const LevelData<FArrayBox>& a_var )
  {
    return sqrt( sumOfSquares( a_var ) );
  }

  inline
  void computeNorm( const LevelData<FArrayBox>& a_var,
                    std::vector<Real>&          a_vec )
  {
    sumOfSquares( a_var, a_vec );
    for (int comp = 0; comp < a_var.nComp(); comp++) {    
      a_vec[comp] = sqrt( a_vec[comp] );
    }
    return;
  }

  inline
  Real computeNorm( const LevelData<FluxBox>& a_var )
  {
    return sqrt( sumOfSquares( a_var ) );
  }

  inline
  void computeNorm( const LevelData<FluxBox>& a_var,
                    std::vector<Real>&        a_vec )
  {
    sumOfSquares( a_var, a_vec );
    for (int comp = 0; comp < a_var.nComp(); comp++) {    
      a_vec[comp] = sqrt( a_vec[comp] );
    }
    return;
  }

  inline
  Real computeNorm( const LevelData<EdgeDataBox>& a_var )
  {
    return sqrt( sumOfSquares( a_var ) );
  }

  inline
  void computeNorm( const LevelData<EdgeDataBox>& a_var,
                    std::vector<Real>&            a_vec )
  {
    sumOfSquares( a_var, a_vec );
    for (int comp = 0; comp < a_var.nComp(); comp++) {    
      a_vec[comp] = sqrt( a_vec[comp] );
    }
    return;
  }

  inline
  Real computeNorm( const LevelData<NodeFArrayBox>& a_var )
  {
    return sqrt( sumOfSquares( a_var ) );
  }

  inline
  void computeNorm( const LevelData<NodeFArrayBox>& a_var,
                    std::vector<Real>&              a_vec )
  {
    sumOfSquares( a_var, a_vec );
    for (int comp = 0; comp < a_var.nComp(); comp++) {    
      a_vec[comp] = sqrt( a_vec[comp] );
    }
    return;
  }

  inline
  int copyFromLevelData( Real* a_dst, const LevelData<FArrayBox>& a_src )
  {
     int offset(0);
     auto dbl( a_src.disjointBoxLayout() );
     for (auto dit( a_src.dataIterator() ); dit.ok(); ++dit) {
        FArrayBox tmp( dbl[dit], a_src.nComp(), ( a_dst + offset ) );
        tmp.copy( a_src[dit] );
        offset += dbl[dit].numPts() * a_src.nComp();
     }
     return offset;
  }
  
  inline
  int copyFromLevelData( Real* a_dst, const LevelData<FluxBox>& a_src )
  {
     int offset(0);
     auto dbl( a_src.disjointBoxLayout() );
     auto phys_domain( dbl.physDomain() );

     for (auto dit( a_src.dataIterator() ); dit.ok(); ++dit) {
        for (int dir=0; dir<SpaceDim; ++dir) {
           
           Box box = surroundingNodes(dbl[dit],dir);

           int idir_bdry = phys_domain.domainBox().bigEnd(dir);
           //bool is_periodic( phys_domain.isPeriodic(dir) );
           //if ((is_periodic) || (box.bigEnd(dir) < idir_bdry)) box.growHi(dir, -1);
           if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);

           for (int n=0; n<a_src.nComp(); n++) {
              FArrayBox tmp( box, 1, ( a_dst + offset ) );
              tmp.copy( a_src[dit][dir], box, n, box, 0, 1);
              offset += box.numPts();
           }
        }
     }
     return offset;
  }
  
  inline
  int copyFromLevelData( Real* a_dst, const LevelData<EdgeDataBox>& a_src )
  {  
     int offset(0);
     auto dbl( a_src.disjointBoxLayout() );
     auto phys_domain( dbl.physDomain() );

     for (auto dit( a_src.dataIterator() ); dit.ok(); ++dit) {
        for (int dir=0; dir<SpaceDim; ++dir) {

           Box box = enclosedCells( surroundingNodes(dbl[dit]), dir );
           for (int adir=0; adir<SpaceDim; ++adir) {
              if (adir != dir) {
                 int idir_bdry = phys_domain.domainBox().bigEnd(adir);
                 //bool is_periodic( phys_domain.isPeriodic(adir) );
                 //if ((is_periodic) || (box.bigEnd(adir) < idir_bdry)) box.growHi(adir, -1);
                 if (box.bigEnd(adir) < idir_bdry) box.growHi(adir, -1);
              }
           }

           for (int n=0; n<a_src.nComp(); n++) {
              FArrayBox tmp( box, 1, ( a_dst + offset ) );
              tmp.copy( a_src[dit][dir], box, n, box, 0, 1 );
              offset += box.numPts();
           }
        }
     }
     return offset;
  }
  
  inline
  int copyFromLevelData( Real* a_dst, const LevelData<NodeFArrayBox>& a_src )
  {
     int offset(0);
     auto dbl( a_src.disjointBoxLayout() );
     auto phys_domain( dbl.physDomain() );

     for (auto dit( a_src.dataIterator() ); dit.ok(); ++dit) {

        Box box = surroundingNodes(dbl[dit]);
        for (int dir=0; dir<SpaceDim; ++dir) {
           int idir_bdry = phys_domain.domainBox().bigEnd(dir);
           //bool is_periodic( phys_domain.isPeriodic(dir) );
           //if ((is_periodic) || (box.bigEnd(dir) < idir_bdry)) box.growHi(dir, -1);
           if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
        }

        for (int n=0; n<a_src.nComp(); n++) {
           FArrayBox tmp( box, 1, ( a_dst + offset ) );
           tmp.copy( a_src[dit].getFab(), box, n, box, 0, 1 );
           offset += box.numPts();
        }
     }
     return offset;
  }
  
  inline
  int copyToLevelData( LevelData<FArrayBox>& a_dst, const Real* const a_src )
  {
     int offset(0);
     auto dbl( a_dst.disjointBoxLayout() );
     for (auto dit( a_dst.dataIterator() ); dit.ok(); ++dit) {
        FArrayBox tmp( dbl[dit], a_dst.nComp(), const_cast<Real*>( a_src + offset ) );
        a_dst[dit].copy( tmp );
        offset += dbl[dit].numPts() * a_dst.nComp();
     }
     a_dst.exchange();
     return offset;
  }
  
  inline
  int copyToLevelData( LevelData<FluxBox>& a_dst, const Real* const a_src )
  {
     int offset(0);
     auto dbl( a_dst.disjointBoxLayout() );
     auto phys_domain( dbl.physDomain() );

     for (auto dit( a_dst.dataIterator() ); dit.ok(); ++dit) {
        for (int dir=0; dir<SpaceDim; ++dir) {

           Box box = surroundingNodes(dbl[dit],dir);
           
           int idir_bdry = phys_domain.domainBox().bigEnd(dir);
           //bool is_periodic( phys_domain.isPeriodic(dir) );
           //if ((is_periodic) || (box.bigEnd(dir) < idir_bdry)) box.growHi(dir, -1);
           if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);

           for (int n=0; n<a_dst.nComp(); n++) {
              FArrayBox tmp( box, 1, const_cast<Real*>( a_src + offset ) );
              a_dst[dit][dir].copy( tmp, box, 0, box, n, 1 );
              offset += box.numPts();
           }
        }
     }
     exchangeFluxBox( a_dst );
     return offset;
  }
  
  inline
  int copyToLevelData( LevelData<EdgeDataBox>& a_dst, const Real* const a_src )
  {
     int offset(0);
     auto dbl( a_dst.disjointBoxLayout() );
     auto phys_domain( dbl.physDomain() );

     for (auto dit( a_dst.dataIterator() ); dit.ok(); ++dit) {
        for (int dir=0; dir<SpaceDim; ++dir) {

           Box box = enclosedCells( surroundingNodes(dbl[dit]), dir );
           for (int adir=0; adir<SpaceDim; ++adir) {
              if (adir != dir) {
                 int idir_bdry = phys_domain.domainBox().bigEnd(adir);
                 //bool is_periodic( phys_domain.isPeriodic(adir) );
                 //if ((is_periodic) || (box.bigEnd(adir) < idir_bdry)) box.growHi(adir, -1);
                 if (box.bigEnd(adir) < idir_bdry) box.growHi(adir, -1);
              }
           }

           for (int n=0; n<a_dst.nComp(); n++) {
              FArrayBox tmp( box, 1, const_cast<Real*>( a_src + offset ) );
              a_dst[dit][dir].copy( tmp, box, 0, box, n, 1 );
              offset += box.numPts();
           }
        }
     }

     exchangeEdgeDataBox(a_dst);
     return offset;
  }
   
  inline
  int copyToLevelData( LevelData<NodeFArrayBox>& a_dst, const Real* const a_src )
  {
     int offset(0);
     auto dbl( a_dst.disjointBoxLayout() );
     auto phys_domain( dbl.physDomain() );

     for (auto dit( a_dst.dataIterator() ); dit.ok(); ++dit) {

        Box box = surroundingNodes(dbl[dit]);
        for (int dir=0; dir<SpaceDim; ++dir) {
           int idir_bdry = phys_domain.domainBox().bigEnd(dir);
           //bool is_periodic( phys_domain.isPeriodic(dir) );
           //if ((is_periodic) || (box.bigEnd(dir) < idir_bdry)) box.growHi(dir, -1);
           if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
        }

        for (int n=0; n<a_dst.nComp(); n++) {
           FArrayBox tmp( box, 1, const_cast<Real*>( a_src + offset ) );
           a_dst[dit].getFab().copy( tmp, box, 0, box, n, 1 );
           offset += box.numPts();
        }
     }

     exchangeNodeFArrayBox( a_dst );
     return offset;
  }

  inline
  int nDOF( const LevelData<FArrayBox>& a_u )
  {
    int retval(0);
    auto dbl( a_u.disjointBoxLayout() );
    for (auto dit( dbl.dataIterator() ); dit.ok(); ++dit) {
      auto box( dbl[dit] );
      retval += box.numPts() * a_u.nComp();
    }
    return retval;
  }
  
  inline
  int nDOF( const LevelData<FluxBox>& a_u )
  {
    int retval(0);
    auto dbl( a_u.disjointBoxLayout() );
    auto phys_domain( dbl.physDomain() );

    for (auto dit( dbl.dataIterator() ); dit.ok(); ++dit) {
      for (int dir(0); dir < SpaceDim; dir++) {

        Box box = surroundingNodes(dbl[dit],dir);
        
        int idir_bdry = phys_domain.domainBox().bigEnd(dir);
        //bool is_periodic( phys_domain.isPeriodic(dir) );
        //if ((is_periodic) || (box.bigEnd(dir) < idir_bdry)) box.growHi(dir, -1);
        if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);

        retval += box.numPts() * a_u.nComp();
      }
    }
    return retval;
  }
  
  inline
  int nDOF( const LevelData<EdgeDataBox>& a_u )
  {
    int retval(0);
    auto dbl( a_u.disjointBoxLayout() );
    auto phys_domain( dbl.physDomain() );

    for (auto dit( dbl.dataIterator() ); dit.ok(); ++dit) {
      for (int dir(0); dir < SpaceDim; dir++) {

        Box box = enclosedCells( surroundingNodes(dbl[dit]), dir );
        for (int adir=0; adir<SpaceDim; ++adir) {
           if (adir != dir) {
              int idir_bdry = phys_domain.domainBox().bigEnd(adir);
              //bool is_periodic( phys_domain.isPeriodic(adir) );
              //if ((is_periodic) || (box.bigEnd(adir) < idir_bdry)) box.growHi(adir, -1);
              if (box.bigEnd(adir) < idir_bdry) box.growHi(adir, -1);
           }
        }

        retval += box.numPts() * a_u.nComp();
      }
    }

    return retval;
  }
  
  inline
  int nDOF( const LevelData<NodeFArrayBox>& a_u )
  {
    int retval(0);
    auto dbl( a_u.disjointBoxLayout() );
    auto phys_domain( dbl.physDomain() );

    for (auto dit( dbl.dataIterator() ); dit.ok(); ++dit) {

      auto box( surroundingNodes(dbl[dit]) );
      for (int dir=0; dir<SpaceDim; ++dir) {
        int idir_bdry = phys_domain.domainBox().bigEnd(dir);
        //bool is_periodic( phys_domain.isPeriodic(dir) );
        //if ((is_periodic) || (box.bigEnd(dir) < idir_bdry)) box.growHi(dir, -1);
        if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
      }
      retval += box.numPts() * a_u.nComp();

    }
    return retval;
  }
  
}

#include "NamespaceFooter.H"
#endif

