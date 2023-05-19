#ifndef _GRID_FUNCTION_H_
#define _GRID_FUNCTION_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"

#include "DomainGrid.H"

#include "NamespaceHeader.H"

/**
 * Generic grid function interface.
*/
class GridFunction
{
   public:

      GridFunction( const int& verbosity );
   
      virtual ~GridFunction() {;}

      // Set the grid function on cell center
      void assign( LevelData<FArrayBox>&  a_data,
             const DomainGrid&            a_mesh,
             const Real&                  a_time ) const;
      
      // Set the grid function on face center
      void assign( LevelData<FluxBox>&  a_data,
             const int                  a_dir,
             const DomainGrid&          a_mesh,
             const Real&                a_time ) const;
      
      // Set the grid function on edge center
      void assign( LevelData<EdgeDataBox>&  a_data,
             const int                      a_dir,
             const DomainGrid&              a_mesh,
             const Real&                    a_time ) const;
      
      // Set the grid function on node center
      void assign( LevelData<NodeFArrayBox>&  a_data,
             const DomainGrid&                a_mesh,
             const Real&                      a_time ) const;
      
      // return function value at location a_X
      virtual Real getValue( const RealVect&  a_X ) const = 0;
      
      // return dir derivative of function value at location a_X
      virtual Real getDerivative( const RealVect&  a_X,
                                  const int        a_dir ) const = 0;

      // Print object parameters.
      virtual void printParameters() const = 0;
   
   protected:
   
      virtual void setPointwise( FArrayBox&   a_dfn,
                           const FArrayBox&   a_real_coords ) const = 0;
      
      bool m_verbosity;
  
};

#include "NamespaceFooter.H"

#endif
