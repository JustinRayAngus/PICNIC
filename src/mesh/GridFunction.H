#ifndef _GRID_FUNCTION_H_
#define _GRID_FUNCTION_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "EdgeDataBox.H"
#include "LevelData.H"

#include "DomainGrid.H"

#include "NamespaceHeader.H"

/**
 * Generic grid function interface.
*/
class GridFunction
{
   public:

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      GridFunction( const int& verbosity );
   
      /// Destructor.
      /**
       */
      virtual ~GridFunction() {;}

      /// Set the grid function.
      /**
       * Sets the grid function on the LevelData for the given geometry
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry.  Cell
       * average values are set unless the optional parameter is false, in
       * which case, pointwise values are set.
       *
       * @param[out] data LevelData to set.
       * @param[in] geometry Corresponding level geometry.
       * @param[in] time Time value at which to evaluate the function.
       * @param[in] as_cell_name Boolean switch for cell average values or
       *                         point values.
       */
      void assign( LevelData<FArrayBox>&  a_data,
             const DomainGrid&            a_mesh,
             const Real&                  a_time ) const;


      /// Set the grid function.
      /**
       * Sets the grid function on the box for the given geometry
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry.  Cell
       * average values are set unless the optional parameter is false, in
       * which case, pointwise values are set.
       *
       * @param[out] data LevelData to set.
       * @param[in] geometry Corresponding level geometry.
       * @param[in] box Box from which to get coordinate system.
       * @param[in] time Time value at which to evaluate the function.
       * @param[in] as_cell_name Boolean switch for cell average values or
       *                         point values.
       */
      void assign( FArrayBox&   a_dfn,
             const DomainGrid&  a_mesh,
             const FArrayBox&   a_real_coords,
             const Real&        a_time ) const;

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const = 0;
   
   protected:
   
      /// Set pointwise values at cell centers.
      /**
       */
      virtual void setPointwise( FArrayBox&   a_dfn,
                           const DomainGrid&  a_mesh,
                           const FArrayBox&   a_real_coords ) const = 0;
   
      bool m_verbosity;
  
};

#include "NamespaceFooter.H"

#endif
