#ifndef _SPECIAL_OPERATOR_H_
#define _SPECIAL_OPERATOR_H_

#include "REAL.H"
#include "PicSpecies.H"
#include "LevelData.H"

#include "DomainGrid.H"

#include "NamespaceHeader.H"

/**
 * Generic special operator interface.
*/
class SpecialOperator
{
   public:

      // Constructor using a database.
      //
      SpecialOperator( const int& verbosity );
   
      // Destructor.
      //
      virtual ~SpecialOperator() {;}

      // apply the special operator to passed pic species
      //
      virtual void applyOp( PicSpecies&  a_picSpecies,
                      const DomainGrid&  a_mesh,
                      const Real&        a_time ) const = 0;

      // update the special operator
      //
      virtual void updateOp( const DomainGrid&  a_mesh,
                             const Real&        a_time ) = 0;

      // Print object parameters.
      //
      virtual void printParameters() const = 0;
   
   protected:
   
      bool m_verbosity;
  
};

#include "NamespaceFooter.H"

#endif
