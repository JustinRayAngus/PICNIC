#ifndef _SPECIAL_OPERATOR_H_
#define _SPECIAL_OPERATOR_H_

#include "REAL.H"
#include "PicChargedSpecies.H"
#include "LevelData.H"
#include "CodeUnits.H"

#include "DomainGrid.H"

#include "NamespaceHeader.H"

/**
 * Generic special operator interface.
*/
class SpecialOperator
{
   public:

      SpecialOperator( const DomainGrid&  a_mesh,
                       const CodeUnits&   a_units,
                       const int&         a_verbosity );
   
      virtual ~SpecialOperator() {;}

      virtual void applyOp( PicChargedSpeciesPtrVect&  a_pic_species_ptr_vect,
                      const Real&               a_dt ) const = 0;

      virtual void updateOp( const Real&  a_dt ) = 0;
      
      virtual Real specialOpsDt() const = 0;

      virtual void printParameters() const = 0;
   
   protected:
   
      const DomainGrid& m_mesh;
      const CodeUnits& m_units;
      bool m_verbosity;
  
};

#include "NamespaceFooter.H"

#endif
