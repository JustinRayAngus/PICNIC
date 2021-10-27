#ifndef  _CODEUNITS_H_
#define  _CODEUNITS_H_

#include <iostream>

#include "REAL.H"
#include "ParmParse.H"
#include "Box.H"

#include "NamespaceHeader.H"

/**
 * code units class.
 *
 * Responsible for converting between dimensionless and dimensional
 * representations.
 *
*/
class CodeUnits
{
   public:

      // Enumerated characteristic scales.
      enum CHAR_SCALE {
         TEMPERATURE,      // [eV]
         NUMBER_DENSITY,   // [num/m^3]
         LENGTH,           // [m]
         AREA,             // [m^2]
         VOLUME,           // [m^3]
         TIME,             // [s]
         ELECTRIC_FIELD,   // [kV/cm]
         MAGNETIC_FIELD,   // [kV/cm/CVAC]
         COUNT};
      

      // Constructor.
      CodeUnits( ParmParse& parm_parse );

      // Default Destructor.
      ~CodeUnits() {;}

      void printParameters( const int a_procID ) const;

      Real getScale( const CHAR_SCALE& type ) const { checkType( type ); return m_scale[type]; }
      
      Real CvacNorm() const { return m_CVAC_NORM; }
      Real NumDenNorm() const { return m_NUM_DEN_NORM; }
     
   private:

      Real m_scale[COUNT];
      
      Real m_CVAC_NORM;
      Real m_NUM_DEN_NORM;

      // Check that the input type is valid.
      void checkType( const CHAR_SCALE& type ) const
      { CH_assert( type>=0 && type<COUNT ) };
      
};

#include "NamespaceFooter.H"

#endif
