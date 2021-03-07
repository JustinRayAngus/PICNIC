#ifndef  _CODEUNITS_H_
#define  _CODEUNITS_H_

#include <iostream>

#include "REAL.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

/**
 * code units class.
 *
 * Responsible for converting between dimensionless and dimensional
 * representations.
 *
 * Required input keys and data types:
 *
 *    -\b temperature
 *      real value specifying the characteristic temperature in [\f$eV\f$]
 *
 *    -\b length
 *      real value specifying the characteristic length in [\f$m\f$]
 *
 *    -\b time
 *      real value specifying the characteristic time in [\f$s\f$]
 *
 *
 * The following represents a sample input entry:
 *
 * \verbatim
 *  units.temperature = 1.0
 *  units.length = 1.0e-2
 *  units.time = 1.0e-9 
 * \endverbatim
*/
class CodeUnits
{
   public:

      /// Enumerated characteristic scales.
      /**
       * The associate mks units are shown.
       */
      enum CHAR_SCALE {
         TEMPERATURE,           /**< [\f$eV\f$] */
         LENGTH,                /**< [\f$m\f$] */
         TIME,                  /**< [\f$s\f$] */
         COUNT};

      // Constructor.
      CodeUnits( ParmParse& parm_parse );

      // Default Destructor.
      ~CodeUnits() {;}

      void printParameters( const int a_procID ) const;

      Real getScale( const CHAR_SCALE& type ) const { checkType( type ); return m_scale[type]; }

   private:

      Real m_scale[COUNT];

      /// Check that the input type is valid.
      /**
       */
      void checkType( const CHAR_SCALE& type ) const
      { CH_assert( type>=0 && type<COUNT ) };
};

#include "NamespaceFooter.H"

#endif
