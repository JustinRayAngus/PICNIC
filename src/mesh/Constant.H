#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#include "GridFunction.H"

#include "NamespaceHeader.H"

/**
 * Generic constant function class.
 * Assigns LevelData of any dimension to a constant.
*/
class Constant : public GridFunction
{
   public:

      // Constructor using a value.
      Constant( const Real& value, const int& verbosity )
         : GridFunction(verbosity), m_value(value) {;}

      // Constructor using a database.
      Constant( ParmParse& pp, const int& verbosity )
         : GridFunction(verbosity), m_value(0.0)
      {
         parseParameters( pp );
      }

      virtual ~Constant() {;}

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Constant grid function parameters:" << std::endl;
            std::cout << "  value: "  << m_value  << std::endl << std::endl;
         }
      }

   private:

      // prohibit copying
      Constant( const Constant& );
      Constant& operator=( const Constant& );

      // Parse the input database for parameters.
      void parseParameters( ParmParse& pp )
      {
         pp.get( "value", m_value );
         if (m_verbosity) {
            printParameters();
         }
      }

      // Set pointwise values at cell centers.
      virtual void setPointwise( FArrayBox&   a_data,
                           const DomainGrid&  a_mesh,
                           const FArrayBox&   a_real_coords ) const
      {
         a_data.setVal(m_value);
      }
      
      Real m_value;
};

#include "NamespaceFooter.H"

#endif

