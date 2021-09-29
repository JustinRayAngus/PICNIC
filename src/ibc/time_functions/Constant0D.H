#ifndef _Constant0D_H_
#define _Constant0D_H_

#include <cmath>

#include "TimeFunction.H"

#include "NamespaceHeader.H"

class Constant0D : public TimeFunction
{
   public:

      Constant0D( const Real&  a_amplitude, 
                  const int&   verbosity )
         : TimeFunction(verbosity), 
           m_amplitude(a_amplitude) {;}

      Constant0D( ParmParse& pp, const int& verbosity )
         : TimeFunction(verbosity), 
           m_amplitude(0.0)
      {
         parseParameters( pp );
      }

      virtual ~Constant0D() {;}

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Constant0D function parameters:" << std::endl;
            std::cout << " amplitude = " << m_amplitude  << std::endl << std::endl;
         }
      }
      
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const
      {
         a_value = m_amplitude;
      }

   private:

      // prohibit copying
      Constant0D( const Constant0D& );
      Constant0D& operator=( const Constant0D& );

      void parseParameters( ParmParse& pp )
      {
         pp.get( "amplitude", m_amplitude );
         if (m_verbosity) {
            printParameters();
         }
      }

      Real m_amplitude;
};

#include "NamespaceFooter.H"

#endif

