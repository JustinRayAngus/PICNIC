#ifndef _SINSQ0D_H_
#define _SINSQ0D_H_

#include <cmath>

#include "TimeFunction.H"

#include "NamespaceHeader.H"

// Sin^2 0D constant function class.

class SinSq0D : public TimeFunction
{
   public:

      SinSq0D( const Real&  a_rise_time,
               const Real&  a_amplitude, 
               const int&   verbosity )
         : TimeFunction(verbosity), 
           m_rise_time(a_rise_time),
           m_exponent(2.0),
           m_amplitude(a_amplitude) {;}

      SinSq0D( ParmParse& pp, const int& verbosity )
         : TimeFunction(verbosity), 
           m_rise_time(0.0),
           m_exponent(2.0),
           m_amplitude(0.0)
      {
         parseParameters( pp );
      }

      virtual ~SinSq0D() {;}

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Sinsq0D function parameters:" << std::endl;
            std::cout << " rise time = " << m_rise_time  << std::endl;
            std::cout << " exponent  = " << m_exponent  << std::endl << std::endl;
            std::cout << " amplitude = " << m_amplitude  << std::endl << std::endl;
         }
      }
      
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const
      {
         a_value = sin(3.14159/2.0*a_time/m_rise_time);
         //a_value = a_value*a_value*m_amplitude;
         a_value = pow(a_value,m_exponent)*m_amplitude;
      }

   private:

      // prohibit copying
      SinSq0D( const SinSq0D& );
      SinSq0D& operator=( const SinSq0D& );

      void parseParameters( ParmParse& pp )
      {
         pp.get( "rise_time", m_rise_time );
         m_exponent = 2.0;
         pp.query( "exponent", m_exponent );
         pp.get( "amplitude", m_amplitude );
         if (m_verbosity) {
            printParameters();
         }
      }

      Real m_rise_time;
      Real m_exponent;
      Real m_amplitude;
};

#include "NamespaceFooter.H"

#endif

