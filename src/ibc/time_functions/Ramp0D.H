#ifndef _RAMP0D_H_
#define _RAMP0D_H_

#include <cmath>

#include "TimeFunction.H"

#include "NamespaceHeader.H"

class Ramp0D : public TimeFunction
{
   public:

      Ramp0D( const Real  a_rise_time,
              const Real  a_steady_time,
              const Real  a_fall_time,
              const Real  a_amplitude,
              const int   verbosity )
         : TimeFunction(verbosity), 
           m_rise_time(a_rise_time),
           m_steady_time(a_steady_time),
           m_fall_time(a_fall_time),
           m_amplitude(a_amplitude) {;}

      Ramp0D( ParmParse& pp, const int verbosity )
         : TimeFunction(verbosity), 
           m_rise_time(0.0),
           m_steady_time(DBL_MAX),
           m_fall_time(0.0),
           m_amplitude(0.0)
      {
         parseParameters( pp );
      }

      virtual ~Ramp0D() {;}

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Ramp0D function parameters:" << std::endl;
            std::cout << "  rise time = " << m_rise_time  << std::endl;
            std::cout << "  steady_time = " << m_steady_time  << std::endl;
            std::cout << "  fall_time = " << m_fall_time  << std::endl;
            std::cout << "  amplitude = " << m_amplitude  << std::endl;
         }
      }
      
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const
      {
         a_value = a_time/m_rise_time*m_amplitude;
         a_value = std::min(a_value,m_amplitude);
         if(a_time>m_rise_time+m_steady_time) {
            a_value = m_amplitude*(1.0-(a_time-m_steady_time-m_rise_time)/m_fall_time);
            a_value = std::max(a_value,0.0);
         }
      }

   private:

      // prohibit copying
      Ramp0D( const Ramp0D& );
      Ramp0D& operator=( const Ramp0D& );

      void parseParameters( ParmParse& pp )
      {
         pp.get( "rise_time", m_rise_time );
         pp.get( "amplitude", m_amplitude );
         if(pp.contains("steady_time")) {
            pp.get( "steady_time", m_steady_time );
            pp.get( "fall_time", m_fall_time );
         }
      }

      Real m_rise_time;
      Real m_steady_time;
      Real m_fall_time;
      Real m_amplitude;
};

#include "NamespaceFooter.H"

#endif

