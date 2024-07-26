#ifndef _GAUSSIANPULSE0D_H_
#define _GAUSSIANPULSE0D_H_

#include <cmath>

#include "TimeFunction.H"

#include "NamespaceHeader.H"

// exp(-((t-shift)^2/(2*width^2))*(A + B*sin(freq*t+phase)) 0D constant function class.

class GaussianPulse0D : public TimeFunction
{
   public:

      GaussianPulse0D( const Real  a_width,
                       const Real  a_shift,
                       const Real  a_freq,
                       const Real  a_phase,
                       const Real  a_A,
                       const Real  a_B,
                       const int   verbosity )
         : TimeFunction(verbosity), 
           m_width(a_width),
           m_shift(a_shift),
           m_freq(a_freq),
           m_phase(a_phase),
           m_A(a_A),
           m_B(a_B) {;}

      GaussianPulse0D( ParmParse& pp, const int& verbosity )
         : TimeFunction(verbosity),
           m_A(0.0),
           m_B(0.0)
      {
         parseParameters( pp );
      }

      virtual ~GaussianPulse0D() {;}

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "GaussianPulse0D function parameters:" << std::endl;
            std::cout << "  f(t) = exp(-0.5*(t-shift)^2/width^2)*(A + B*sin(freq*t+phase))" << std::endl;
            std::cout << "  width = " << m_width  << std::endl;
            std::cout << "  shift = " << m_shift  << std::endl;
            std::cout << "  freq  = " << m_freq  << std::endl;
            std::cout << "  phase = " << m_phase  << std::endl;
            std::cout << "  A = " << m_A  << std::endl;
            std::cout << "  B = " << m_B  << std::endl;
         }
      }
      
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const
      {
         const Real exp_arg = -0.5*std::pow((a_time-m_shift)/m_width,2);
         const Real sin_arg = m_freq*a_time + m_phase;
         //cout << "JRA: sin_arg = " << sin_arg << std::endl;
         //cout << "JRA: m_B = " << m_B << std::endl;;
         //cout << "JRA: sin(sin_arg) = " << std::sin(sin_arg) << std::endl;
         a_value = std::exp(exp_arg)*(m_A + m_B*std::sin(sin_arg));
      }

   private:

      // prohibit copying
      GaussianPulse0D( const GaussianPulse0D& );
      GaussianPulse0D& operator=( const GaussianPulse0D& );

      void parseParameters( ParmParse& pp )
      {
         pp.get( "width", m_width );
         pp.get( "shift", m_shift );
         pp.get( "freq",  m_freq );
         pp.get( "phase", m_phase );
         pp.get( "A", m_A );
         pp.get( "B", m_B );
      }

      Real m_width = 0.0;
      Real m_shift = 0.0;
      Real m_freq = 0.0;
      Real m_phase = 0.0;
      Real m_A = 0.0;
      Real m_B = 0.0;
};

#include "NamespaceFooter.H"

#endif

