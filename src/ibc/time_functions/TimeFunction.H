#ifndef _TIME_FUNCTION_H_
#define _TIME_FUNCTION_H_

#include "REAL.H"

#include "NamespaceHeader.H"

class TimeFunction
{
   public:

      TimeFunction( const int& verbosity );
   
      virtual ~TimeFunction() {;}

      virtual void printParameters() const = 0;
   
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const = 0;
   
   protected:
   
      bool m_verbosity;
  
};

#include "NamespaceFooter.H"

#endif
