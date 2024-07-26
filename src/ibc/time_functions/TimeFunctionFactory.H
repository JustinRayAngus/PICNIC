#ifndef _TIME_FUNCTION_FACTORY_
#define _TIME_FUNCTION_FACTORY_

#include "ParmParse.H"
#include "TimeFunction.H"

#include "Ramp0D.H"
#include "Table0D.H"
#include "SinSq0D.H"
#include "GaussianPulse0D.H"
#include "Constant0D.H"

#include "NamespaceHeader.H"

/**
 * Time (0D) function factory class.
 * Allows calling routine to create a function of abstract class
 * TimeFunction with no knowledge of the specific type.  In this way,
 * the calling code is insulted from modifications due to the addition of
 * new functions.
*/
class TimeFunctionFactory
{
   public:

      TimeFunctionFactory() {;}

      ~TimeFunctionFactory() {;}

      inline RefCountedPtr<TimeFunction> create( ParmParse& pp,
                                                 const int& verbosity = 0 );
};


inline
RefCountedPtr<TimeFunction> TimeFunctionFactory::create(
   ParmParse& a_pp,
   const int& a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "type", model );

   RefCountedPtr<TimeFunction> result(NULL);
   if (model == "ramp" || model == "RAMP") {
      result = RefCountedPtr<TimeFunction>( new Ramp0D( a_pp, a_verbosity ) );
   }
   else if (model == "table" || model == "TABLE") {
      result = RefCountedPtr<TimeFunction>( new Table0D( a_pp, a_verbosity ) );
   }
   else if (model == "sinsq" || model == "SINSQ") {
      result = RefCountedPtr<TimeFunction>( new SinSq0D( a_pp, a_verbosity ) );
   }
   else if (model == "gaussian_pulse") {
      result = RefCountedPtr<TimeFunction>( new GaussianPulse0D( a_pp, a_verbosity ) );
   }
   else if (model == "constant") {
      result = RefCountedPtr<TimeFunction>( new Constant0D( a_pp, a_verbosity ) );
   }
   else if (model == "zero") {
      result = RefCountedPtr<TimeFunction>( new Constant0D( 0.0, a_verbosity ) );
   }
   else {
      MayDay::Error( "TimeFunctionFactory: Unimplemented time function type selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
