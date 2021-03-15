#ifndef _SCATTERING_FACTORY_
#define _SCATTERING_FACTORY_

#include "ParmParse.H"
#include "Scattering.H"

#include "HardSphere.H"
#include "VariableHardSphere.H"
#include "TakizukaAbe.H"

#include "NamespaceHeader.H"

/**
 * scattering factory class.
 * Allows calling routine to create a operator of abstract class
 * scattering models with no knowledge of the specific type.  In this way,
 * the calling code is insulated from modifications due to the addition of
 * new functions.
*/
class ScatteringFactory
{
   public:

      /// Constructor.
      /**
       */
      ScatteringFactory() {;}

      /// Destructor.
      /**
       */
      ~ScatteringFactory() {;}

      /// Create new Scattering object
      /**
       * Returns a dynamically allocated new Scattering object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       */
      inline RefCountedPtr<Scattering> create( ParmParse&   a_pp,
                                         const int&         a_verbosity = 0 );
};


inline
RefCountedPtr<Scattering> ScatteringFactory::create( ParmParse&  a_pp,
                                               const int&        a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "model", model );

   //std::string name( "foo" );

   RefCountedPtr<Scattering> result(NULL);
   if (model == "hs" || model=="HS") {
      result = RefCountedPtr<Scattering>( new HardSphere( a_pp, a_verbosity ) );
   }
   else if (model == "vhs" || model=="VHS") {
      result = RefCountedPtr<Scattering>( new VariableHardSphere( a_pp, a_verbosity ) );
   }
   else if (model == "ta" || model=="TA") {
      result = RefCountedPtr<Scattering>( new TakizukaAbe( a_pp, a_verbosity ) );
   }
   else {
      cout << "model = " << model << endl;
      MayDay::Error( "ScatteringFactory: Unimplemented scattering model selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
