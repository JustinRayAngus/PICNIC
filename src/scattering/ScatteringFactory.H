#ifndef _SCATTERING_FACTORY_
#define _SCATTERING_FACTORY_

#include "ParmParse.H"
#include "Scattering.H"

#include "HardSphere.H"
#include "VariableHardSphere.H"
#include "TakizukaAbe.H"
#include "Coulomb.H"
#include "Elastic.H"
#include "Ionization.H"
#include "ChargeExchange.H"
#include "MonteCarloNull.H"

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
      inline RefCountedPtr<Scattering> create( ParmParse&    a_pp,
                                         const WeightMethod  a_weight_method,
                                         const int           a_verbosity = 0 );
      
      inline RefCountedPtr<Scattering> createCoulomb( const int   a_sp1,
                                                      const int   a_sp2,
                                                      const Real  a_Clog,
                                                      const std::string  a_model,
                                                      const std::string  a_angular_model,
                                                      const std::string  a_weight_method,
                                                      const int   a_verbosity = 0 );

};


inline
RefCountedPtr<Scattering> ScatteringFactory::create( ParmParse&    a_pp,
                                               const WeightMethod  a_weight_method,
                                               const int           a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "model", model );

   //std::string name( "foo" );

   RefCountedPtr<Scattering> result(NULL);
   if (model == "hs" || model=="HS") {
      result = RefCountedPtr<Scattering>( new HardSphere( a_pp, a_weight_method, a_verbosity ) );
   }
   else if (model == "vhs" || model=="VHS") {
      result = RefCountedPtr<Scattering>( new VariableHardSphere( a_pp, a_verbosity ) );
   }
   else if (model == "ta" || model=="TA") {
      result = RefCountedPtr<Scattering>( new TakizukaAbe( a_pp, a_verbosity ) );
   }
   else if (model == "coulomb" || model=="COULOMB") {
      result = RefCountedPtr<Scattering>( new Coulomb( a_pp,  a_weight_method, a_verbosity ) );
   }
   else if (model == "ela" || model=="ELA") {
      result = RefCountedPtr<Scattering>( new Elastic( a_pp, a_weight_method, a_verbosity ) );
   }
   //else if (model == "exc" || model=="EXC") {
   //   result = RefCountedPtr<Scattering>( new Excitation( a_pp, a_verbosity ) );
   //}
   //else if (model == "dis" || model=="DIS") {
   //   result = RefCountedPtr<Scattering>( new Dissociation( a_pp, a_verbosity ) );
   //}
   else if (model == "izn" || model=="IZN") {
      result = RefCountedPtr<Scattering>( new Ionization( a_pp, a_verbosity ) );
   }
   else if (model == "cex" || model=="CEX") {
      result = RefCountedPtr<Scattering>( new ChargeExchange( a_pp, a_weight_method, a_verbosity ) );
   }
   else if (model == "mcn" || model=="MCN") {
      result = RefCountedPtr<Scattering>( new MonteCarloNull( a_pp, a_weight_method, a_verbosity ) );
   }
   else {
      cout << "model = " << model << endl;
      MayDay::Error( "ScatteringFactory: Unimplemented scattering model selected" );
   }
   return result;
}

inline
RefCountedPtr<Scattering> ScatteringFactory::createCoulomb( const int   a_sp1,
                                                            const int   a_sp2,
                                                            const Real  a_Clog,
                                                            const std::string  a_model,
                                                            const std::string  a_angular_model,
                                                            const std::string  a_weight_method,
                                                            const int   a_verbosity )
{
   RefCountedPtr<Scattering> result(NULL);
   if (a_model == "ta" || a_model=="TA") {
      result = RefCountedPtr<Scattering>( new TakizukaAbe( a_sp1, a_sp2, a_Clog, a_verbosity ) );
   }
   else if (a_model == "coulomb" || a_model=="COULOMB") {
      result = RefCountedPtr<Scattering>( new Coulomb( a_sp1, a_sp2, a_Clog, a_angular_model, 
			                               a_weight_method, a_verbosity ) );
   }
   else {
      cout << "model = " << a_model << endl;
      MayDay::Error( "ScatteringFactory: Unimplemented scattering model selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
