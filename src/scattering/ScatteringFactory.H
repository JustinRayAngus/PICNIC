#ifndef _SCATTERING_FACTORY_
#define _SCATTERING_FACTORY_

#include "ParmParse.H"
#include "Scattering.H"

#include "HardSphere.H"
#include "VariableHardSphere.H"
#include "TakizukaAbe.H"
#include "Coulomb.H"
#include "Fusion.H"
#include "Bremsstrahlung.H"
#include "IBremsstrahlung.H"
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
						      const ParmParse&  a_pp,
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
   else if (model == "fusion" || model == "FUSION") {
      result = RefCountedPtr<Scattering>( new Fusion( a_pp, a_verbosity ) );
   }
   else if (model == "bremsstrahlung") {
      result = RefCountedPtr<Scattering>( new Bremsstrahlung( a_pp, a_verbosity ) );
   }
   else if (model == "ibremsstrahlung") {
      result = RefCountedPtr<Scattering>( new IBremsstrahlung( a_pp, a_verbosity ) );
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
							    const ParmParse&  a_pp,
                                                            const int   a_verbosity )
{
   // query the Coulomb collision model type       
   std::string model = "TA";
   a_pp.query("model",model);

   RefCountedPtr<Scattering> result(NULL);
   if (model == "ta" || model=="TA") {
      Real Clog = 0.0;
      if(a_pp.contains("coulomb_logarithm")) {
         a_pp.get( "coulomb_logarithm", Clog );
         CH_assert(Clog>=2.0);
      }
      result = RefCountedPtr<Scattering>( new TakizukaAbe( a_sp1, a_sp2, Clog, a_verbosity ) );
   }
   else if (model == "coulomb" || model=="COULOMB") {
      result = RefCountedPtr<Scattering>( new Coulomb( a_sp1, a_sp2, a_pp, a_verbosity ) );
   }
   else {
      cout << "model = " << model << endl;
      MayDay::Error( "ScatteringFactory: Unimplemented scattering model selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
