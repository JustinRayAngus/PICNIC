#ifndef _Special_Operator_FACTORY_
#define _Special_Operator_FACTORY_

#include "ParmParse.H"
#include "SpecialOperator.H"

#include "Piston.H"

#include "NamespaceHeader.H"

/**
 * special operator factory class.
 * Allows calling routine to create a operator of abstract class
 * SpecialOperator with no knowledge of the specific type.  In this way,
 * the calling code is insulated from modifications due to the addition of
 * new functions.
*/
class SpecialOperatorFactory
{
   public:

      SpecialOperatorFactory() {;}

      ~SpecialOperatorFactory() {;}

      /// Create new SpecialOperator object
      /**
       * Returns a dynamically allocated new SpecialOperator object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       */
      inline RefCountedPtr<SpecialOperator> create( ParmParse& pp,
                                              const DomainGrid&  a_mesh,
                                              const CodeUnits&  a_units,
                                              const int& verbosity = 0 );
};


inline
RefCountedPtr<SpecialOperator> SpecialOperatorFactory::create(
   ParmParse&  a_pp,
   const DomainGrid&  a_mesh,
   const CodeUnits&  a_units,
   const int&  a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "model", model );

   //std::string name( "foo" );

   RefCountedPtr<SpecialOperator> result(NULL);
   if (model == "piston" || model == "piston" ) {
      result = RefCountedPtr<SpecialOperator>( new Piston( a_pp, a_mesh, a_units, a_verbosity ) );
   }
   else {
      MayDay::Error( "SpecialOperatorFactory: Unimplemented special operator type selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
