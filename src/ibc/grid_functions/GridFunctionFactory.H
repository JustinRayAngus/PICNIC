#ifndef _Grid_Function_FACTORY_
#define _Grid_Function_FACTORY_

#include "ParmParse.H"
#include "GridFunction.H"
#include "DomainGrid.H"

#include "Cosine.H"
#include "Heavyside.H"
#include "Polar.H"
#include "Bennett.H"
#include "Constant.H"

#include "NamespaceHeader.H"

/**
 * Configuration space function factory class.
 * Allows calling routine to create a function of abstract class
 * GridFunction with no knowledge of the specific type.  In this way,
 * the calling code is insulated from modifications due to the addition of
 * new functions.
*/
class GridFunctionFactory
{
   public:

      GridFunctionFactory() {;}

      ~GridFunctionFactory() {;}

      /// Create new GridFunction object
      /**
       * Returns a dynamically allocated new GridFunction object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       */
      inline RefCountedPtr<GridFunction> create( ParmParse&   a_pp,
                                           const DomainGrid&  a_mesh,
                                           const int          a_verbosity = 0 );
      
      inline RefCountedPtr<GridFunction> createZero() 
      {
         RefCountedPtr<GridFunction> result(NULL);
         result = RefCountedPtr<GridFunction>( new Constant( 0.0, 0 ) );
         return result;
      }

};


inline
RefCountedPtr<GridFunction> GridFunctionFactory::create(
         ParmParse&   a_pp,
   const DomainGrid&  a_mesh,
   const int          a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "type", model );

   RefCountedPtr<GridFunction> result(NULL);
   if (model == "cosine") {
      result = RefCountedPtr<GridFunction>( new Cosine( a_pp, a_mesh, a_verbosity ) );
   }
   else if (model == "heavyside") {
      result = RefCountedPtr<GridFunction>( new Heavyside( a_pp, a_mesh, a_verbosity ) );
   }
   else if (model == "polar") {
      result = RefCountedPtr<GridFunction>( new Polar( a_pp, a_verbosity ) );
   }
   else if (model == "bennett") {
      result = RefCountedPtr<GridFunction>( new Bennett( a_pp, a_verbosity ) );
   }
   else if (model == "constant") {
      result = RefCountedPtr<GridFunction>( new Constant( a_pp, a_verbosity ) );
   }
   else if (model == "zero") {
      result = RefCountedPtr<GridFunction>( new Constant( 0.0, a_verbosity ) );
   }
   else {
      MayDay::Error( "GridFunctionFactory: grid function type not implemented" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
