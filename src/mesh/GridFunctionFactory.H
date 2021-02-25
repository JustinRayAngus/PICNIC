#ifndef _Grid_Function_FACTORY_
#define _Grid_Function_FACTORY_

#include "ParmParse.H"
#include "GridFunction.H"

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

      /// Constructor.
      /**
       */
      GridFunctionFactory() {;}

      /// Destructor.
      /**
       */
      ~GridFunctionFactory() {;}

      /// Create new GridFunction object
      /**
       * Returns a dynamically allocated new GridFunction object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       *
       * The following represents a sample input entry:
       *
       * \verbatim
       * type = "localized"
       */
      inline RefCountedPtr<GridFunction> create( ParmParse& pp,
                                                 const int& verbosity = 0 );
};


inline
RefCountedPtr<GridFunction> GridFunctionFactory::create(
   ParmParse& a_pp,
   const int& a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "type", model );

   //std::string name( "foo" );

   RefCountedPtr<GridFunction> result(NULL);
   if (model == "constant") {
      result = RefCountedPtr<GridFunction>( new Constant( a_pp, a_verbosity ) );
   }
   else if (model == "zero") {
      result = RefCountedPtr<GridFunction>( new Constant( 0.0, a_verbosity ) );
   }
   else {
      MayDay::Error( "GridFunctionFactory: Unimplemented grid function type selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
