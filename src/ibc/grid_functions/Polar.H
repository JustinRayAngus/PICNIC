#ifndef _POLAR_H_
#define _POLAR_H_

#include "GridFunction.H"
#include "Constants.H"

#include "NamespaceHeader.H"

/**
 * Polar function class.
 * f(x,z) = (A*x + B*z)/sqrt(x^2 + z^2)
*/
class Polar : public GridFunction
{
   public:

      Polar( ParmParse&  a_pp, 
       const int         a_verbosity )
         : GridFunction(a_verbosity)
      {
         CH_assert(SpaceDim==2);
         parseParameters( a_pp );
      }

      virtual ~Polar() {;}
      
      // return function value at location a_X
      virtual Real getValue( const RealVect&  a_X ) const {
         
         Real X0, Z0, R0, value; 
         X0 = a_X[0];
         Z0 = a_X[2];
         R0 = sqrt(X0*X0 + Z0*Z0);
         value = (m_a*X0 + m_b*Z0)/R0;
         return value;

      }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Polar grid function parameters:" << std::endl;
            std::cout << " f = (a*x + b*z)/sqrt(x^2 + z^2)" << std::endl;
            std::cout << " a = " << m_a << std::endl;
            std::cout << " b = " << m_b << std::endl << std::endl;
         }
      }

   private:

      // prohibit copying
      Polar( const Polar& );
      Polar& operator=( const Polar& );

      void parseParameters( ParmParse& a_pp )
      {
   
         a_pp.get( "a", m_a );
         a_pp.get( "b", m_b );

         if (m_verbosity) {
            printParameters();
         }

      }

      // Set pointwise values at cell centers.
      virtual void setPointwise( FArrayBox&  a_data,
                           const FArrayBox&  a_real_coords ) const
      {
         const Box data_box = a_data.box();
         BoxIterator bit(data_box); 
         IntVect ib;
         Real value;
         RealVect X;
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            X[0] = a_real_coords.get(ib,0);
            X[1] = a_real_coords.get(ib,1);
            value = getValue(X);
            a_data.set(ib,0,value);
         }
      }
      
      Real m_a;
      Real m_b;

};

#include "NamespaceFooter.H"

#endif

