#ifndef _BENNETT_H_
#define _BENNETT_H_

#include "GridFunction.H"
#include "PicnicConstants.H"

#include "NamespaceHeader.H"

/**
 * Bennett B-field function class.
 * f(X) = 2*A*X[0]*a/(a^2 + X[0]^2)
*/
class Bennett : public GridFunction
{
   public:

      Bennett( ParmParse&  a_pp, 
       const int         a_verbosity )
         : GridFunction(a_verbosity)
      {
         CH_assert(SpaceDim==2);
         parseParameters( a_pp );
      }

      virtual ~Bennett() {;}
      
      // return function value at location a_X
      virtual Real getValue( const RealVect&  a_X ) const {
         
         Real X0, denom, value; 
         X0 = a_X[0];
         denom = m_a*m_a + X0*X0;
         value = 2.0*m_A*X0/denom;

         return value;

      }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Bennett grid function parameters:" << std::endl;
            std::cout << " f = 2*A*X[0]/(a^2 + X[0]^2)" << std::endl;
            std::cout << " a = " << m_a << std::endl;
            std::cout << " A = " << m_A << std::endl << std::endl;
         }
      }

   private:

      // prohibit copying
      Bennett( const Bennett& );
      Bennett& operator=( const Bennett& );

      void parseParameters( ParmParse& a_pp )
      {
   
         a_pp.get( "a", m_a );
         a_pp.get( "A", m_A );

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
      Real m_A;

};

#include "NamespaceFooter.H"

#endif

