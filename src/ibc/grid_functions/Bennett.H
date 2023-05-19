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
      
      virtual Real getDerivative( const RealVect&  a_X,
                                  const int        a_dir ) const {
         CH_assert(a_dir<SpaceDim);
         
         Real deriv = 0.0;
         if(a_dir==0) {
            Real X0 = a_X[0];
            Real denom = m_a*m_a + X0*X0;
            deriv = 2.0*m_A*(1.0 - 2.0*X0*X0/denom)/denom;
         }

         return deriv;
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
            for(int n=0; n<SpaceDim; n++) {
               X[n] = a_real_coords.get(ib,n);
            }
            value = getValue(X);
            a_data.set(ib,0,value);
         }
      }
      
      Real m_a;
      Real m_A;

};

#include "NamespaceFooter.H"

#endif

