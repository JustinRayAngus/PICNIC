#ifndef _EXPONENTIAL_H_
#define _EXPONENTIAL_H_

#include "GridFunction.H"
#include "PicnicConstants.H"

#include "NamespaceHeader.H"

/**
 * Exponential grid function class.
 * f(X) = (A + B*exp(-x^p))/(C + x^(p-1)), where x = X[0]/a
*/
class Exponential : public GridFunction
{
   public:

      Exponential( ParmParse&  a_pp, 
       const int         a_verbosity )
         : GridFunction(a_verbosity)
      {
         parseParameters( a_pp );
      }

      virtual ~Exponential() {;}
      
      // return function value at location a_X
      virtual Real getValue( const RealVect&  a_X ) const {
         
         Real X0, denom, value; 
         X0 = a_X[0]/m_a;
         denom = m_C + std::pow(X0,m_p-1.);
         value = (m_A + m_B*std::exp(-std::pow(X0,m_p)))/denom;

         return value;

      }
      
      virtual Real getDerivative( const RealVect&  a_X,
                                  const int        a_dir ) const {
         CH_assert(a_dir<SpaceDim);
         
         Real deriv = 0.0;
         if(a_dir==0) {
            Real X0 = a_X[0]/m_a;
            Real denom = m_C + std::pow(X0,m_p-1.);
            Real numer1 = -denom*m_p*m_B*std::exp(-std::pow(X0,m_p-1.));
            Real numer2 = (m_A + m_B*std::exp(-std::pow(X0,m_p)))*(m_C + (m_p-1.0)*std::pow(X0,m_p-2.));
            deriv = (numer1 - numer2)/denom/denom;
         }

         return deriv;
      }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Exponential grid function parameters:" << std::endl;
            std::cout << " f = (A + B*exp(-x^p))/(C + x^(p-1)) with x = X[0]/a" << std::endl;
            std::cout << " a = " << m_a << std::endl;
            std::cout << " p = " << m_p << std::endl;
            std::cout << " A = " << m_A << std::endl;
            std::cout << " B = " << m_B << std::endl;
            std::cout << " C = " << m_C << std::endl << std::endl;
         }
      }

   private:

      // prohibit copying
      Exponential( const Exponential& );
      Exponential& operator=( const Exponential& );

      void parseParameters( ParmParse& a_pp )
      {
   
         a_pp.get( "a", m_a );
         a_pp.get( "p", m_p );
         a_pp.get( "A", m_A );
         a_pp.get( "B", m_B );
         a_pp.get( "C", m_C );

         CH_assert(m_p>=1);

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
      
      Real m_a, m_p;
      Real m_A, m_B, m_C;

};

#include "NamespaceFooter.H"

#endif

