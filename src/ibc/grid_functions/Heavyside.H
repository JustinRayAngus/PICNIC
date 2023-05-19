#ifndef _HEAVYSIDE_H_
#define _HEAVYSIDE_H_

#include "GridFunction.H"
#include "PicnicConstants.H"

#include "NamespaceHeader.H"

/**
 * Heavyside function class.
 * f(x) = C + A*H(x-x0)
 * H(y) = 0.0 for y < 0
 *      = 0.5 for y = 0
 *      = 1.0 for y > 0
*/
class Heavyside : public GridFunction
{
   public:

      Heavyside( ParmParse&   a_pp, 
        const DomainGrid&  a_mesh, 
        const int          a_verbosity )
         : GridFunction(a_verbosity)
      {
         parseParameters( a_pp );
    
         // define epsilon for equating x=x0
         const RealVect& dX = a_mesh.getdX();  
         //const Real epsilon = std::numeric_limits<Real>::epsilon();
         const Real epsilon = 1.0e-4;
         for (int dir=0; dir<SpaceDim; dir++) {
            m_epsilon[dir] = epsilon*dX[dir]; 
         }

      }

      virtual ~Heavyside() {;}
      
      // return function value at location a_X
      virtual Real getValue( const RealVect&  a_X ) const
      {   
         // f(x) = \prod_{i=1}^D ( C[i] + A[i]*H(x[i]-x0[i]) )
         Real arg, H;
         RealVect value;
         for (int dir=0; dir<SpaceDim; dir++) {
            arg = a_X[dir] - m_X0[dir];
            if(arg < 0.0) H = 0.0;
            else H = 1.0;
            if(abs(arg)<m_epsilon[dir]) H = 0.5;
            value[dir] = (m_C[dir] + m_A[dir]*H);
         }
         return value.product();

      }

      virtual Real getDerivative( const RealVect&  a_X,
                                  const int        a_dir ) const {
         CH_assert(a_dir<SpaceDim);
         return 0.0;         
      }

      virtual void printParameters()  const {
         if (procID()>0) return;
            std::cout << "Heavyside grid function parameters:" << std::endl;
            std::cout << " f = C + A*H(x-x0)" << std::endl;
         for (int dir=0; dir<SpaceDim; dir++) {
            std::cout << "dir = " << dir << std::endl;
            std::cout << " C   = " << m_C[dir] << std::endl;
            std::cout << " A   = " << m_A[dir] << std::endl;
            std::cout << " x0  = " << m_X0[dir] << std::endl;
         }
      }

   private:

      // prohibit copying
      Heavyside( const Heavyside& );
      Heavyside& operator=( const Heavyside& );

      void parseParameters( ParmParse& a_pp )
      {
   
         Vector<Real> temp(SpaceDim);
         temp.assign( 0.0 );
         a_pp.getarr( "constant", temp, 0, SpaceDim );
         m_C = RealVect( temp );

         Vector<Real> temp2(SpaceDim);
         temp2.assign( 0.0 );
         a_pp.getarr( "amplitude", temp2, 0, SpaceDim );
         m_A = RealVect( temp2 );
         
         Vector<Real> temp3(SpaceDim);
         temp3.assign( 0.0 );
         a_pp.getarr( "shift", temp3, 0, SpaceDim );
         m_X0 = RealVect( temp3 );
         
         if (m_verbosity) {
            printParameters();
         }

      }

      // Set pointwise values at cell centers.
      virtual void setPointwise( FArrayBox&  a_data,
                           const FArrayBox&  a_real_coords ) const
      {
         // f(x) = \prod_{i=1}^D ( C[i] + A[i]*H(x[i]-x0[i]) )
       
         a_data.setVal(0.0);
            
         const Box data_box = a_data.box();
         BoxIterator bit(data_box); 
         IntVect ib;
         RealVect this_X;
         Real value;
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            for(int dir=0; dir<SpaceDim; dir++) {
               this_X[dir] = a_real_coords.get(ib,dir);
            }
            value = getValue(this_X);
            a_data.set(ib,0,value);
         }
      
      }
      
      RealVect m_C;
      RealVect m_A;
      RealVect m_X0;
      RealVect m_epsilon;

};

#include "NamespaceFooter.H"

#endif

