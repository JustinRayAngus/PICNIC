#ifndef _COSINE_H_
#define _COSINE_H_

#include "GridFunction.H"
#include "Constants.H"

#include "NamespaceHeader.H"

/**
 * Cosine function class.
 * f = C + A \cdot \prod_{i=1}^D \cos(2\pi mode_i x_i / L_i + phi_i*pi)
*/
class Cosine : public GridFunction
{
   public:

      Cosine( ParmParse&   a_pp, 
        const DomainGrid&  a_mesh, 
        const int          a_verbosity )
         : GridFunction(a_verbosity)
      {
         parseParameters( a_pp );
         const RealVect& Xmin = a_mesh.getXmin();  
         const RealVect& Xmax = a_mesh.getXmax();  
         m_L = Xmax - Xmin;
      }

      virtual ~Cosine() {;}
      
      // return function value at location a_X
      virtual Real getValue( const RealVect&  a_X ) const {
         
         Real X0, L0, arg; 
         Real value = m_amplitude;
         for (int dir=0; dir<SpaceDim; dir++) {
            X0 = a_X[dir];
            L0 = m_L[dir];
            arg = Constants::TWOPI*m_mode[dir]*X0/L0 + m_phase[dir]*Constants::PI;
            arg = std::fmod(arg,Constants::TWOPI); // cos(pi/2) != cos(3*pi/2) != cos(5*pi/2)...
            value = value*cos(arg);
         }
         value = value + m_constant;
         return value;

      }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Cosine grid function parameters:" << std::endl;
            std::cout << " f = C + A cdot prod_{i=1}^D cos(2pi mode_i x_i / L_i + phi_i*pi)" << std::endl;
            std::cout << " C   = " << m_constant << std::endl;
            std::cout << " A   = " << m_amplitude << std::endl;
            std::cout << " k   = " << m_mode << std::endl;
            std::cout << " phi = " << m_phase << std::endl << std::endl;
         }
      }

   private:

      // prohibit copying
      Cosine( const Cosine& );
      Cosine& operator=( const Cosine& );

      void parseParameters( ParmParse& a_pp )
      {
   
         a_pp.get( "constant", m_constant );
         a_pp.get( "amplitude", m_amplitude );

         //Vector<int> temp(SpaceDim);
         //temp.assign( 0 );
         std::vector<int> temp(SpaceDim,0);
         a_pp.getarr( "mode", temp, 0, SpaceDim );
         m_mode = IntVect( temp );

         Vector<Real> temp2(SpaceDim);
         temp2.assign( 0.0 );
         a_pp.getarr( "phase", temp2, 0, SpaceDim );
         m_phase = RealVect( temp2 );
         
         if (m_verbosity) {
            printParameters();
         }

      }

      // Set pointwise values at cell centers.
      virtual void setPointwise( FArrayBox&  a_data,
                           const FArrayBox&  a_real_coords ) const
      {
         // f = C + A \cdot \prod_{i=1}^D \cos(2\pi mode_i x_i / L_i + phi_i*pi)
       
         a_data.setVal(m_amplitude);
            
         const Box data_box = a_data.box();
         BoxIterator bit(data_box); 
         IntVect ib;
         Real X0, L0, arg, value;
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            for (int dir=0; dir<SpaceDim; dir++) {
               X0 = a_real_coords.get(ib,dir);
               L0 = m_L[dir];
               arg = Constants::TWOPI*m_mode[dir]*X0/L0 + m_phase[dir]*Constants::PI;
               arg = std::fmod(arg,Constants::TWOPI); // cos(pi/2) != cos(3*pi/2) != cos(5*pi/2)...
               value = a_data.get(ib,0);
               value = value*cos(arg);
               a_data.set(ib,0,value);
            }
         }
         a_data.plus(m_constant);
      
      }
      
      Real     m_constant;
      Real     m_amplitude;
      IntVect  m_mode;
      RealVect m_phase;
      RealVect m_L;

};

#include "NamespaceFooter.H"

#endif

