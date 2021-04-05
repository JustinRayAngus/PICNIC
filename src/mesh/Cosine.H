#ifndef _COSINE_H_
#define _COSINE_H_

#include "GridFunction.H"
#include "Constants.H"

#include "NamespaceHeader.H"

/**
 * Cosine function class.
 * f = C + A \cdot \prod_{i=1}^D \cos(2\pi mode_i x_i / L_i + phi_i)
*/
class Cosine : public GridFunction
{
   public:

      Cosine( ParmParse& pp, const int& verbosity )
         : GridFunction(verbosity)
      {
         parseParameters( pp );
      }

      virtual ~Cosine() {;}

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Cosine grid function parameters:" << std::endl;
            std::cout << " f = C + A cdot prod_{i=1}^D cos(2pi mode_i x_i / L_i + phi_i)" << std::endl;
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
                           const DomainGrid&  a_mesh,
                           const FArrayBox&  a_real_coords ) const
      {
         // f = C + A \cdot \prod_{i=1}^D \cos(2\pi mode_i x_i / L_i + phi_i)
       
         const RealVect& Xmin = a_mesh.getXmin();  
         const RealVect& Xmax = a_mesh.getXmax();  
 
         a_data.setVal(m_amplitude);
            
         const Box data_box = a_data.box();
         BoxIterator bit(data_box); 
         IntVect ib;
         Real X0, L0, arg, value;
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            for (int dir=0; dir<SpaceDim; dir++) {
               X0 = a_real_coords.get(ib,dir);
               L0 = Xmax[dir] - Xmin[dir];
               arg = Constants::TWOPI*m_mode[dir]*X0/L0 + m_phase[dir];
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

};

#include "NamespaceFooter.H"

#endif

