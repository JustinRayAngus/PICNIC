#ifndef _PISTON_H_
#define _PISTON_H_

#include "SpecialOperator.H"
#include "JustinsParticle.H"
#include "ParticleData.H"
#include "PicnicConstants.H"

#include "NamespaceHeader.H"

/**
 * piston special operator class.
*/
class Piston : public SpecialOperator
{
   public:

      Piston( ParmParse&   a_pp, 
        const DomainGrid&  a_mesh, 
        const CodeUnits&   a_units, 
        const int&         a_verbosity )
         : SpecialOperator(a_mesh, a_units, a_verbosity)
      {
         parseParameters( a_pp );
      }

      virtual ~Piston() {;}

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Piston special operator parameters:" << std::endl;
            std::cout << " initial piston position: "  << m_piston_position  << std::endl;
            std::cout << " piston direction: " << m_piston_dir << std::endl;
            std::cout << " piston velocity:  " << m_piston_velocity << " m/s" << std::endl;
            std::cout << " piston dtFrac:    " << m_dtFrac  << std::endl;
            if(m_inverse_r_velocity) std::cout << " piston velocity ~ 1/r " << std::endl;
         }
      }

      virtual void applyOp( PicSpeciesPtrVect&  a_pic_species_ptr_vect, 
                      const Real&               a_dt ) const
      {
         CH_TIME("Piston::applyOp()");
         
         const Real cnormDt = m_units.CvacNorm()*a_dt;
      
         for (int sp=0; sp<a_pic_species_ptr_vect.size(); sp++) {
            PicSpeciesPtr species(a_pic_species_ptr_vect[sp]);
   
            // This special operator does its own push, so check to make sure the species
            // motion is turned off
            bool pMotion = species->motion();
            CH_assert(!pMotion);
    
            const ProblemDomain& domain(m_mesh.getDomain());
   
            // compute domain extent and get total num parts and sim volume
            const IntVect domainDimensions = domain.size();
            const RealVect& Xmin(m_mesh.getXmin());
            const RealVect& Xmax(m_mesh.getXmax());
            const RealVect Lbox = Xmax - Xmin; 
 
            // these are used for cylindrical coords
            Real beta_x, beta_y, x_car, y_car; 
            int th_dir;

            // Each proc loops over its own boxes
            ParticleData<JustinsParticle>& Pdata = species->partData();
            const BoxLayout& BL = Pdata.getBoxes();
            DataIterator dit(BL);
            const Real beta_piston = m_piston_velocity/Constants::CVAC;
            for(dit.begin(); dit.ok(); ++dit) {

               List<JustinsParticle>& pList = Pdata[dit].listItems();
               ListIterator<JustinsParticle> li(pList);
               for(li.begin(); li.ok(); ++li) {

                  RealVect&  xp = li().position();
                  std::array<Real,3>&  betap = li().velocity(); // actually beta

                  // update particle position
                  for(int dir=0; dir<SpaceDim; dir++) xp[dir] += betap[dir]*cnormDt;

                  // correct velocity/positions for axisymmetric               
                  if(m_mesh.axisymmetric()) {
                     th_dir = 1;
                     if(m_mesh.anticyclic()) th_dir = 2;
                     beta_x = betap[0];
                     beta_y = betap[th_dir];
                     x_car = xp[0];
                     y_car = beta_y*cnormDt;
                     xp[0] = sqrt(x_car*x_car + y_car*y_car);
                     if(xp[0] > 0.0) {
                        betap[0] =  (x_car*beta_x + y_car*beta_y)/xp[0];
                        betap[th_dir] = (-y_car*beta_x + x_car*beta_y)/xp[0];
                     }
                  }

                  // check to see if particle crosses piston boundary
                  if(xp[m_piston_dir]>=m_piston_position) { // specular reflection
                     Real dr = xp[m_piston_dir]-m_piston_position;
                     if(m_mesh.axisymmetric()) { // conserve angular momentum
                        th_dir = 1;
                        if(m_mesh.anticyclic()) th_dir = 2;
                        Real rpold = xp[m_piston_dir];
                        betap[th_dir] = betap[th_dir]*rpold/(rpold-dr);
                     }
                     xp[m_piston_dir] = m_piston_position-dr;
                     betap[m_piston_dir] = -(betap[m_piston_dir]-beta_piston) + beta_piston;
                  }
         
                  // set particle boundary conditions here for now
                  for(int dir=0; dir<SpaceDim; dir++) {
                     if( domain.isPeriodic(dir) ) {
                        if( xp[dir]<Xmin[dir] ) xp[dir] = xp[dir] + Lbox[dir];
                        if( xp[dir]>=Xmax[dir] ) xp[dir] = xp[dir] - Lbox[dir];
                     }
                     else { // symmetry BCs
                        if(xp[dir]<=Xmin[dir]) {
                           xp[dir] = 2.*Xmin[dir] - xp[dir];
                           betap[dir] = -betap[dir];
                        }
                        if(xp[dir]>=Xmax[dir]) {
                           xp[dir] = 2.*Xmax[dir] - xp[dir];
                           betap[dir] = -betap[dir];
                        }
                     }
                  }

               } // end loop over particle list
      
            } // end loop over boxes
            Pdata.gatherOutcast();
            Pdata.remapOutcast();
            CH_assert(Pdata.isClosed());

         }

      }
      
      virtual Real specialOpsDt() const { return m_piston_dt; }
 
      virtual void updateOp( const Real&  a_dt )
      {
         const Real cnormDt = m_units.CvacNorm()*a_dt;
         const RealVect& dX(m_mesh.getdX());
         const RealVect& Xmin(m_mesh.getXmin());
         Real beta_piston;         
        
         // update piston velocity and position
         if(m_inverse_r_velocity) { // drp/dt = v0*r0/rp = vp
            Real rpn = m_piston_position;
            m_piston_position = sqrt(rpn*rpn + 2.0*m_vp0*m_rp0*cnormDt/Constants::CVAC);
            m_piston_velocity = m_vp0*m_rp0/m_piston_position;
            beta_piston = m_piston_velocity/Constants::CVAC;
         }
         else {
            beta_piston = m_piston_velocity/Constants::CVAC;
            m_piston_position = m_piston_position + cnormDt*beta_piston;
         }

         // check if piston is too close to the boundary        
         if(m_piston_position<dX[m_piston_dir]) {
            m_piston_position = Xmin[m_piston_dir] + dX[m_piston_dir];
            m_piston_velocity = -m_piston_velocity;
            if(!procID()) cout << "piston too close to boundary !!! " << endl;
            exit(EXIT_FAILURE);
         }

         // set the time step for piston special op
         m_piston_dt = m_dtFrac*dX[m_piston_dir]/abs(beta_piston*m_units.CvacNorm());

      }


   private:

      // prohibit copying
      //Piston( const Constant& );
      //Piston& operator=( const Constant& );

      /// Parse the input database for parameters.
      /**
       */
      void parseParameters( ParmParse& a_pp )
      {
         a_pp.get( "position", m_piston_position );
         a_pp.get( "velocity", m_piston_velocity );
         m_rp0 = m_piston_position;
         m_vp0 = m_piston_velocity;
         a_pp.get( "direction", m_piston_dir );
         a_pp.query( "dtFrac", m_dtFrac );
         a_pp.query( "inverse_r_velocity", m_inverse_r_velocity );
         if (m_verbosity) {
            printParameters();
         }
      }
  
      Real m_piston_position;
      Real m_piston_velocity;
      int m_piston_dir;
      Real m_dtFrac=0.5;
      Real m_piston_dt=DBL_MAX;
      bool m_inverse_r_velocity=false;

      Real m_rp0, m_vp0;
       
};

#include "NamespaceFooter.H"

#endif

