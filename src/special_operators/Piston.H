#ifndef _PISTON_H_
#define _PISTON_H_

#include "SpecialOperator.H"
#include "JustinsParticle.H"
#include "ParticleData.H"
#include "Constants.H"

#include "NamespaceHeader.H"

/**
 * piston special operator class.
*/
class Piston : public SpecialOperator
{
   public:

      // Constructor
      //
      Piston( ParmParse& a_pp, const int& verbosity )
         : SpecialOperator(verbosity)
      {
         parseParameters( a_pp );
      }

      // Destructor.
      //
      virtual ~Piston() {;}

      // Print object parameters.
      //
      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Piston special operator parameters:" << std::endl;
            std::cout << " initial piston position: "  << m_piston_position  << std::endl;
            std::cout << " piston direction: " << m_piston_dir << std::endl;
            std::cout << " piston velocity:  " << m_piston_velocity << " m/s" << std::endl;
            std::cout << " piston dtFrac:    " << m_dtFrac  << std::endl;
         }
      }
      
      virtual void applyOp( PicSpecies&  a_picSpecies, 
                      const DomainGrid&  a_mesh,
                      const Real&        a_cnormDt ) const
      {
         CH_TIME("Piston::applyOp()");
   
         // This special operator does its own push, so check to make sure the species
         // motion is turned off
         //
         bool pMotion = a_picSpecies.motion();
         CH_assert(!pMotion);
    
         // get some Grid info
         //
         const ProblemDomain& domain(a_mesh.getDomain());
         const RealVect& dX(a_mesh.getdX());
   
         // compute domain extent and get total num parts and sim volume
         //
         const IntVect domainDimensions = domain.size();
         const RealVect& Xmin(a_mesh.getXmin());
         const RealVect& Xmax(a_mesh.getXmax());
         const RealVect Lbox = Xmax - Xmin; 
   
         // Each proc loops over its own boxes
         //
         ParticleData<JustinsParticle>& Pdata = a_picSpecies.partData();
         const BoxLayout& BL = Pdata.getBoxes();
         DataIterator dit(BL);
         const Real beta_piston = m_piston_velocity/Constants::CVAC;
         for(dit.begin(); dit.ok(); ++dit) {

            //ListBox<JustinsParticle>& box_list = Pdata[dit];
            List<JustinsParticle>& pList = Pdata[dit].listItems();
            ListIterator<JustinsParticle> li(pList);
            for(li.begin(); li.ok(); ++li) {

               RealVect&  xp = li().position();
               std::array<Real,3>&  betap = li().velocity(); // actually beta

               // update particle position
               //xp += betap*a_cnormDt;
               for(int dir=0; dir<SpaceDim; dir++) {
                  xp[dir] += betap[dir]*a_cnormDt;
               }

               // check to see if particle crosses piston boundary
               //
               if(xp[m_piston_dir]>=m_piston_position) { // specular reflection
                  Real dr = xp[m_piston_dir]-m_piston_position;
                  xp[m_piston_dir] = m_piston_position-dr;
                  //betap[m_piston_dir] = -2.0*abs(betap[m_piston_dir] + beta_piston);
                  betap[m_piston_dir] = -(betap[m_piston_dir]-beta_piston) + beta_piston;
               }
         
               // set particle boundary conditions here for now
               //
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
               //
               //
               //////////////////////////////////////////////////

            } // end loop over particle list
      
         } // end loop over boxes
         Pdata.gatherOutcast();
         Pdata.remapOutcast();
         CH_assert(Pdata.isClosed());

      }
      
      virtual Real specialOpsDt() const { return m_piston_dt; }
 
      virtual void updateOp( const DomainGrid&  a_mesh,
                             const CodeUnits&   a_units,
                             const Real&        a_cnormDt )
      {
         const ProblemDomain& domain(a_mesh.getDomain());
         const RealVect& dX(a_mesh.getdX());
         const RealVect& Xmin(a_mesh.getXmin());
         
         const Real beta_piston = m_piston_velocity/Constants::CVAC;
         m_piston_position = m_piston_position + a_cnormDt*beta_piston;
         if(m_piston_position<dX[m_piston_dir]) {
            m_piston_position = Xmin[m_piston_dir] + dX[m_piston_dir];
            m_piston_velocity = -m_piston_velocity;
            if(!procID()) cout << "piston too close to boundary !!! " << endl;
            exit(EXIT_FAILURE);
         }

         // set the time step for piston special op
         m_piston_dt = m_dtFrac*dX[m_piston_dir]/abs(beta_piston*a_units.CvacNorm());

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
         a_pp.get( "direction", m_piston_dir );
         a_pp.query( "dtFrac", m_dtFrac );
         
         if (m_verbosity) {
            printParameters();
         }
      }

      Real m_piston_position;
      Real m_piston_velocity;
      int m_piston_dir;
      Real m_dtFrac=0.5;
      Real m_piston_dt=DBL_MAX;

};

#include "NamespaceFooter.H"

#endif

