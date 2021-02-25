#ifndef _PISTON_H_
#define _PISTON_H_

#include "SpecialOperator.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

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
            std::cout << " piston direction: "  << m_piston_dir  << std::endl;
            std::cout << " piston velocity: "  << m_piston_velocity  << std::endl << std::endl;
         }
      }
      
      virtual void applyOp( PicSpecies&  a_picSpecies, 
                      const DomainGrid&  a_mesh,
                      const Real&        a_dt ) const
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
         for(dit.begin(); dit.ok(); ++dit) {

            //ListBox<JustinsParticle>& box_list = Pdata[dit];
            List<JustinsParticle>& pList = Pdata[dit].listItems();
            ListIterator<JustinsParticle> li(pList);
            for(li.begin(); li.ok(); ++li) {

               RealVect&  xp = li().position();
               RealVect&  vp = li().velocity();

               // update particle position
               xp += vp*a_dt;

               // check to see if particle crosses piston boundary
               //
               if(xp[m_piston_dir]>=m_piston_position) { // specular reflection
                  Real dr = xp[m_piston_dir]-m_piston_position;
                  xp[m_piston_dir] = m_piston_position-dr;
                  //vp[m_piston_dir] = -2.0*abs(vp[m_piston_dir]+m_piston_velocity);
                  vp[m_piston_dir] = -(vp[m_piston_dir]-m_piston_velocity) + m_piston_velocity;
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
                        vp[dir] = -vp[dir];
                     }
                     if(xp[dir]>=Xmax[dir]) {
                        xp[dir] = 2.*Xmax[dir] - xp[dir];
                        vp[dir] = -vp[dir];
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
 
      virtual void updateOp( const DomainGrid&  a_mesh,
                             const Real&        a_dt )
      {
         const ProblemDomain& domain(a_mesh.getDomain());
         const RealVect& dX(a_mesh.getdX());
         const RealVect& Xmin(a_mesh.getXmin());
         
         m_piston_position = m_piston_position + a_dt*m_piston_velocity;
         if(m_piston_position<dX[m_piston_dir]) {
            m_piston_position = Xmin[m_piston_dir] + dX[m_piston_dir];
            m_piston_velocity = -m_piston_velocity;
            if(!procID()) cout << "piston too close to boundary !!! " << endl;
            exit(EXIT_FAILURE);
         }
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
         /*
         std::vector<Real> input_piston_position( SpaceDim );
         std::vector<Real> input_piston_velocity( SpaceDim );
         a_pp.getarr( "position", input_piston_position, 0, SpaceDim );
         a_pp.getarr( "velocity", input_piston_velocity, 0, SpaceDim );
         for (int dim=0; dim<SpaceDim; dim++)  {
            m_piston_position[dim] = input_piston_position[dim];
            m_piston_velocity[dim] = input_piston_velocity[dim];
         }
         */
         a_pp.query( "position", m_piston_position );
         a_pp.query( "velocity", m_piston_velocity );
         a_pp.query( "direction", m_piston_dir );
         
         if (m_verbosity) {
            printParameters();
         }
      }

      Real m_piston_position;
      Real m_piston_velocity;
      int m_piston_dir;

};

#include "NamespaceFooter.H"

#endif

