#include "PicSpeciesBC.H"
#include "LoHiSide.H"
#include "BCUtils.H"

#include "CodeUnits.H"

#include "NamespaceHeader.H"


PicSpeciesBC::PicSpeciesBC( const std::string&  a_species_name,
                            const Real&         a_species_mass,
                            const DomainGrid&   a_mesh,
                            const CodeUnits&    a_units,
                            const int           a_verbosity )
   : m_species_name(a_species_name),
     m_species_mass(a_species_mass),
     m_mesh(a_mesh),
     m_verbosity(a_verbosity)
{
   // parse the input file
   string pp_prefix = "BC." + m_species_name;
   ParmParse pp(pp_prefix.c_str());
   parseParameters(pp,a_units);

   if(m_verbosity) printParameters();

}

PicSpeciesBC::~PicSpeciesBC()
{
   m_bdry_name.resize(0);
   m_bc_type.resize(0);
   m_bc_binary.resize(0); 
   m_periodic_bc_binary.resize(0); 
   m_outflow_list_vector.resize(0);
   m_inflow_list_vector.resize(0);
}

void PicSpeciesBC::parseParameters( ParmParse&  a_pp,
                              const CodeUnits&  a_units )
{
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   m_bc_type.resize(bdry_layout.size());
   m_bdry_name.resize(bdry_layout.size());
 
   m_outflow_list_vector.resize(bdry_layout.size());  
   m_inflow_list_vector.resize(bdry_layout.size());  
 
   m_inflow_bc_ptr_vect.resize(bdry_layout.size());

   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      m_bdry_name[b] = this_bdry_layout.name();
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[b];
      ParmParse fpp( prefix.c_str() );
      m_bc_type[b] = "outflow";
      if( fpp.contains("type") ) {
         fpp.query( "type", m_bc_type[b] );
         CH_assert( m_bc_type[b] == "symmetry" ||    
                    m_bc_type[b] == "outflow"  ||
                    m_bc_type[b] == "inflow_outflow" );
         if(m_bc_type[b]=="inflow_outflow") {
            InflowBC* inflow_bc = new InflowBC( b, m_bdry_name[b], 
                                                this_bdry_layout.dir(), this_bdry_layout.side(),
                                                m_species_name, m_species_mass, m_mesh, a_units );
            m_inflow_bc_ptr_vect[b] = InflowBCPtr(inflow_bc);
         }
      }
   }
   
   // check if this processor includes physical boundaries
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const Box& domain_box = domain.domainBox();
   
   m_bc_binary.resize(bdry_layout.size());
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const Side::LoHiSide bdry_side( this_bdry_layout.side() );
  
      bool is_bdry_box(false); 
      DataIterator dit(grids);
      for(dit.begin(); dit.ok(); ++dit) {
         is_bdry_box |= BCUtils::touchesPhysicalBoundary( domain_box, grids[dit],
                                                           bdry_dir, bdry_side );
      }
      m_bc_binary[b] = is_bdry_box;
   }
   
   const BoundaryBoxLayoutPtrVect& periodic_bdry_layout = m_mesh.getPeriodicBoundaryLayout();
   m_periodic_bc_binary.resize(periodic_bdry_layout.size());
   for (int b(0); b<periodic_bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(periodic_bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const Side::LoHiSide bdry_side( this_bdry_layout.side() );
  
      bool is_bdry_box(false); 
      DataIterator dit(grids);
      for(dit.begin(); dit.ok(); ++dit) {
         is_bdry_box |= BCUtils::touchesPhysicalBoundary( domain_box, grids[dit],
                                                           bdry_dir, bdry_side );
      }
      m_periodic_bc_binary[b] = is_bdry_box;
      if(is_bdry_box) {
      //   cout << "JRA: periodic_bc: procID = " << procID() << endl;
      //   cout << "JRA: periodic_bc: dir  = " << bdry_dir << endl;
      //   cout << "JRA: periodic_bc: side = " << bdry_side << endl;
      }
   }
   
}

void PicSpeciesBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "PicSpeciesBC =======================================" << std::endl;
      std::cout << "  " << m_species_name << std::endl;
      for (int i(0); i<m_bc_type.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << m_bc_type[i] << std::endl;
         if(m_bc_type[i]=="inflow_outflow") {
            if(m_inflow_bc_ptr_vect[i]!=NULL) m_inflow_bc_ptr_vect[i]->printParameters();
         }
      }
      std::cout << "===============================================" << std::endl << std::endl;
   }
}

void PicSpeciesBC::apply( List<JustinsParticle>&  a_outcast_list,
                    const bool&                   a_intermediate_advance,
                    const Real&                   a_time )
{
   CH_TIME("PicSpeciesBC::apply()");
   
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());

   // first enforce periodic BCs
   for (int dir(0); dir<SpaceDim; dir++) {
      if(domain.isPeriodic(dir)) {
         enforcePeriodic( a_outcast_list, dir, Xmin[dir], Xmax[dir]); 
      }
   }
   
   // second, loop over non-periodic boundaries and apply BCs
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];
      List<JustinsParticle>& inflow_list = m_inflow_list_vector[b];
      List<JustinsParticle>& outflow_list = m_outflow_list_vector[b];
      if(bdry_side==0 && m_bc_binary[b]) applyBC_Lo( a_outcast_list, inflow_list, outflow_list,
                                                     this_bc, bdry_dir, Xmin[bdry_dir], a_intermediate_advance);      
      if(bdry_side==1 && m_bc_binary[b]) applyBC_Hi( a_outcast_list, inflow_list, outflow_list,
                                                     this_bc, bdry_dir, Xmax[bdry_dir], a_intermediate_advance);      
   }
   
}

void PicSpeciesBC::createInflowParticles( const Real&  a_time,
                                          const Real&  a_cnormDt )
{
   CH_TIME("PicSpeciesBC::createInflowParticles()");
   
   // loop over non-periodic boundaries and create inflow particles
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   const LevelData<EdgeDataBox>& Xphys_ec = m_mesh.getXec();

   for (int b(0); b<bdry_layout.size(); b++) {
      if(m_inflow_bc_ptr_vect[b]==NULL) continue;
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
         
      List<JustinsParticle>& inflow_list = m_inflow_list_vector[b];
      CH_assert(inflow_list.length()==0);
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         const FArrayBox& this_Xec( Xphys_ec[interior_dit][bdry_dir] );
        
         Box bdry_box = bdry_grids[dit];
         if(bdry_side==0) bdry_box.setSmall(bdry_dir,bdry_box.bigEnd(bdry_dir));
         if(bdry_side==1) bdry_box.setBig(bdry_dir,bdry_box.smallEnd(bdry_dir));

         if(m_bc_binary[b]) {
            m_inflow_bc_ptr_vect[b]->apply( inflow_list, this_Xec,
                                            bdry_box, a_cnormDt );
         }

      }
   }
       
}
   
void PicSpeciesBC::removeOutflowParticles()
{
   CH_TIME("PicSpeciesBC::removeOutflowParticles()");
   
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());

   // loop over non-periodic boundaries and remove particles that are out of bounds
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];
      List<JustinsParticle>&  outflow_list = m_outflow_list_vector[b];
      if(bdry_side==0 && m_bc_binary[b]) remove_Lo(outflow_list,bdry_dir,Xmin[bdry_dir]);
      if(bdry_side==1 && m_bc_binary[b]) remove_Hi(outflow_list,bdry_dir,Xmax[bdry_dir]);
      CH_assert(outflow_list.length()==0); 
   }
   
}

void PicSpeciesBC::repositionOutflowParticles( List<JustinsParticle>&  a_outcast_list )
{
   CH_TIME("PicSpeciesBC::repositionOutflowParticles()");
   
   // This is only called during iterative particle advance for implicit schemes.
   // loop over non-periodic boundaries and reposition particles that are out of bounds
   // to be just inside the physical domain   

   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect& dX(m_mesh.getdX());

   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      List<JustinsParticle>&  outflow_list = m_outflow_list_vector[b];
      if(bdry_side==0 && m_bc_binary[b]) { 
         ListIterator<JustinsParticle> lit(outflow_list);
         for(lit.begin(); lit.ok();) {
            JustinsParticle& this_particle = lit();
            RealVect& this_x = this_particle.position();
            if(this_x[bdry_dir] < Xmin[bdry_dir]) {
               this_x[bdry_dir] = Xmin[bdry_dir] + 1.0e-4*dX[bdry_dir];
               a_outcast_list.transfer(lit);
            }
            else {
               ++lit;
            } 
         }
      }
      if(bdry_side==1 && m_bc_binary[b]) { 
         ListIterator<JustinsParticle> lit(outflow_list);
         for(lit.begin(); lit.ok();) {
            JustinsParticle& this_particle = lit();
            RealVect& this_x = this_particle.position();
            if(this_x[bdry_dir] >= Xmax[bdry_dir]) {
               this_x[bdry_dir] = Xmax[bdry_dir] - 1.0e-4*dX[bdry_dir]; 
               a_outcast_list.transfer(lit);
            } 
            else {
               ++lit;
            } 
         }
      }
   }
   
}

void PicSpeciesBC::repositionInflowParticles( ParticleData<JustinsParticle>&  a_data )
{
   CH_TIME("PicSpeciesBC::repositionInflowParticles()");
   
   // This is only called during iterative particle advance for implicit schemes.
   // loop over non-periodic boundaries and move particles that have their
   // old position out of bounds (inflow particles) back to the inflow list

   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());

   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      
      List<JustinsParticle>& inflow_list = m_inflow_list_vector[b];
      
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         ListBox<JustinsParticle>& box_list = a_data[interior_dit];
         List<JustinsParticle>& valid_list = box_list.listItems();
         ListIterator<JustinsParticle> lit(valid_list);
      
         if(bdry_side==0 && m_bc_binary[b]) { 
            for(lit.begin(); lit.ok();) {
               JustinsParticle& this_particle = lit();
               const RealVect& this_x_old = this_particle.position_old();
               if(this_x_old[bdry_dir] < Xmin[bdry_dir]) {
                  RealVect& this_x = this_particle.position(); // this_x = Xhalf when this is called
                  this_x[bdry_dir] = 2.0*this_x[bdry_dir] - this_x_old[bdry_dir];
                  inflow_list.transfer(lit);
               } 
               else {
                  ++lit;
               } 
            }
         }
         if(bdry_side==1 && m_bc_binary[b]) { 
            for(lit.begin(); lit.ok();) {
               JustinsParticle& this_particle = lit();
               const RealVect& this_x_old = this_particle.position_old();
               if(this_x_old[bdry_dir] >= Xmax[bdry_dir]) {
                  RealVect& this_x = this_particle.position(); // this_x = Xhalf when this is called
                  this_x[bdry_dir] = 2.0*this_x[bdry_dir] - this_x_old[bdry_dir]; 
                  inflow_list.transfer(lit);
               } 
               else {
                  ++lit;
               } 
            }
         }

      }

   }
   
}

void PicSpeciesBC::enforcePeriodicForIterativeSolver( List<JustinsParticle>&  a_list)
{
   CH_TIME("PicSpeciesBC::enforcePeriodicForIterativeSolver");

   // WARNING. Do not update Xold here. I originally did this in order for the step 
   // norm calc to work. However, I get 2D numerical energy test with collisions 
   // crashing at step 7847 (time step = 0.1). It does not crashes with optimization on,
   // but don't update Xold here anyway!

   // this function is called from the iterative particle position/velocity solver
   // used to achieve a converged solution for each particle xp and vp. 
   // If this is not called, then a situation can occur where the particle lives
   // on a proc touching the lower boundary, but the old position is at the upper
   // boundary (or vice-versa). Then, the position update, which occurs before
   // interpolating the fields to the particle in an iterative loop can have a 
   // position outside the domain of the proc that the particle lives in and will crash.
   // Note that remap is not called until after the iterative loop.
    
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect& dX(m_mesh.getdX());
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getPeriodicBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      if(!m_periodic_bc_binary[b]) continue;
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int dir = this_bdry_layout.dir();
      const int side(this_bdry_layout.side());
      
      Real Lbox = Xmax[dir] - Xmin[dir];       
      if(side==0) {
         Real Lhi = Xmax[dir] - dX[dir];       
         ListIterator<JustinsParticle> lit(a_list);
         for(lit.begin(); lit.ok(); ++lit) {      
            JustinsParticle& this_particle = lit();
            RealVect& this_x = this_particle.position();
            if(this_x[dir] >= Lhi) this_x[dir] = this_x[dir] - Lbox;
         }
      }
      if(side==1) {
         Real Llo = Xmin[dir] + dX[dir];       
         ListIterator<JustinsParticle> lit(a_list);
         for(lit.begin(); lit.ok(); ++lit) {      
            JustinsParticle& this_particle = lit();
            RealVect& this_x = this_particle.position();
            if(this_x[dir] <= Llo) this_x[dir] = this_x[dir] + Lbox;
         }
      }
   }

}

void PicSpeciesBC::enforcePeriodic( List<JustinsParticle>&  a_list,
                              const int&      a_dir,
                              const Real&     a_leftEdge,
                              const Real&     a_rightEdge )
{

   CH_TIME("PicSpeciesBC::enforcePeriodic");

   Real Lbox = a_rightEdge - a_leftEdge;

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      //RealVect& this_xold = this_particle.position_old();
      
      if (this_x[a_dir] < a_leftEdge) {
         this_x[a_dir] = this_x[a_dir] + Lbox;
         //this_xold[a_dir] = this_xold[a_dir] + Lbox;
      }
      if (this_x[a_dir] >= a_rightEdge) {
         this_x[a_dir] = this_x[a_dir] - Lbox;
         //this_xold[a_dir] = this_xold[a_dir] - Lbox;
      }
   }

}

void PicSpeciesBC::applyBC_Lo( List<JustinsParticle>&  a_outcast_list,
                               List<JustinsParticle>&  a_inflow_list,
                               List<JustinsParticle>&  a_outflow_list,
                         const std::string&            a_bc_type,
                         const int&                    a_dir,
                         const Real&                   a_leftEdge,
                         const bool&                   a_intermediate_advance )
{

   CH_TIME("PicSpeciesBC::applyBC_Lo");
   
   if(a_bc_type=="symmetry") {
      symmetry_Lo(a_outcast_list,a_dir,a_leftEdge);
   }
   if(a_bc_type=="outflow") {
      outflow_Lo(a_outcast_list,a_outflow_list,a_dir,a_leftEdge);
   }
   if(a_bc_type=="inflow_outflow") {
      outflow_Lo(a_outcast_list,a_outflow_list,a_dir,a_leftEdge);
      inflow_Lo(a_outcast_list,a_inflow_list,a_dir,a_leftEdge,a_intermediate_advance);
   }

}

void PicSpeciesBC::applyBC_Hi( List<JustinsParticle>&  a_outcast_list,
                               List<JustinsParticle>&  a_inflow_list,
                               List<JustinsParticle>&  a_outflow_list,
                         const std::string&            a_bc_type,
                         const int&                    a_dir,
                         const Real&                   a_rightEdge,
                         const bool&                   a_intermediate_advance )
{

   CH_TIME("PicSpeciesBC::applyBC_Hi");
   
   if(a_bc_type=="symmetry") {
      symmetry_Hi(a_outcast_list,a_dir,a_rightEdge);
   }
   if(a_bc_type=="outflow") {
      outflow_Hi(a_outcast_list,a_outflow_list,a_dir,a_rightEdge);
   }
   if(a_bc_type=="inflow_outflow") {
      outflow_Hi(a_outcast_list,a_outflow_list,a_dir,a_rightEdge);
      inflow_Hi(a_outcast_list,a_inflow_list,a_dir,a_rightEdge,a_intermediate_advance);
   }

}

void PicSpeciesBC::symmetry_Lo( List<JustinsParticle>&  a_list,
                          const int&                    a_dir,
                          const Real&                   a_leftEdge )
{

   CH_TIME("PicSpeciesBC::symmetry_Lo");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_xold = this_particle.position_old();
      RealVect& this_x = this_particle.position();
      std::array<Real,3>&  this_vold = this_particle.velocity_old(); // actually beta
      std::array<Real,3>&  this_v = this_particle.velocity(); // actually beta
      
      if (this_x[a_dir] <= a_leftEdge) {
         this_x[a_dir] = 2.*a_leftEdge - this_x[a_dir];
         this_v[a_dir] = -this_v[a_dir];

         // update of old values needed for split/iterative particle advance
         this_xold[a_dir] = 2.*a_leftEdge - this_xold[a_dir];
         this_vold[a_dir] = -this_vold[a_dir];
      }

   }

}

void PicSpeciesBC::symmetry_Hi( List<JustinsParticle>&  a_list,
                          const int&                    a_dir,
                          const Real&                   a_rightEdge )
{

   CH_TIME("PicSpeciesBC::symmetry_Hi");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_xold = this_particle.position_old();
      RealVect& this_x = this_particle.position();
      std::array<Real,3>&  this_vold = this_particle.velocity_old(); // actually beta
      std::array<Real,3>&  this_v = this_particle.velocity(); // actually beta
      
      if (this_x[a_dir] >= a_rightEdge) {
         this_x[a_dir] = 2.*a_rightEdge - this_x[a_dir];
         this_v[a_dir] = -this_v[a_dir];
         
         // update of old values needed for split/iterative particle advance
         this_xold[a_dir] = 2.*a_rightEdge - this_xold[a_dir]; // needed for split/iterative particle advance
         this_vold[a_dir] = -this_vold[a_dir];
         
         // what do I do here if particle is right on upper boundary? valid domain is [Xmin,Xmax)
         // will it ever be right on to machine precision?
         if(this_x[a_dir]==a_rightEdge) this_x[a_dir] = 0.999999999*a_rightEdge; // total hack!
      }

   }

}

void PicSpeciesBC::outflow_Lo( List<JustinsParticle>&  a_outcast_list,
                               List<JustinsParticle>&  a_outflow_list,
                         const int&                    a_dir,
                         const Real&                   a_leftEdge )
{

   CH_TIME("PicSpeciesBC::outflow_Lo");

   ListIterator<JustinsParticle> lit(a_outcast_list);
   for(lit.begin(); lit.ok();) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_dir] < a_leftEdge) {
         a_outflow_list.transfer(lit); // this updates iterator
      }
      else {
         ++lit;
      }

   }

}

void PicSpeciesBC::outflow_Hi( List<JustinsParticle>&  a_outcast_list,
                               List<JustinsParticle>&  a_outflow_list,
                         const int&                    a_dir,
                         const Real&                   a_rightEdge )
{

   CH_TIME("PicSpeciesBC::outflow_Hi");

   ListIterator<JustinsParticle> lit(a_outcast_list);
   for(lit.begin(); lit.ok();) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_dir] >= a_rightEdge) {
         a_outflow_list.transfer(lit); // this updates iterator
      }
      else {
         ++lit;
      }

   }

}

void PicSpeciesBC::inflow_Lo( List<JustinsParticle>&  a_outcast_list,
                              List<JustinsParticle>&  a_inflow_list,
                        const int&                    a_dir,
                        const Real&                   a_leftEdge,
                        const bool&                   a_intermediate_advance )
{

   CH_TIME("PicSpeciesBC::inflow_Lo");

   ListIterator<JustinsParticle> lit(a_inflow_list);
      
   if(a_intermediate_advance) {
   
      for(lit.begin(); lit.ok();) {
         JustinsParticle& this_particle = lit();
         RealVect& this_xp= this_particle.position();
         Real this_xp_dir = this_xp[a_dir];
         Real this_xp_dir_old = this_particle.position_old(a_dir);
         Real this_xp_dir_half = (this_xp_dir + this_xp_dir_old)/2.0;
         if(this_xp_dir_half >= a_leftEdge) {
            this_xp[a_dir] = this_xp_dir_half;
            a_outcast_list.transfer(lit); // this updates iterator
         }
         else {
            ++lit;
         }
      }   

   }
   else {
      
      for(lit.begin(); lit.ok();) {
         JustinsParticle& this_particle = lit();
         RealVect& this_xp= this_particle.position();
         Real this_xp_dir = this_xp[a_dir];
         if(this_xp_dir >= a_leftEdge) {
            a_outcast_list.transfer(lit); // this updates iterator
         }
         else {
            ++lit;
         }
      }   

   }

}

void PicSpeciesBC::inflow_Hi( List<JustinsParticle>&  a_outcast_list,
                              List<JustinsParticle>&  a_inflow_list,
                        const int&                    a_dir,
                        const Real&                   a_rightEdge,
                        const bool&                   a_intermediate_advance )
{

   CH_TIME("PicSpeciesBC::inflow_Hi");

   if(a_intermediate_advance) {
   
      ListIterator<JustinsParticle> lit(a_inflow_list);
      for(lit.begin(); lit.ok();) {
         JustinsParticle& this_particle = lit();
         RealVect& this_xp= this_particle.position();
         Real this_xp_dir = this_xp[a_dir];
         Real this_xp_dir_old = this_particle.position_old(a_dir);
         Real this_xp_dir_half = (this_xp_dir + this_xp_dir_old)/2.0;
         if(this_xp_dir_half < a_rightEdge) {
            this_xp[a_dir] = this_xp_dir_half;
            a_outcast_list.transfer(lit); // this updates iterator
         }
         else {
            ++lit;
         }
      }   

   }
   else {
      
      ListIterator<JustinsParticle> lit(a_inflow_list);
      for(lit.begin(); lit.ok();) {
         JustinsParticle& this_particle = lit();
         RealVect& this_xp= this_particle.position();
         Real this_xp_dir = this_xp[a_dir];
         if(this_xp_dir < a_rightEdge) {
            a_outcast_list.transfer(lit); // this updates iterator
         }
         else {
            ++lit;
         }
      }   

   }

}

void PicSpeciesBC::remove_Lo( List<JustinsParticle>&  a_list,
                        const int&                    a_dir,
                        const Real&                   a_leftEdge )
{

   CH_TIME("PicSpeciesBC::remove_Lo");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok();) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_dir] < a_leftEdge) {
         a_list.remove(lit); // this updates iterator
      }
      else {
         ++lit;
      }

   }

}

void PicSpeciesBC::remove_Hi( List<JustinsParticle>&  a_list,
                        const int&                    a_dir,
                        const Real&                   a_rightEdge )
{

   CH_TIME("PicSpeciesBC::remove_Hi");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok();) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_dir] >= a_rightEdge) {
         a_list.remove(lit); // this updates iterator
      }
      else {
         ++lit;
      }

   }

}


#include "NamespaceFooter.H"
