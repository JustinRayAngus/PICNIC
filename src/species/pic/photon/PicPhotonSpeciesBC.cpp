#include "PicPhotonSpeciesBC.H"
#include "LoHiSide.H"
#include "BCUtils.H"
#include "FieldBCUtils.H"

#include "CodeUnits.H"

#include "NamespaceHeader.H"


PicPhotonSpeciesBC::PicPhotonSpeciesBC( const std::string&  a_species_name,
                                        const DomainGrid&   a_mesh,
                                        const CodeUnits&    a_units,
                                        const int           a_verbosity )
   : m_species_name(a_species_name),
     m_mesh(a_mesh),
     m_verbosity(a_verbosity)
{

   // parse the input file
   string pp_prefix = "BC." + m_species_name;
   ParmParse pp(pp_prefix.c_str());
   parseParameters(pp,a_units);

   zeroDeltas();
   if (m_verbosity) { printParameters(); }

}

PicPhotonSpeciesBC::~PicPhotonSpeciesBC()
{
   m_bdry_name.resize(0);
   m_bc_type.resize(0);
   m_bc_binary.resize(0); 
   m_periodic_bc_binary.resize(0); 
   m_outflow_list_vector.resize(0);
}

void PicPhotonSpeciesBC::parseParameters( ParmParse&  a_pp,
                              const CodeUnits&  a_units )
{
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   m_bc_type.resize(bdry_layout.size());
   m_bdry_name.resize(bdry_layout.size());
   m_outflow_list_vector.resize(bdry_layout.size());  
   
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

   for (int b(0); b<bdry_layout.size(); b++) {

      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      m_bdry_name[b] = this_bdry_layout.name();
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[b];
      ParmParse fpp( prefix.c_str() );

      m_bc_type[b] = "outflow";
      if (fpp.contains("type") ) {
         fpp.query( "type", m_bc_type[b] );
         CH_assert( m_bc_type[b] == "symmetry" ||    
                    m_bc_type[b] == "axis"  ||
                    m_bc_type[b] == "outflow" );
         if (m_bc_type[b]=="axis" && !m_mesh.axisymmetric()) { m_bc_type[b] = "symmetry"; }
         if (m_bc_type[b]=="symmetry" && m_mesh.axisymmetric()) {
            int bdry_side(this_bdry_layout.side());
            const RealVect& Xmin(m_mesh.getXmin());
            if (bdry_side==0 && Xmin[0]==0.0) { m_bc_type[b] = "axis"; }
         }
      }
   
      // below is needed for charge-conserving interp      
      if (m_bc_type[b]=="outflow") {
         const int bdry_dir = this_bdry_layout.dir();
         const int bdry_side(this_bdry_layout.side());
         if(bdry_side==0) { m_isOutflowBC_lo[bdry_dir] = (int)m_bc_binary[b]; }
         if(bdry_side==1) { m_isOutflowBC_hi[bdry_dir] = (int)m_bc_binary[b]; }
      }

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
   }
   
}

void PicPhotonSpeciesBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "PicPhotonSpeciesBC =======================================" << std::endl;
      std::cout << "  " << m_species_name << std::endl;
      for (int i(0); i<m_bc_type.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << m_bc_type[i] << std::endl;
      }
      std::cout << "===============================================" << std::endl << std::endl;
   }
}

void PicPhotonSpeciesBC::apply( List<PhotonParticle>&     a_outcast_list,
                          const bool                      a_intermediate_advance,
                          const Real                      a_time )
{
   CH_TIME("PicPhotonSpeciesBC::apply()");
   
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());

   // first, loop over non-periodic boundaries and apply symmetry BCs
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];
      if (bdry_side==0 && m_bc_binary[b]) {
         if (this_bc=="axis") { axis( a_outcast_list, bdry_dir ); }
         if (this_bc=="symmetry") { symmetry_Lo( a_outcast_list, bdry_dir, Xmin[bdry_dir] ); }
      }
      if (bdry_side==1 && m_bc_binary[b]) {
         if (this_bc=="symmetry") { symmetry_Hi( a_outcast_list, bdry_dir, Xmax[bdry_dir] ); }
      }
   }
   
   // second, loop over non-periodic boundaries again and apply outflow BCs
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];
      List<PhotonParticle>& outflow_list = m_outflow_list_vector[b];
      if (bdry_side==0 && m_bc_binary[b]) {
         if (this_bc=="outflow") {
            outflow_Lo( a_outcast_list, outflow_list, bdry_dir, Xmin[bdry_dir] );
         }
      }
      if (bdry_side==1 && m_bc_binary[b]) {
         if (this_bc=="outflow") {
            outflow_Hi( a_outcast_list, outflow_list, bdry_dir, Xmax[bdry_dir] );
         }
      }
   }
   
   // third, enforce periodic BCs
   for (int dir(0); dir<SpaceDim; dir++) {
      if (domain.isPeriodic(dir)) {
         // should be using m_periodic_bc_binary here...
         enforcePeriodic( a_outcast_list, dir, Xmin[dir], Xmax[dir]); 
      }
   }
   
   
}

void PicPhotonSpeciesBC::removeOutflowParticles()
{
   CH_TIME("PicPhotonSpeciesBC::removeOutflowParticles()");
   
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());

   // loop over non-periodic boundaries and remove particles that are out of bounds
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];
      List<PhotonParticle>&  outflow_list = m_outflow_list_vector[b];
      if (bdry_side==0 && m_bc_binary[b]) {
	remove_Lo( outflow_list, bdry_dir, Xmin[bdry_dir] );
      }
      if (bdry_side==1 && m_bc_binary[b]) {
	remove_Hi( outflow_list, bdry_dir, Xmax[bdry_dir] );
      }
      CH_assert(outflow_list.length()==0); 
   }
   
}

void PicPhotonSpeciesBC::zeroDeltas( )
{
   for (int dir=0; dir<SpaceDim; dir++) {
      m_delta_MassOut_lo[dir] = 0.0;
      m_delta_MassOut_hi[dir] = 0.0;
      m_delta_MomXOut_lo[dir] = 0.0;
      m_delta_MomXOut_hi[dir] = 0.0;
      m_delta_MomYOut_lo[dir] = 0.0;
      m_delta_MomYOut_lo[dir] = 0.0;
      m_delta_MomZOut_hi[dir] = 0.0;
      m_delta_MomZOut_hi[dir] = 0.0;
      m_delta_EnergyOut_lo[dir] = 0.0; 
      m_delta_EnergyOut_hi[dir] = 0.0; 
   }
}   

void PicPhotonSpeciesBC::enforcePeriodic( List<PhotonParticle>&  a_list,
                                    const int                    a_dir,
                                    const Real                   a_leftEdge,
                                    const Real                   a_rightEdge )
{

   CH_TIME("PicPhotonSpeciesBC::enforcePeriodic");

   Real Lbox = a_rightEdge - a_leftEdge;

   ListIterator<PhotonParticle> lit(a_list);
   for (lit.begin(); lit.ok(); ++lit) {

      PhotonParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();

      if (this_x[a_dir] < a_leftEdge) {
         this_x[a_dir] = this_x[a_dir] + Lbox;
      }
      if (this_x[a_dir] >= a_rightEdge) {
         this_x[a_dir] = this_x[a_dir] - Lbox;
      }
   }

}

void PicPhotonSpeciesBC::axis( List<PhotonParticle>&  a_list,
                          const int                   a_dir )
{
   CH_TIME("PicPhotonSpeciesBC::axis");

   ListIterator<PhotonParticle> lit(a_list);
   for (lit.begin(); lit.ok(); ++lit ) {

      PhotonParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      std::array<Real,3>&  this_v = this_particle.velocity();

      if (this_x[a_dir] <= 0.0) {
         this_x[a_dir] = -this_x[a_dir];
         this_v[a_dir] = -this_v[a_dir];
      }

   }

}

void PicPhotonSpeciesBC::symmetry_Lo( List<PhotonParticle>&  a_list,
                          const int                     a_dir,
                          const Real                    a_leftEdge )
{

   CH_TIME("PicPhotonSpeciesBC::symmetry_Lo");

   ListIterator<PhotonParticle> lit(a_list);
   for (lit.begin(); lit.ok(); ++lit) {

      PhotonParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      std::array<Real,3>&  this_v = this_particle.velocity();

      if (this_x[a_dir] <= a_leftEdge) {
         this_x[a_dir] = 2.*a_leftEdge - this_x[a_dir];
         this_v[a_dir] = -this_v[a_dir];
      }

   }

}

void PicPhotonSpeciesBC::symmetry_Hi( List<PhotonParticle>&  a_list,
                          const int                     a_dir,
                          const Real                    a_rightEdge )
{

   CH_TIME("PicPhotonSpeciesBC::symmetry_Hi");

   ListIterator<PhotonParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      PhotonParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      std::array<Real,3>&  this_v = this_particle.velocity();
      
      if (this_x[a_dir] >= a_rightEdge) {
         this_x[a_dir] = 2.*a_rightEdge - this_x[a_dir];
         this_v[a_dir] = -this_v[a_dir];

         // what do I do here if particle is right on upper boundary?
         // will it ever be right on to machine precision?
	 // valid domain is [Xmin,Xmax). See isEnclosed() in ListBoxI.H
         if (this_x[a_dir]==a_rightEdge) { this_x[a_dir] = 0.999999999*a_rightEdge; }
      }

   }

}

void PicPhotonSpeciesBC::outflow_Lo( List<PhotonParticle>&  a_outcast_list,
                               List<PhotonParticle>&  a_outflow_list,
                         const int                     a_dir,
                         const Real                    a_leftEdge )
{

   CH_TIME("PicPhotonSpeciesBC::outflow_Lo");

   ListIterator<PhotonParticle> lit(a_outcast_list);
   for (lit.begin(); lit.ok();) {

      PhotonParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if (this_x[a_dir] < a_leftEdge) {
         a_outflow_list.transfer(lit); // this updates iterator
      }
      else {
         ++lit;
      }

   }

}

void PicPhotonSpeciesBC::outflow_Hi( List<PhotonParticle>&  a_outcast_list,
                               List<PhotonParticle>&  a_outflow_list,
                         const int                     a_dir,
                         const Real                    a_rightEdge )
{

   CH_TIME("PicPhotonSpeciesBC::outflow_Hi");

   ListIterator<PhotonParticle> lit(a_outcast_list);
   for (lit.begin(); lit.ok();) {

      PhotonParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if (this_x[a_dir] >= a_rightEdge) {
         a_outflow_list.transfer(lit); // this updates iterator
      }
      else {
         ++lit;
      }

   }

}

void PicPhotonSpeciesBC::remove_Lo( List<PhotonParticle>&  a_list,
                        const int                     a_dir,
                        const Real                    a_leftEdge )
{

   CH_TIME("PicPhotonSpeciesBC::remove_Lo");

   Real gbpsq, energy;

   ListIterator<PhotonParticle> lit(a_list);
   for (lit.begin(); lit.ok();) {

      PhotonParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if (this_x[a_dir] < a_leftEdge) {

         const Real& wp = lit().weight();
         const std::array<Real,3>& up = lit().velocity();
         m_delta_MassOut_lo[a_dir] += wp;
         m_delta_MomXOut_lo[a_dir] += wp*up[0];
         m_delta_MomYOut_lo[a_dir] += wp*up[1];
         m_delta_MomZOut_lo[a_dir] += wp*up[2];
         
         gbpsq = up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
         energy = std::sqrt(gbpsq);
         m_delta_EnergyOut_lo[a_dir] += wp*energy;

         a_list.remove(lit); // this updates iterator

      }
      else {
         ++lit;
      }

   }

}

void PicPhotonSpeciesBC::remove_Hi( List<PhotonParticle>&  a_list,
                        const int                     a_dir,
                        const Real                    a_rightEdge )
{

   CH_TIME("PicPhotonSpeciesBC::remove_Hi");

   Real gbpsq, energy;

   ListIterator<PhotonParticle> lit(a_list);
   for (lit.begin(); lit.ok();) {

      PhotonParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if (this_x[a_dir] >= a_rightEdge) {

         const Real& wp = lit().weight();
         const std::array<Real,3>& up = lit().velocity();
         m_delta_MassOut_hi[a_dir] += wp;
         m_delta_MomXOut_hi[a_dir] += wp*up[0];
         m_delta_MomYOut_hi[a_dir] += wp*up[1];
         m_delta_MomZOut_hi[a_dir] += wp*up[2];

         gbpsq = up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
         energy = std::sqrt(gbpsq);
         m_delta_EnergyOut_hi[a_dir] += wp*energy;

         a_list.remove(lit); // this updates iterator

      }
      else {
         ++lit;
      }

   }

}


#include "NamespaceFooter.H"
