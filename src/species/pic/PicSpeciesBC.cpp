#include "PicSpeciesBC.H"
#include "LoHiSide.H"
#include "BCUtils.H"
#include "FieldBCUtils.H"

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

   zeroDeltas();
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
                    m_bc_type[b] == "axis"  ||
                    m_bc_type[b] == "outflow"  ||
                    m_bc_type[b] == "inflow_outflow" );
         if(m_bc_type[b]=="inflow_outflow") {
            InflowBC* inflow_bc = new InflowBC( b, m_bdry_name[b], 
                                                this_bdry_layout.dir(), this_bdry_layout.side(),
                                                m_species_name, m_species_mass, m_mesh, a_units );
            m_inflow_bc_ptr_vect[b] = InflowBCPtr(inflow_bc);
         }
         if(m_bc_type[b]=="axis" && !m_mesh.axisymmetric()) m_bc_type[b] = "symmetry";
      }
   
      // below is needed for charge-conserving interp      
      if(m_bc_type[b]=="inflow_outflow" || m_bc_type[b]=="outflow") {
         const int bdry_dir = this_bdry_layout.dir();
         const int bdry_side(this_bdry_layout.side());
         if(bdry_side==0) m_isOutflowBC_lo[bdry_dir] = 1;
         if(bdry_side==1) m_isOutflowBC_hi[bdry_dir] = 1;
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
      //if(is_bdry_box) {
      //   cout << "JRA: periodic_bc: procID = " << procID() << endl;
      //   cout << "JRA: periodic_bc: dir  = " << bdry_dir << endl;
      //   cout << "JRA: periodic_bc: side = " << bdry_side << endl;
      //}
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

   // first, enforce periodic BCs
   for (int dir(0); dir<SpaceDim; dir++) {
      if(domain.isPeriodic(dir)) {
         // should be using m_periodic_bc_binary here...
         enforcePeriodic( a_outcast_list, dir, Xmin[dir], Xmax[dir]); 
      }
   }
   
   // second, loop over non-periodic boundaries and apply symmetry BCs
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];
      if(bdry_side==0 && m_bc_binary[b]) {
         if(this_bc=="axis") axis( a_outcast_list, bdry_dir );
         if(this_bc=="symmetry") symmetry_Lo( a_outcast_list, bdry_dir, Xmin[bdry_dir] );
      }
      if(bdry_side==1 && m_bc_binary[b]) {
         if(this_bc=="symmetry") symmetry_Hi( a_outcast_list, bdry_dir, Xmax[bdry_dir] );
      }
   }
   
   // third, loop over non-periodic boundaries again and apply inflow/outflow BCs
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];
      List<JustinsParticle>& inflow_list = m_inflow_list_vector[b];
      List<JustinsParticle>& outflow_list = m_outflow_list_vector[b];
      if(bdry_side==0 && m_bc_binary[b]) {
         if(this_bc=="outflow") {
            outflow_Lo( a_outcast_list, outflow_list, bdry_dir, Xmin[bdry_dir] );
         }
         if(this_bc=="inflow_outflow") {
            outflow_Lo( a_outcast_list, outflow_list, bdry_dir, Xmin[bdry_dir] );
            inflow_Lo( a_outcast_list, inflow_list, bdry_dir, Xmin[bdry_dir], a_intermediate_advance );
         }
      }
      if(bdry_side==1 && m_bc_binary[b]) {
         if(this_bc=="outflow") {
            outflow_Hi( a_outcast_list, outflow_list, bdry_dir, Xmax[bdry_dir] );
         }
         if(this_bc=="inflow_outflow") {
            outflow_Hi( a_outcast_list, outflow_list, bdry_dir, Xmax[bdry_dir] );
            inflow_Hi( a_outcast_list, inflow_list, bdry_dir, Xmax[bdry_dir], a_intermediate_advance );
         }
      }
   }
   
}

void PicSpeciesBC::applyToJ( LevelData<EdgeDataBox>&    a_J_inPlane,
                             LevelData<NodeFArrayBox>&  a_J_virtual )
{
   CH_TIME("PicSpeciesBC::applyToJ()");
 
   FieldBCUtils::applyToJ_PIC( a_J_inPlane,
                               a_J_virtual, 
                               m_mesh,
                               m_bc_type );

}

void PicSpeciesBC::applyToRho( LevelData<FArrayBox>&  a_Rho )
{
   CH_TIME("PicSpeciesBC::applyToRho() (cells)");
 
   // loop over non-periodic boundaries and apply BCs to Rho
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {

      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];

      
      for(DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         const Box bdry_box( bdry_grids[dit] );
          
         // collapse cell_box to 1 cell thick in bdry_dir direction                  
         Box cell_box = bdry_box;
         if(bdry_side==0) cell_box.setSmall(bdry_dir,cell_box.bigEnd(bdry_dir));
         if(bdry_side==1) cell_box.setBig(bdry_dir,cell_box.smallEnd(bdry_dir));
                  
         FArrayBox& this_Rho = a_Rho[interior_dit];
         Box dst_box = cell_box;
         Box src_box = cell_box;
         const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
         for(int n=0; n<nG; n++) {
            if(bdry_side==0) dst_box.shift(bdry_dir,1);
            if(bdry_side==1) dst_box.shift(bdry_dir,-1);
               this_Rho.plus(this_Rho,src_box,dst_box,0,0,1);
            if(bdry_side==0) src_box.shift(bdry_dir,-1);
            if(bdry_side==1) src_box.shift(bdry_dir,1);
         }

      }

   }
   
}

void PicSpeciesBC::applyToRho( LevelData<FluxBox>&  a_Rho )
{
   CH_TIME("PicSpeciesBC::applyToRho() (faces)");
 
   // loop over non-periodic boundaries and apply BCs to Rho
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {

      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];

      for(DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         const Box bdry_box( bdry_grids[dit] );
         const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
          
         for (int dir(0); dir<SpaceDim; dir++) {
            
            // collapse face_box to 1 cell thick in bdry_dir direction                  
            Box face_box = surroundingNodes(bdry_box,dir);
            if(bdry_side==0) face_box.setSmall(bdry_dir,face_box.bigEnd(bdry_dir));
            if(bdry_side==1) face_box.setBig(bdry_dir,face_box.smallEnd(bdry_dir));
                  
            FArrayBox& this_Rho( a_Rho[interior_dit][dir] );
            if(dir!=bdry_dir) { // cell centered
               Box dst_box = face_box;
               Box src_box = face_box;
               for (int n=0; n<nG; n++) {
                  if(bdry_side==0) dst_box.shift(bdry_dir,1);
                  if(bdry_side==1) dst_box.shift(bdry_dir,-1);
                  this_Rho.plus(this_Rho,src_box,dst_box,0,0,1);
                  if(bdry_side==0) src_box.shift(bdry_dir,-1);
                  if(bdry_side==1) src_box.shift(bdry_dir,1);
               }
            }
            else { // node centered
               if(this_bc!="axis") this_Rho.mult(2.0,face_box,0,this_Rho.nComp());
               Box dst_box = face_box;
               Box src_box = face_box;
               for(int n=0; n<nG; n++) {
                  dst_box.shift(bdry_dir,1 - 2*bdry_side);
                  src_box.shift(bdry_dir,2*bdry_side - 1);
                  this_Rho.plus(this_Rho,src_box,dst_box,0,0,this_Rho.nComp());
               }
            }
         }
               
      }
                  
   }
   
}

void PicSpeciesBC::applyToRho( LevelData<NodeFArrayBox>&  a_Rho )
{
   CH_TIME("PicSpeciesBC::applyToRho() (nodes)");
 
   // The particles near a boundary will have some rho deposited to the ghost cells.
   // This rho is added back to the interior. For symmetry boundaries, this is equivalent
   // to the rho that would come from mirror particles across the boundary. For other BCs,
   // adding the rho deposited to the ghost cells back to the interior amounts to using
   // a lower-order scheme near the boundary is makes is such that interior rho is conserved.
   // Right on the boundary, a factor of 2 is way to account for a 1/2 factor in dV.

   // loop over non-periodic boundaries and apply BCs to Rho
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {

      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];

      for(DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         const Box bdry_box( bdry_grids[dit] );
          
         // collapse node_box to 1 cell thick in bdry_dir direction                  
         Box node_box = surroundingNodes(bdry_box);
         if(bdry_side==0) node_box.setSmall(bdry_dir,node_box.bigEnd(bdry_dir));
         if(bdry_side==1) node_box.setBig(bdry_dir,node_box.smallEnd(bdry_dir));
               
         FArrayBox& this_Rho( a_Rho[interior_dit].getFab() );
         if(this_bc!="axis") this_Rho.mult(2.0,node_box,0,this_Rho.nComp());
         Box dst_box = node_box;
         Box src_box = node_box;
         const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
         for(int n=0; n<nG; n++) {
            dst_box.shift(bdry_dir,1 - 2*bdry_side);
            src_box.shift(bdry_dir,2*bdry_side - 1);
            this_Rho.plus(this_Rho,src_box,dst_box,0,0,this_Rho.nComp());
         }

      }

   }
   
}

void PicSpeciesBC::createInflowParticles( const Real&  a_time,
                                          const Real&  a_cnormDt,
                                          const ParticleData<JustinsParticle>&  a_data )
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
         const ListBox<JustinsParticle>& this_listbox = a_data[interior_dit];
         const List<JustinsParticle>& valid_list = this_listbox.listItems();
        
         Box bdry_box = bdry_grids[dit];
         if(bdry_side==0) bdry_box.setSmall(bdry_dir,bdry_box.bigEnd(bdry_dir));
         if(bdry_side==1) bdry_box.setBig(bdry_dir,bdry_box.smallEnd(bdry_dir));

         if(m_bc_binary[b]) {
            m_inflow_bc_ptr_vect[b]->apply( inflow_list, valid_list, this_Xec,
                                            bdry_box, a_time, a_cnormDt );
         }

      }
   }
       
}

void PicSpeciesBC::removeOutflowParticles( )
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
      if(bdry_side==0 && m_bc_binary[b]) remove_Lo( outflow_list, bdry_dir, Xmin[bdry_dir] );
      if(bdry_side==1 && m_bc_binary[b]) remove_Hi( outflow_list, bdry_dir, Xmax[bdry_dir] );
      CH_assert(outflow_list.length()==0); 
   }
   
}

void PicSpeciesBC::zeroDeltas( )
{
   for(int dir=0; dir<SpaceDim; dir++) {
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
      //
      m_delta_MassIn_lo[dir] = 0.0;
      m_delta_MassIn_hi[dir] = 0.0;
      m_delta_MomXIn_lo[dir] = 0.0;
      m_delta_MomXIn_hi[dir] = 0.0;
      m_delta_MomYIn_lo[dir] = 0.0;
      m_delta_MomYIn_lo[dir] = 0.0;
      m_delta_MomZIn_hi[dir] = 0.0;
      m_delta_MomZIn_hi[dir] = 0.0;
      m_delta_EnergyIn_lo[dir] = 0.0; 
      m_delta_EnergyIn_hi[dir] = 0.0; 
   }
}   

void PicSpeciesBC::repositionOutflowParticles( List<JustinsParticle>&  a_outcast_list )
{
   CH_TIME("PicSpeciesBC::repositionOutflowParticles()");
   
   // This is only called during iterative particle advance for implicit schemes.
   // loop over non-periodic boundaries and reposition outflow particles that are 
   // out of bounds to be just inside the physical domain and move them to the
   // outcast list.   

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

void PicSpeciesBC::injectInflowParticles( ParticleData<JustinsParticle>&  a_data )
{
   CH_TIME("PicSpeciesBC::injectInflowParticles()");
 
   // this is used by iterative implicit solvers to inject all inflow particles
   // at the beginning of the time step
 
   // will need to modify this for 2D if the spatial region along the
   // boundary is covered by multiple boxes per processor, in which case
   // a position check transverse to the boundary is needed to put the
   // particle in the correct box
 
   // loop over non-periodic boundaries and inject inflow particles
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      
      const std::string this_bc = m_bc_type[b];
      if(this_bc!="inflow_outflow") continue;
      if(!m_bc_binary[b]) continue;

      List<JustinsParticle>& inflow_list = m_inflow_list_vector[b];
      ListIterator<JustinsParticle> lit(inflow_list);
      
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         ListBox<JustinsParticle>& this_listbox = a_data[interior_dit];
         List<JustinsParticle>& valid_list = this_listbox.listItems();

         if(bdry_side==0) {
            for(lit.begin(); lit.ok();) {
               const Real& wp = lit().weight();
               //const std::array<Real,3>& vp = lit().velocity(); // actually beta
               const std::array<Real,3>& vp = lit().velocity_old(); // actually beta
               m_delta_MassIn_lo[bdry_dir] += wp;
               m_delta_MomXIn_lo[bdry_dir] += wp*vp[0];
               m_delta_MomYIn_lo[bdry_dir] += wp*vp[1];
               m_delta_MomZIn_lo[bdry_dir] += wp*vp[2];
               m_delta_EnergyIn_lo[bdry_dir] += 0.5*wp*(vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]);
               valid_list.transfer(lit);
            }
         }
         if(bdry_side==1) {
            for(lit.begin(); lit.ok();) {
               const Real& wp = lit().weight();
               //const std::array<Real,3>& vp = lit().velocity(); // actually beta
               const std::array<Real,3>& vp = lit().velocity_old(); // actually beta
               m_delta_MassIn_hi[bdry_dir] += wp;
               m_delta_MomXIn_hi[bdry_dir] += wp*vp[0];
               m_delta_MomYIn_hi[bdry_dir] += wp*vp[1];
               m_delta_MomZIn_hi[bdry_dir] += wp*vp[2];
               m_delta_EnergyIn_hi[bdry_dir] += 0.5*wp*(vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]);
               valid_list.transfer(lit);
            }
         }

      }

   }
   
}

void PicSpeciesBC::depositInflowOutflowJ( LevelData<EdgeDataBox>&    a_currentDensity,
                                          LevelData<NodeFArrayBox>&  a_currentDensity_virtual, 
                                    const MeshInterp&                a_meshInterp, 
                                    const InterpType                 a_interpJToGrid )
{
   CH_TIME("PicSpeciesBC::depositInflowOutflowJ()");
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );

      List<JustinsParticle>&  outflow_list = m_outflow_list_vector[b];
      List<JustinsParticle>&  inflow_list = m_inflow_list_vector[b];
   
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
      
         EdgeDataBox& J_inPlane = a_currentDensity[interior_dit];
#if CH_SPACEDIM<3
         FArrayBox& J_virtual = a_currentDensity_virtual[interior_dit].getFab();
#endif
         // deposit outflow particles
         a_meshInterp.depositCurrent( J_inPlane[0],
#if CH_SPACEDIM>=2
                                      J_inPlane[1],
#else
                                      J_virtual,
#endif
#if CH_SPACEDIM==3
                                      J_inPlane[2],
#else
                                      J_virtual,
#endif
                                      outflow_list,
                                      a_interpJToGrid );

         // deposit inflow particles
         a_meshInterp.depositCurrent( J_inPlane[0],
#if CH_SPACEDIM>=2
                                      J_inPlane[1],
#else
                                      J_virtual,
#endif
#if CH_SPACEDIM==3
                                      J_inPlane[2],
#else
                                      J_virtual,
#endif
                                      inflow_list,
                                      a_interpJToGrid );

      } 
      
   }
     
}

void PicSpeciesBC::depositInflowJ( LevelData<EdgeDataBox>&    a_inflowJ,
                                   LevelData<NodeFArrayBox>&  a_inflowJ_virtual, 
                             const MeshInterp&                a_meshInterp, 
                             const Real                       a_cnormDt )
{
   CH_TIME("PicSpeciesBC::depositInflowJ()");
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());

      List<JustinsParticle>&  inflow_list = m_inflow_list_vector[b];
   
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
      
         EdgeDataBox& J_inPlane = a_inflowJ[interior_dit];
#if CH_SPACEDIM<3
         FArrayBox& J_virtual = a_inflowJ_virtual[interior_dit].getFab();
#endif

         // deposit inflow particles
         a_meshInterp.depositInflowCurrent( J_inPlane[0],
#if CH_SPACEDIM>=2
                                            J_inPlane[1],
#else
                                            J_virtual,
#endif
#if CH_SPACEDIM==3
                                            J_inPlane[2],
#else
                                            J_virtual,
#endif
                                            inflow_list,
                                            bdry_dir, 
                                            bdry_side,
                                            a_cnormDt );

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
      RealVect& this_xold = this_particle.position_old();
      
      if (this_x[a_dir] < a_leftEdge) {
         this_x[a_dir] = this_x[a_dir] + Lbox;
         this_xold[a_dir] = this_xold[a_dir] + Lbox;
      }
      if (this_x[a_dir] >= a_rightEdge) {
         this_x[a_dir] = this_x[a_dir] - Lbox;
         this_xold[a_dir] = this_xold[a_dir] - Lbox;
      }
   }

}

void PicSpeciesBC::axis( List<JustinsParticle>&  a_list,
                   const int                     a_dir )
{
   CH_TIME("PicSpeciesBC::axis");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_xold = this_particle.position_old();
      RealVect& this_x = this_particle.position();
      std::array<Real,3>&  this_vold = this_particle.velocity_old(); // actually beta
      std::array<Real,3>&  this_v = this_particle.velocity(); // actually beta
      
      if (this_x[a_dir] < 0.0) {
         this_x[a_dir] = -this_x[a_dir];
         this_v[a_dir] = -this_v[a_dir];

         // update of old values needed for split/iterative particle advance
         this_xold[a_dir] = -this_xold[a_dir];
         this_vold[a_dir] = -this_vold[a_dir];
      }
      if (this_x[a_dir] == 0.0) {
         int th_dir = 1;  
         if(m_mesh.anticyclic()) th_dir = 2;
         Real Vsq = this_v[a_dir]*this_v[a_dir] + this_v[th_dir]*this_v[th_dir];
         this_v[a_dir] = -sqrt(Vsq);
         this_v[th_dir] = 0.0;
      }

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
         const JustinsParticle& this_particle = lit();
         const RealVect& this_xp = this_particle.position();
         if(this_xp[a_dir] >= a_leftEdge) {
            const Real& wp = lit().weight();
            const std::array<Real,3>& vpold = lit().velocity_old();
            m_delta_MassIn_lo[a_dir] += wp;
            m_delta_MomXIn_lo[a_dir] += wp*vpold[0];
            m_delta_MomYIn_lo[a_dir] += wp*vpold[1];
            m_delta_MomZIn_lo[a_dir] += wp*vpold[2];
            m_delta_EnergyIn_lo[a_dir] += 0.5*wp*(vpold[0]*vpold[0] + vpold[1]*vpold[1] + vpold[2]*vpold[2]);
            a_outcast_list.transfer(lit); // this updates iterator
         }
         else ++lit; // should never be here
      }   

   }
   else {
      
      for(lit.begin(); lit.ok();) {
         
         // convert vpbar and xpbar to vpnew and xpnew 
         JustinsParticle& this_particle = lit();
         RealVect& xp = this_particle.position();
         std::array<Real,3>& vp = lit().velocity();
         const RealVect& xpold = this_particle.position_old();
         const std::array<Real,3>& vpold = lit().velocity_old();
         for(int n=0; n<SpaceDim; n++) xp[n] = 2.0*xp[n] - xpold[n];
         for(int n=0; n<3; n++) vp[n] = 2.0*vp[n] - vpold[n];
        
         // update probes and transfer to outcast_list 
         if(xp[a_dir] >= a_leftEdge) {
            const Real& wp = lit().weight();
            m_delta_MassIn_lo[a_dir] += wp;
            m_delta_MomXIn_lo[a_dir] += wp*vpold[0];
            m_delta_MomYIn_lo[a_dir] += wp*vpold[1];
            m_delta_MomZIn_lo[a_dir] += wp*vpold[2];
            m_delta_EnergyIn_lo[a_dir] += 0.5*wp*(vpold[0]*vpold[0] + vpold[1]*vpold[1] + vpold[2]*vpold[2]);
            a_outcast_list.transfer(lit); // this updates iterator
         }
         else { // possible to be here when using sub-orbit model for inflow particles
            a_inflow_list.remove(lit); // this updates iterator
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
         const JustinsParticle& this_particle = lit();
         const RealVect& this_xp = this_particle.position();
         if(this_xp[a_dir] < a_rightEdge) {
            const Real& wp = lit().weight();
            const std::array<Real,3>& vpold = lit().velocity_old();
            m_delta_MassIn_hi[a_dir] += wp;
            m_delta_MomXIn_hi[a_dir] += wp*vpold[0];
            m_delta_MomYIn_hi[a_dir] += wp*vpold[1];
            m_delta_MomZIn_hi[a_dir] += wp*vpold[2];
            m_delta_EnergyIn_hi[a_dir] += 0.5*wp*(vpold[0]*vpold[0] + vpold[1]*vpold[1] + vpold[2]*vpold[2]);
            a_outcast_list.transfer(lit); // this updates iterator
         }
         else ++lit;
      }   

   }
   else {
      
      ListIterator<JustinsParticle> lit(a_inflow_list);
      for(lit.begin(); lit.ok();) {

         // convert vpbar and xpbar to vpnew and xpnew 
         JustinsParticle& this_particle = lit();
         RealVect& xp = this_particle.position();
         std::array<Real,3>& vp = lit().velocity();
         const RealVect& xpold = this_particle.position_old();
         const std::array<Real,3>& vpold = lit().velocity_old();
         for(int n=0; n<SpaceDim; n++) xp[n] = 2.0*xp[n] - xpold[n];
         for(int n=0; n<3; n++) vp[n] = 2.0*vp[n] - vpold[n];

         // update probes and transfer to outcast_list 
         if(xp[a_dir] < a_rightEdge) {
            const Real& wp = lit().weight();
            m_delta_MassIn_hi[a_dir] += wp;
            m_delta_MomXIn_hi[a_dir] += wp*vpold[0];
            m_delta_MomYIn_hi[a_dir] += wp*vpold[1];
            m_delta_MomZIn_hi[a_dir] += wp*vpold[2];
            m_delta_EnergyIn_hi[a_dir] += 0.5*wp*(vpold[0]*vpold[0] + vpold[1]*vpold[1] + vpold[2]*vpold[2]);
            a_outcast_list.transfer(lit); // this updates iterator
         }
         else { // possible to be here when using sub-orbit model for inflow particles
            a_inflow_list.remove(lit); // this updates iterator
         }

      }   

   }

}

void PicSpeciesBC::remove_Lo( List<JustinsParticle>&  a_list,
                        const int                     a_dir,
                        const Real                    a_leftEdge )
{

   CH_TIME("PicSpeciesBC::remove_Lo");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok();) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_dir] < a_leftEdge) {
         const Real& wp = lit().weight();
         const std::array<Real,3>& vp = lit().velocity(); // actually beta
         m_delta_MassOut_lo[a_dir] += wp;
         m_delta_MomXOut_lo[a_dir] += wp*vp[0];
         m_delta_MomYOut_lo[a_dir] += wp*vp[1];
         m_delta_MomZOut_lo[a_dir] += wp*vp[2];
         m_delta_EnergyOut_lo[a_dir] += 0.5*wp*(vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]);
         a_list.remove(lit); // this updates iterator
      }
      else {
         ++lit;
      }

   }

}

void PicSpeciesBC::remove_Hi( List<JustinsParticle>&  a_list,
                        const int                     a_dir,
                        const Real                    a_rightEdge )
{

   CH_TIME("PicSpeciesBC::remove_Hi");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok();) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_dir] >= a_rightEdge) {
         const Real& wp = lit().weight();
         const std::array<Real,3>& vp = lit().velocity(); // actually beta
         m_delta_MassOut_hi[a_dir] += wp;
         m_delta_MomXOut_hi[a_dir] += wp*vp[0];
         m_delta_MomYOut_hi[a_dir] += wp*vp[1];
         m_delta_MomZOut_hi[a_dir] += wp*vp[2];
         m_delta_EnergyOut_hi[a_dir] += 0.5*wp*(vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]);
         a_list.remove(lit); // this updates iterator
      }
      else {
         ++lit;
      }

   }

}


#include "NamespaceFooter.H"