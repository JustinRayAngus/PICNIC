#include "PicSpeciesBC.H"
#include "LoHiSide.H"
#include "BCUtils.H"
#include "FieldBCUtils.H"

#include "CodeUnits.H"

#include "NamespaceHeader.H"


PicSpeciesBC::PicSpeciesBC( const std::string&  a_species_name,
                            const Real          a_species_mass,
                            const int           a_species_charge,
                            const InterpType&   a_interpRhoToGrid,
                            const DomainGrid&   a_mesh,
                            const CodeUnits&    a_units,
                            const int           a_verbosity )
   : m_species_name(a_species_name),
     m_species_mass(a_species_mass),
     m_species_charge(a_species_charge),
     m_interp_order(1),
     m_mesh(a_mesh),
     m_verbosity(a_verbosity)
{
   if(a_interpRhoToGrid==TSC) { m_interp_order = 2; }
   
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
         if(m_bc_type[b]=="symmetry" && m_mesh.axisymmetric()) {
            int bdry_side(this_bdry_layout.side());
            const RealVect& Xmin(m_mesh.getXmin());
            if(bdry_side==0 && Xmin[0]==0.0) m_bc_type[b] = "axis";
         }
      }
   
      // below is needed for charge-conserving interp      
      if(m_bc_type[b]=="inflow_outflow") {
         const int bdry_dir = this_bdry_layout.dir();
         const int bdry_side(this_bdry_layout.side());
         if(bdry_side==0) { m_isInflowBC_lo[bdry_dir] = (int)m_bc_binary[b]; }
         if(bdry_side==1) { m_isInflowBC_hi[bdry_dir] = (int)m_bc_binary[b]; }
      }
      if(m_bc_type[b]=="inflow_outflow" || m_bc_type[b]=="outflow") {
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

void PicSpeciesBC::apply( List<JustinsParticle>&     a_outcast_list,
                          LevelData<NodeFArrayBox>&  a_surfaceCharge,
                    const bool&                      a_intermediate_advance,
                    const Real&                      a_time )
{
   CH_TIME("PicSpeciesBC::apply()");
   
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
      if(bdry_side==0 && m_bc_binary[b]) {
         if(this_bc=="axis") axis( a_outcast_list, bdry_dir );
         if(this_bc=="symmetry") symmetry_Lo( a_outcast_list, bdry_dir, Xmin[bdry_dir] );
      }
      if(bdry_side==1 && m_bc_binary[b]) {
         if(this_bc=="symmetry") symmetry_Hi( a_outcast_list, bdry_dir, Xmax[bdry_dir] );
      }
   }
   
   // second, loop over non-periodic boundaries again and apply inflow/outflow BCs
   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
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
            for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
               const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
	       FArrayBox&  this_sigma( a_surfaceCharge[interior_dit].getFab() );
               const Box& bdry_box( bdry_grids[dit] );
	       extractSurfaceChargeLo( this_sigma, bdry_box, inflow_list, 
	   		               bdry_dir, Xmin[bdry_dir], a_intermediate_advance );
	    } 
            inflow_Lo( a_outcast_list, inflow_list, bdry_dir, Xmin[bdry_dir], a_intermediate_advance );
         }
      }
      if(bdry_side==1 && m_bc_binary[b]) {
         if(this_bc=="outflow") {
            outflow_Hi( a_outcast_list, outflow_list, bdry_dir, Xmax[bdry_dir] );
         }
         if(this_bc=="inflow_outflow") {
            outflow_Hi( a_outcast_list, outflow_list, bdry_dir, Xmax[bdry_dir] );
            for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
               const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
	       FArrayBox&  this_sigma( a_surfaceCharge[interior_dit].getFab() );
               const Box& bdry_box( bdry_grids[dit] );
	       extractSurfaceChargeHi( this_sigma, bdry_box, inflow_list, 
	   		               bdry_dir, Xmax[bdry_dir], a_intermediate_advance );
	    } 
            inflow_Hi( a_outcast_list, inflow_list, bdry_dir, Xmax[bdry_dir], a_intermediate_advance );
         }
      }
   }
   
   // third, enforce periodic BCs
   for (int dir(0); dir<SpaceDim; dir++) {
      if(domain.isPeriodic(dir)) {
         // should be using m_periodic_bc_binary here...
         enforcePeriodic( a_outcast_list, dir, Xmin[dir], Xmax[dir]); 
      }
   }
   
   
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
         const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
	 
	 // set cell_box for bdry and grow ghost cells in transverse direction
         Box cell_box( bdry_box );
         IntVect grow_vect = a_Rho.ghostVect();
         grow_vect[bdry_dir] = 0;
         cell_box.grow(grow_vect);
          
         // collapse cell_box to 1 cell thick in bdry_dir direction                  
         if(bdry_side==0) cell_box.setSmall(bdry_dir,cell_box.bigEnd(bdry_dir));
         if(bdry_side==1) cell_box.setBig(bdry_dir,cell_box.smallEnd(bdry_dir));
                  
         FArrayBox& this_Rho = a_Rho[interior_dit];
         Box dst_box = cell_box;
         Box src_box = cell_box;
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
            
	    // set face_box for bdry and grow ghost cells in transverse direction
            Box face_box = surroundingNodes(bdry_box,dir);
            IntVect grow_vect = a_Rho.ghostVect();
            grow_vect[bdry_dir] = 0;
            face_box.grow(grow_vect);
            
	    // collapse face_box to 1 cell thick in bdry_dir direction                  
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
   // a lower-order scheme near the boundary and makes it such that interior rho is conserved.

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
	       
	 // convert bdry_box to node_box and grow ghost cells in transverse direction
         Box node_box = surroundingNodes(bdry_box);
         IntVect grow_vect = a_Rho.ghostVect();
         grow_vect[bdry_dir] = 0;
         node_box.grow(grow_vect);

         // collapse node_box to 1 cell thick in bdry_dir direction                  
         if(bdry_side==0) node_box.setSmall(bdry_dir,node_box.bigEnd(bdry_dir));
         if(bdry_side==1) node_box.setBig(bdry_dir,node_box.smallEnd(bdry_dir));
         
         FArrayBox& this_Rho( a_Rho[interior_dit].getFab() );
         if(this_bc=="inflow_outflow" || this_bc=="outflow") {
            Box dst_box = node_box;
            Box src_box = node_box;
            if(bdry_side==0) src_box.shift(bdry_dir,-1);
            if(bdry_side==1) src_box.shift(bdry_dir,1);
	    const Real scale = 2.0;
            this_Rho.plus(this_Rho,src_box,dst_box,scale,0,0,1);
            if(bdry_side==0) dst_box.shift(bdry_dir,1);
            if(bdry_side==1) dst_box.shift(bdry_dir,-1);
            this_Rho.minus(this_Rho,src_box,dst_box,0,0,1);
	 }
	 else {
            Box dst_box = node_box;
            Box src_box = node_box;
            for(int n=0; n<nG; n++) {
               dst_box.shift(bdry_dir,1 - 2*bdry_side);
               src_box.shift(bdry_dir,2*bdry_side - 1);
               this_Rho.plus(this_Rho,src_box,dst_box,0,0,this_Rho.nComp());
            }
	 }

      }

   }
   
}

void PicSpeciesBC::applyToRhoInGhosts( LevelData<NodeFArrayBox>&  a_Rho )
{
   CH_TIME("PicSpeciesBC::applyToRhoInGhosts()");
 
   // this function fills the charge density container in the ghost cells
   // this is only needed when filtering is being used

   // loop over non-periodic boundaries and apply BCs to Rho
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {

      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const std::string this_bc = m_bc_type[b];
               
      std::string sub_bc_type = "zero";
      if(this_bc=="axis" || this_bc=="symmetry") sub_bc_type = "even"; 

      for(DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         const Box bdry_box( bdry_grids[dit] );
               
	 // convert to node type and grow box to include tranverse ghosts 
         const Box fill_box( bdry_box );
         Box fill_box_grown = surroundingNodes(fill_box);
         IntVect grow_vect = a_Rho.ghostVect();
         grow_vect[this_bdry_layout.dir()] = 0;
         fill_box_grown.grow(grow_vect);
               
	 // dont change value on face in setBC()
         const int ISIDE(this_bdry_layout.side());
         if(ISIDE==0) fill_box_grown.growHi(bdry_dir,-1);
         if(ISIDE==1) fill_box_grown.growLo(bdry_dir,-1);
            
         FArrayBox& this_Rho( a_Rho[interior_dit].getFab() );
	 for (int n(0); n<this_Rho.nComp(); n++) {
            BoundaryConditions::setBC( this_Rho,
                                       fill_box_grown,
                                       n,
                                       sub_bc_type,
                                       bdry_dir,
                                       this_bdry_layout.side() );
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
   const LevelData<NodeFArrayBox>& Xphys_nc = m_mesh.getXnc();

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
         const FArrayBox& this_Xnc( Xphys_nc[interior_dit].getFab() );
         const ListBox<JustinsParticle>& this_listbox = a_data[interior_dit];
         const List<JustinsParticle>& valid_list = this_listbox.listItems();
        
         Box bdry_box = bdry_grids[dit];
         if(bdry_side==0) { bdry_box.setSmall(bdry_dir,bdry_box.bigEnd(bdry_dir)); }
         if(bdry_side==1) { bdry_box.setBig(bdry_dir,bdry_box.smallEnd(bdry_dir)); }

         if(m_bc_binary[b]) {
            m_inflow_bc_ptr_vect[b]->apply( inflow_list, valid_list, this_Xnc,
                                            bdry_box, a_time, a_cnormDt );
         }

      }
   }
       
}

void PicSpeciesBC::removeOutflowParticles( LevelData<NodeFArrayBox>&  a_surfaceCharge )
{
   CH_TIME("PicSpeciesBC::removeOutflowParticles()");
   
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
      List<JustinsParticle>&  outflow_list = m_outflow_list_vector[b];
      if(bdry_side==0 && m_bc_binary[b]) {
        for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
           const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
	   FArrayBox&  this_sigma( a_surfaceCharge[interior_dit].getFab() );
           const Box& bdry_box( bdry_grids[dit] );
	   depositSurfaceChargeLo( this_sigma, bdry_box, outflow_list, 
			           bdry_dir, Xmin[bdry_dir] );
	}
	remove_Lo( outflow_list, bdry_dir, Xmin[bdry_dir] );
      }
      if(bdry_side==1 && m_bc_binary[b]) {
        for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
           const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
	   FArrayBox&  this_sigma( a_surfaceCharge[interior_dit].getFab() );
           const Box& bdry_box( bdry_grids[dit] );
	   depositSurfaceChargeHi( this_sigma, bdry_box, outflow_list, 
			           bdry_dir, Xmax[bdry_dir] );
	}
	remove_Hi( outflow_list, bdry_dir, Xmax[bdry_dir] );
      }
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

void PicSpeciesBC::injectInflowParticles( ParticleData<JustinsParticle>&  a_data,
                                          LevelData<NodeFArrayBox>&       a_surfaceCharge )
{
   CH_TIME("PicSpeciesBC::injectInflowParticles()");
 
   // this is used by iterative implicit solvers to inject all inflow particles
   // at the beginning of the time step
 
   // will need to modify this for 2D if the spatial region along the
   // boundary is covered by multiple boxes per processor, in which case
   // a position check transverse to the boundary is needed to put the
   // particle in the correct box
   
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
 
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
            
	 FArrayBox&  this_sigma( a_surfaceCharge[interior_dit].getFab() );
         const Box& bdry_box( bdry_grids[dit] );
   
         Real gbpsq, gammap;

         if(bdry_side==0) {
	    extractSurfaceChargeLo( this_sigma, bdry_box, inflow_list, 
	     		            bdry_dir, Xmin[bdry_dir], false );
            for(lit.begin(); lit.ok();) {

               const Real& wp = lit().weight();
               const std::array<Real,3>& up = lit().velocity_old();
               m_delta_MassIn_lo[bdry_dir] += wp;
               m_delta_MomXIn_lo[bdry_dir] += wp*up[0];
               m_delta_MomYIn_lo[bdry_dir] += wp*up[1];
               m_delta_MomZIn_lo[bdry_dir] += wp*up[2];

               gbpsq = up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
#ifdef RELATIVISTIC_PARTICLES
               gammap = sqrt(1.0 + gbpsq);
#else
               gammap = 1.0;
#endif
               m_delta_EnergyIn_lo[bdry_dir] += wp*gbpsq/(gammap + 1.0);
               valid_list.transfer(lit);

            }
         }
         if(bdry_side==1) {
	    extractSurfaceChargeHi( this_sigma, bdry_box, inflow_list, 
	     		            bdry_dir, Xmax[bdry_dir], false );
            for(lit.begin(); lit.ok();) {

               const Real& wp = lit().weight();
               const std::array<Real,3>& up = lit().velocity_old();
               m_delta_MassIn_hi[bdry_dir] += wp;
               m_delta_MomXIn_hi[bdry_dir] += wp*up[0];
               m_delta_MomYIn_hi[bdry_dir] += wp*up[1];
               m_delta_MomZIn_hi[bdry_dir] += wp*up[2];
               
               gbpsq = up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
#ifdef RELATIVISTIC_PARTICLES
               gammap = sqrt(1.0 + gbpsq);
#else
               gammap = 1.0;
#endif
               m_delta_EnergyIn_hi[bdry_dir] += wp*gbpsq/(gammap + 1.0);
               valid_list.transfer(lit);

            }
         }

      }

   }
   
}

void PicSpeciesBC::depositInflowOutflowJ( LevelData<EdgeDataBox>&    a_currentDensity,
                                          LevelData<NodeFArrayBox>&  a_currentDensity_virtual, 
                                    const MeshInterp&                a_meshInterp, 
                                    const InterpType                 a_interpJToGrid,
                                    const Real                       a_cnormDt,
                                    const bool                       a_from_explicit_solver )
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
                                      a_interpJToGrid,
                                      a_cnormDt,
                                      false,
                                      a_from_explicit_solver );

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
                                      a_interpJToGrid,
                                      a_cnormDt,
                                      false,
                                      a_from_explicit_solver );

      } 
      
   }
     
}

void PicSpeciesBC::enforcePeriodic( List<JustinsParticle>&  a_list,
                              const int                     a_dir,
                              const Real                    a_leftEdge,
                              const Real                    a_rightEdge )
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
      std::array<Real,3>&  this_vold = this_particle.velocity_old();
      std::array<Real,3>&  this_v = this_particle.velocity();
      
      if (this_x[a_dir] < 0.0) {
         this_x[a_dir] = -this_x[a_dir];
         this_v[a_dir] = -this_v[a_dir];

         // update of old values needed for split/iterative particle advance
         this_xold[a_dir] = -this_xold[a_dir];
         this_vold[a_dir] = -this_vold[a_dir];
      }
      if (this_x[a_dir] == 0.0) {
	 if(m_mesh.geomType()=="sph_R") {
            Real Vsq = this_v[0]*this_v[0] + this_v[1]*this_v[1] + this_v[2]*this_v[2];
	    this_v[0] = -std::sqrt(Vsq);
	    this_v[1] = 0.0;
	    this_v[2] = 0.0;
	 }
	 else{
            int th_dir = 1;  
            if(m_mesh.anticyclic()) th_dir = 2;
            Real Vsq = this_v[a_dir]*this_v[a_dir] + this_v[th_dir]*this_v[th_dir];
            this_v[a_dir] = -std::sqrt(Vsq);
            this_v[th_dir] = 0.0;
	 }
      }

   }

}

void PicSpeciesBC::symmetry_Lo( List<JustinsParticle>&  a_list,
                          const int                     a_dir,
                          const Real                    a_leftEdge )
{

   CH_TIME("PicSpeciesBC::symmetry_Lo");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_xold = this_particle.position_old();
      RealVect& this_x = this_particle.position();
      std::array<Real,3>&  this_vold = this_particle.velocity_old();
      std::array<Real,3>&  this_v = this_particle.velocity();
      
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
                          const int                     a_dir,
                          const Real                    a_rightEdge )
{

   CH_TIME("PicSpeciesBC::symmetry_Hi");

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_xold = this_particle.position_old();
      RealVect& this_x = this_particle.position();
      std::array<Real,3>&  this_vold = this_particle.velocity_old();
      std::array<Real,3>&  this_v = this_particle.velocity();
      
      if (this_x[a_dir] >= a_rightEdge) {
         this_x[a_dir] = 2.*a_rightEdge - this_x[a_dir];
         this_v[a_dir] = -this_v[a_dir];
         
         // update of old values needed for split/iterative particle advance
         this_xold[a_dir] = 2.*a_rightEdge - this_xold[a_dir]; // needed for split/iterative particle advance
         this_vold[a_dir] = -this_vold[a_dir];
         
         // what do I do here if particle is right on upper boundary?
         // will it ever be right on to machine precision?
	 // valid domain is [Xmin,Xmax). See isEnclosed() in ListBoxI.H
         if(this_x[a_dir]==a_rightEdge) this_x[a_dir] = 0.999999999*a_rightEdge; // total hack!
      }

   }

}

void PicSpeciesBC::outflow_Lo( List<JustinsParticle>&  a_outcast_list,
                               List<JustinsParticle>&  a_outflow_list,
                         const int                     a_dir,
                         const Real                    a_leftEdge )
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
                         const int                     a_dir,
                         const Real                    a_rightEdge )
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
                        const int                     a_dir,
                        const Real                    a_leftEdge,
                        const bool                    a_intermediate_advance )
{

   CH_TIME("PicSpeciesBC::inflow_Lo");
   
   Real gbpsq, gammap;

   ListIterator<JustinsParticle> lit(a_inflow_list);
      
   if(a_intermediate_advance) {

      for(lit.begin(); lit.ok();) {
         const JustinsParticle& this_particle = lit();
         const RealVect& this_xp = this_particle.position();
         if(this_xp[a_dir] >= a_leftEdge) {

            const Real& wp = lit().weight();
            const std::array<Real,3>& upold = lit().velocity_old();
            m_delta_MassIn_lo[a_dir] += wp;
            m_delta_MomXIn_lo[a_dir] += wp*upold[0];
            m_delta_MomYIn_lo[a_dir] += wp*upold[1];
            m_delta_MomZIn_lo[a_dir] += wp*upold[2];
            
            gbpsq = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
#ifdef RELATIVISTIC_PARTICLES
            gammap = sqrt(1.0 + gbpsq);
#else
            gammap = 1.0;
#endif
            m_delta_EnergyIn_lo[a_dir] += wp*gbpsq/(gammap + 1.0);
            a_outcast_list.transfer(lit); // this updates iterator

         }
         else ++lit; // should never be here
      }   

   }
   else {
      
      for(lit.begin(); lit.ok();) {
         
         // convert upbar and xpbar to upnew and xpnew 
         JustinsParticle& this_particle = lit();
         RealVect& xp = this_particle.position();
         std::array<Real,3>& up = lit().velocity();
         const RealVect& xpold = this_particle.position_old();
         const std::array<Real,3>& upold = lit().velocity_old();
         for(int n=0; n<SpaceDim; n++) xp[n] = 2.0*xp[n] - xpold[n];
         for(int n=0; n<3; n++) up[n] = 2.0*up[n] - upold[n];
        
         // update probes and transfer to outcast_list 
         if(xp[a_dir] >= a_leftEdge) {

            const Real& wp = lit().weight();
            m_delta_MassIn_lo[a_dir] += wp;
            m_delta_MomXIn_lo[a_dir] += wp*upold[0];
            m_delta_MomYIn_lo[a_dir] += wp*upold[1];
            m_delta_MomZIn_lo[a_dir] += wp*upold[2];
            
            gbpsq = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
#ifdef RELATIVISTIC_PARTICLES
            gammap = sqrt(1.0 + gbpsq);
#else
            gammap = 1.0;
#endif
            m_delta_EnergyIn_lo[a_dir] += wp*gbpsq/(gammap + 1.0);
            a_outcast_list.transfer(lit); // this updates iterator

         }
         else { // possible to be here when using sub-orbit model for inflow particles
	    //cout << "JRA: reflected particle with:" << endl;
	    //cout << "xpnew[0] = " << xp[0] << endl;   
	    //cout << "upnew[0] = " << up[0] << endl;   
            a_inflow_list.remove(lit); // this updates iterator
         }

      }   

   }

}

void PicSpeciesBC::inflow_Hi( List<JustinsParticle>&  a_outcast_list,
                              List<JustinsParticle>&  a_inflow_list,
                        const int                     a_dir,
                        const Real                    a_rightEdge,
                        const bool                    a_intermediate_advance )
{

   CH_TIME("PicSpeciesBC::inflow_Hi");
   
   Real gbpsq, gammap;

   if(a_intermediate_advance) {
   
      ListIterator<JustinsParticle> lit(a_inflow_list);
      for(lit.begin(); lit.ok();) {
         const JustinsParticle& this_particle = lit();
         const RealVect& this_xp = this_particle.position();
         if(this_xp[a_dir] < a_rightEdge) {

            const Real& wp = lit().weight();
            const std::array<Real,3>& upold = lit().velocity_old();
            m_delta_MassIn_hi[a_dir] += wp;
            m_delta_MomXIn_hi[a_dir] += wp*upold[0];
            m_delta_MomYIn_hi[a_dir] += wp*upold[1];
            m_delta_MomZIn_hi[a_dir] += wp*upold[2];
            
            gbpsq = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
#ifdef RELATIVISTIC_PARTICLES
            gammap = sqrt(1.0 + gbpsq);
#else
            gammap = 1.0;
#endif
            m_delta_EnergyIn_hi[a_dir] += wp*gbpsq/(gammap + 1.0);
            a_outcast_list.transfer(lit); // this updates iterator

         }
         else ++lit;
      }   

   }
   else {
      
      ListIterator<JustinsParticle> lit(a_inflow_list);
      for(lit.begin(); lit.ok();) {

         // convert upbar and xpbar to upnew and xpnew 
         JustinsParticle& this_particle = lit();
         RealVect& xp = this_particle.position();
         std::array<Real,3>& up = lit().velocity();
         const RealVect& xpold = this_particle.position_old();
         const std::array<Real,3>& upold = lit().velocity_old();
         for(int n=0; n<SpaceDim; n++) xp[n] = 2.0*xp[n] - xpold[n];
         for(int n=0; n<3; n++) up[n] = 2.0*up[n] - upold[n];

         // update probes and transfer to outcast_list 
         if(xp[a_dir] < a_rightEdge) {

            const Real& wp = lit().weight();
            m_delta_MassIn_hi[a_dir] += wp;
            m_delta_MomXIn_hi[a_dir] += wp*upold[0];
            m_delta_MomYIn_hi[a_dir] += wp*upold[1];
            m_delta_MomZIn_hi[a_dir] += wp*upold[2];
         
            gbpsq = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
#ifdef RELATIVISTIC_PARTICLES
            gammap = sqrt(1.0 + gbpsq);
#else
            gammap = 1.0;
#endif
            m_delta_EnergyIn_hi[a_dir] += wp*gbpsq/(gammap + 1.0);
            a_outcast_list.transfer(lit); // this updates iterator

         }
         else { // possible to be here when using sub-orbit model for inflow particles
            a_inflow_list.remove(lit); // this updates iterator
         }

      }   

   }

}

void PicSpeciesBC::depositSurfaceChargeLo( FArrayBox&              a_sigma,
                                     const Box&                    a_bdry_box,
		                     const List<JustinsParticle>&  a_list,
                                     const int                     a_bdry_dir,
                                     const Real                    a_leftEdge )
{

   if(m_species_charge==0) { return; }
   CH_TIME("PicSpeciesBC::depositSurfaceChargeLo");
   
   // initialize the index for deposit on the boundary
   Box node_box = surroundingNodes(a_bdry_box);
   IntVect index;
   index[a_bdry_dir] = node_box.bigEnd(a_bdry_dir);

#if CH_SPACEDIM==2
   const RealVect& dX(m_mesh.getdX());
   const RealVect& Xmin(m_mesh.getXmin());
   RealVect dXp, xpnew0;
   Real slope, deltap0;
   int par_dir, index_par;
   if(a_bdry_dir==0) { par_dir = 1; }
   else if(a_bdry_dir==1) { par_dir = 0; }
   std::array<Real,3> weight_vec;
#endif
   
   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_bdry_dir] < a_leftEdge) {
         const Real& wp = lit().weight();
#if CH_SPACEDIM==1
         a_sigma(index,0) += wp*m_species_charge;
#elif CH_SPACEDIM==2
         const RealVect& xpold = lit().position_old();
         const RealVect& xpnew = lit().position();
	 dXp[0] = xpnew[0]-xpold[0];
	 dXp[1] = xpnew[1]-xpold[1];
	 slope = dXp[par_dir]/dXp[a_bdry_dir];
	 xpnew0[a_bdry_dir] = a_leftEdge;
	 xpnew0[par_dir] = xpold[par_dir] + slope*( xpnew0[a_bdry_dir] - xpold[a_bdry_dir] );
         if(m_interp_order==2) { 
	    index_par = std::floor( (xpnew0[par_dir] - 0.5*dX[par_dir] - Xmin[par_dir])/dX[par_dir] );
            deltap0 = ( xpnew0[par_dir] - ( (index_par + 1.0)*dX[par_dir] + Xmin[par_dir] ) )/dX[par_dir];
            weight_vec[0] = 0.5*std::pow((0.5 - deltap0),2);
            weight_vec[1] = 0.75 - deltap0*deltap0;
            weight_vec[2] = 0.5*std::pow((0.5 + deltap0),2);
	 }
	 else{
	    index_par = std::floor( (xpnew0[par_dir] - Xmin[par_dir])/dX[par_dir] );
            deltap0 = ( xpnew0[par_dir] - ( index_par*dX[par_dir] + Xmin[par_dir] ) )/dX[par_dir];
            weight_vec[1] = deltap0;
            weight_vec[0] = 1.0 - weight_vec[1];
	 }
	 for (int ii=0; ii<=m_interp_order; ++ii) {
	    index[par_dir] = index_par + ii;
	    a_sigma(index,0) += wp*m_species_charge*weight_vec[ii];
	 }
#endif
      }

   }
         
}

void PicSpeciesBC::depositSurfaceChargeHi( FArrayBox&              a_sigma,
                                     const Box&                    a_bdry_box,
		                     const List<JustinsParticle>&  a_list,
                                     const int                     a_bdry_dir,
                                     const Real                    a_rightEdge )
{
   if(m_species_charge==0) { return; }
   CH_TIME("PicSpeciesBC::depositSurfaceChargeHi");
   
   // initialize the index for deposit on the boundary
   Box node_box = surroundingNodes(a_bdry_box);
   IntVect index;
   index[a_bdry_dir] = node_box.smallEnd(a_bdry_dir);

#if CH_SPACEDIM==2
   const RealVect& dX(m_mesh.getdX());
   const RealVect& Xmin(m_mesh.getXmin());
   RealVect dXp, xpnew0;
   Real slope, deltap0;
   int par_dir, index_par;
   if(a_bdry_dir==0) { par_dir = 1; }
   else if(a_bdry_dir==1) { par_dir = 0; }
   std::array<Real,3> weight_vec;
#endif

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_bdry_dir] >= a_rightEdge) {
         const Real& wp = lit().weight();
#if CH_SPACEDIM==1
         a_sigma(index,0) += wp*m_species_charge;
#elif CH_SPACEDIM==2
         const RealVect& xpold = lit().position_old();
         const RealVect& xpnew = lit().position();
	 dXp[0] = xpnew[0]-xpold[0];
	 dXp[1] = xpnew[1]-xpold[1];
	 slope = dXp[par_dir]/dXp[a_bdry_dir];
	 xpnew0[a_bdry_dir] = a_rightEdge;
	 xpnew0[par_dir] = xpold[par_dir] + slope*( xpnew0[a_bdry_dir] - xpold[a_bdry_dir] );
         if(m_interp_order==2) { 
	    index_par = std::floor( (xpnew0[par_dir] - 0.5*dX[par_dir] - Xmin[par_dir])/dX[par_dir] );
            deltap0 = ( xpnew0[par_dir] - ( (index_par + 1.0)*dX[par_dir] + Xmin[par_dir] ) )/dX[par_dir];
            weight_vec[0] = 0.5*std::pow((0.5 - deltap0),2);
            weight_vec[1] = 0.75 - deltap0*deltap0;
            weight_vec[2] = 0.5*std::pow((0.5 + deltap0),2);
	 }
	 else{
	    index_par = std::floor( (xpnew0[par_dir] - Xmin[par_dir])/dX[par_dir] );
            deltap0 = ( xpnew0[par_dir] - ( index_par*dX[par_dir] + Xmin[par_dir] ) )/dX[par_dir];
            weight_vec[1] = deltap0;
            weight_vec[0] = 1.0 - weight_vec[1];
	 }
	 for (int ii=0; ii<=m_interp_order; ++ii) {
	    index[par_dir] = index_par + ii;
	    a_sigma(index,0) += wp*m_species_charge*weight_vec[ii];
	 }
#endif
      }

   }

}

void PicSpeciesBC::extractSurfaceChargeLo( FArrayBox&              a_sigma,
                                     const Box&                    a_bdry_box,
		                     const List<JustinsParticle>&  a_list,
                                     const int                     a_bdry_dir,
                                     const Real                    a_leftEdge,
                                     const bool                    a_intermediate_advance )
{

   if(m_species_charge==0) { return; }
   CH_TIME("PicSpeciesBC::extractSurfaceChargeLo");
   
   // initialize the index for deposit on the boundary
   Box node_box = surroundingNodes(a_bdry_box);
   IntVect index;
   index[a_bdry_dir] = node_box.bigEnd(a_bdry_dir);

#if CH_SPACEDIM==2
   const RealVect& dX(m_mesh.getdX());
   const RealVect& Xmin(m_mesh.getXmin());
   RealVect dXp, xpnew0;
   Real slope, deltap0;
   int par_dir, index_par;
   if(a_bdry_dir==0) { par_dir = 1; }
   else if(a_bdry_dir==1) { par_dir = 0; }
   std::array<Real,3> weight_vec;
#endif
   
   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();
      
      const RealVect& xpbar = this_particle.position();
      const RealVect& xpold = this_particle.position_old();
      RealVect this_xp = xpbar;
      if(!a_intermediate_advance) { // compute xpnew 
         for(int n=0; n<SpaceDim; n++) { this_xp[n] = 2.0*xpbar[n] - xpold[n]; }
      }

      if(this_xp[a_bdry_dir] >= a_leftEdge) {
         const Real& wp = lit().weight();
#if CH_SPACEDIM==1
         a_sigma(index,0) -= wp*m_species_charge;
#elif CH_SPACEDIM==2
	 dXp[0] = this_xp[0]-xpold[0];
	 dXp[1] = this_xp[1]-xpold[1];
	 slope = dXp[par_dir]/dXp[a_bdry_dir];
	 xpnew0[a_bdry_dir] = a_leftEdge;
	 xpnew0[par_dir] = xpold[par_dir] + slope*( xpnew0[a_bdry_dir] - xpold[a_bdry_dir] );
         if(m_interp_order==2) { 
	    index_par = std::floor( (xpnew0[par_dir] - 0.5*dX[par_dir] - Xmin[par_dir])/dX[par_dir] );
            deltap0 = ( xpnew0[par_dir] - ( (index_par + 1.0)*dX[par_dir] + Xmin[par_dir] ) )/dX[par_dir];
            weight_vec[0] = 0.5*std::pow((0.5 - deltap0),2);
            weight_vec[1] = 0.75 - deltap0*deltap0;
            weight_vec[2] = 0.5*std::pow((0.5 + deltap0),2);
	 }
	 else{
	    index_par = std::floor( (xpnew0[par_dir] - Xmin[par_dir])/dX[par_dir] );
            deltap0 = ( xpnew0[par_dir] - ( index_par*dX[par_dir] + Xmin[par_dir] ) )/dX[par_dir];
            weight_vec[1] = deltap0;
            weight_vec[0] = 1.0 - weight_vec[1];
	 }
	 for (int ii=0; ii<=m_interp_order; ++ii) {
	    index[par_dir] = index_par + ii;
	    a_sigma(index,0) -= wp*m_species_charge*weight_vec[ii];
	 }
#endif
      }

   }
         
}

void PicSpeciesBC::extractSurfaceChargeHi( FArrayBox&              a_sigma,
                                     const Box&                    a_bdry_box,
		                     const List<JustinsParticle>&  a_list,
                                     const int                     a_bdry_dir,
                                     const Real                    a_rightEdge,
                                     const bool                    a_intermediate_advance )
{
   if(m_species_charge==0) { return; }
   CH_TIME("PicSpeciesBC::extractSurfaceChargeHi");
   
   // initialize the index for deposit on the boundary
   Box node_box = surroundingNodes(a_bdry_box);
   IntVect index;
   index[a_bdry_dir] = node_box.smallEnd(a_bdry_dir);

#if CH_SPACEDIM==2
   const RealVect& dX(m_mesh.getdX());
   const RealVect& Xmin(m_mesh.getXmin());
   RealVect dXp, xpnew0;
   Real slope, deltap0;
   int par_dir, index_par;
   if(a_bdry_dir==0) { par_dir = 1; }
   else if(a_bdry_dir==1) { par_dir = 0; }
   std::array<Real,3> weight_vec;
#endif

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok(); ++lit) {
      
      JustinsParticle& this_particle = lit();

      const RealVect& xpbar = this_particle.position();
      const RealVect& xpold = this_particle.position_old();
      RealVect this_xp = xpbar;
      if(!a_intermediate_advance) { // compute xpnew 
         for(int n=0; n<SpaceDim; n++) { this_xp[n] = 2.0*xpbar[n] - xpold[n]; }
      }
      
      if(this_xp[a_bdry_dir] < a_rightEdge) {
         const Real& wp = lit().weight();
#if CH_SPACEDIM==1
         a_sigma(index,0) -= wp*m_species_charge;
#elif CH_SPACEDIM==2
	 dXp[0] = this_xp[0]-xpold[0];
	 dXp[1] = this_xp[1]-xpold[1];
	 slope = dXp[par_dir]/dXp[a_bdry_dir];
	 xpnew0[a_bdry_dir] = a_rightEdge;
	 xpnew0[par_dir] = xpold[par_dir] + slope*( xpnew0[a_bdry_dir] - xpold[a_bdry_dir] );
         if(m_interp_order==2) { 
	    index_par = std::floor( (xpnew0[par_dir] - 0.5*dX[par_dir] - Xmin[par_dir])/dX[par_dir] );
            deltap0 = ( xpnew0[par_dir] - ( (index_par + 1.0)*dX[par_dir] + Xmin[par_dir] ) )/dX[par_dir];
            weight_vec[0] = 0.5*std::pow((0.5 - deltap0),2);
            weight_vec[1] = 0.75 - deltap0*deltap0;
            weight_vec[2] = 0.5*std::pow((0.5 + deltap0),2);
	 }
	 else{
	    index_par = std::floor( (xpnew0[par_dir] - Xmin[par_dir])/dX[par_dir] );
            deltap0 = ( xpnew0[par_dir] - ( index_par*dX[par_dir] + Xmin[par_dir] ) )/dX[par_dir];
            weight_vec[1] = deltap0;
            weight_vec[0] = 1.0 - weight_vec[1];
	 }
	 for (int ii=0; ii<=m_interp_order; ++ii) {
	    index[par_dir] = index_par + ii;
	    a_sigma(index,0) -= wp*m_species_charge*weight_vec[ii];
	 }
#endif
      }

   }

}

void PicSpeciesBC::remove_Lo( List<JustinsParticle>&  a_list,
                        const int                     a_dir,
                        const Real                    a_leftEdge )
{

   CH_TIME("PicSpeciesBC::remove_Lo");
   
   Real gbpsq, gammap;

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok();) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_dir] < a_leftEdge) {

         const Real& wp = lit().weight();
         const std::array<Real,3>& up = lit().velocity();
         m_delta_MassOut_lo[a_dir] += wp;
         m_delta_MomXOut_lo[a_dir] += wp*up[0];
         m_delta_MomYOut_lo[a_dir] += wp*up[1];
         m_delta_MomZOut_lo[a_dir] += wp*up[2];
         
         gbpsq = up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
#ifdef RELATIVISTIC_PARTICLES
         gammap = sqrt(1.0 + gbpsq);
#else
         gammap = 1.0;
#endif
         m_delta_EnergyOut_lo[a_dir] += wp*gbpsq/(gammap + 1.0);
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

   Real gbpsq, gammap;

   ListIterator<JustinsParticle> lit(a_list);
   for(lit.begin(); lit.ok();) {
      
      JustinsParticle& this_particle = lit();
      RealVect& this_x = this_particle.position();
      if(this_x[a_dir] >= a_rightEdge) {

         const Real& wp = lit().weight();
         const std::array<Real,3>& up = lit().velocity();
         m_delta_MassOut_hi[a_dir] += wp;
         m_delta_MomXOut_hi[a_dir] += wp*up[0];
         m_delta_MomYOut_hi[a_dir] += wp*up[1];
         m_delta_MomZOut_hi[a_dir] += wp*up[2];

         gbpsq = up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
#ifdef RELATIVISTIC_PARTICLES
         gammap = sqrt(1.0 + gbpsq);
#else
         gammap = 1.0;
#endif
         m_delta_EnergyOut_hi[a_dir] += wp*gbpsq/(gammap + 1.0);
         a_list.remove(lit); // this updates iterator

      }
      else {
         ++lit;
      }

   }

}


#include "NamespaceFooter.H"
