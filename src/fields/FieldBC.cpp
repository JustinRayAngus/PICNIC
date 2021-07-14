#include "FieldBC.H"
//#include "CodimBC.H"
//#include "FieldBCUtils.H"

#include "NamespaceHeader.H"


FieldBC::FieldBC( const std::string&  a_variable_name,
                  const DomainGrid&   a_mesh,
                  const int           a_verbosity )
   : m_variable_name(a_variable_name),
     m_mesh(a_mesh),
     //m_insulator_conductor_bc(false),
     m_verbosity(a_verbosity)
{
   m_bc_type.resize(NUM_BOUNDARIES);
   m_bdry_name.resize(NUM_BOUNDARIES);

   setNames();

   // parse to set m_bc_type (includes all boundaries)
   string pp_prefix = "BC.em_field." + m_variable_name;
   ParmParse pp(pp_prefix.c_str());
   parseParameters(pp);
   
   // set m_domain_bc_type and m_domain_bdry_names (does not include periodic directions)
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   setDomainTypesAndNames( bdry_layout );

   if(m_verbosity) printParameters();

}

FieldBC::~FieldBC()
{
   m_bdry_name.resize(0);
   m_bc_type.resize(0);
}

void
FieldBC::setNames()
{
   m_bdry_name[DIR0_LOWER] = "dir0_lower";
   m_bdry_name[DIR0_UPPER] = "dir0_upper";
#if CH_SPACEDIM>=2
   m_bdry_name[DIR1_LOWER] = "dir1_lower";
   m_bdry_name[DIR1_UPPER] = "dir1_upper";
#endif
#if CH_SPACEDIM==3
   m_bdry_name[DIR2_LOWER] = "dir2_lower";
   m_bdry_name[DIR2_UPPER] = "dir2_upper";
#endif
}

std::string
FieldBC::getBCType( const int  a_dir,
                    const int  a_side )
{
   std::string bc_type;

   if ( a_dir == 0 ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[DIR0_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[DIR0_UPPER];
      }
      else {
         MayDay::Error("FieldBC::getBCType(): Invalid side argument");
      }
   }
#if CH_SPACEDIM>=2
   else if ( a_dir == 1 ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[DIR1_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[DIR1_UPPER];
      }
      else {
         MayDay::Error("FieldBC::getBCType(): Invalid side argument");
      }
   }
#endif
#if CH_SPACEDIM==3
   else if ( a_dir == 2 ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[DIR2_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[DIR2_UPPER];
      }
      else {
         MayDay::Error("FieldBC::getBCType(): Invalid side argument");
      }
   }
#endif
   else {
      cout << " a_dir = " << a_dir << " a_side = " << a_dir << endl;
      MayDay::Error("FieldBC::getBCType(): Invalid direction argument");
   }

   return bc_type;
}

std::string
FieldBC::getBdryName( const int  a_dir,
                      const int  a_side )
{
   std::string bc_name;

   if ( a_dir == 0 ) {
      if ( a_side == 0 ) {
         bc_name = m_bdry_name[DIR0_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_name = m_bdry_name[DIR0_UPPER];
      }
      else {
         MayDay::Error("FieldBC::getBdryName(): Invalid side argument");
      }
   }
#if CH_SPACEDIM>=2
   else if ( a_dir == 1 ) {
      if ( a_side == 0 ) {
         bc_name = m_bdry_name[DIR1_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_name = m_bdry_name[DIR1_UPPER];
      }
      else {
         MayDay::Error("FieldBC::getBdryName(): Invalid side argument");
      }
   }
#endif
#if CH_SPACEDIM==3
   else if ( a_dir == 2 ) {
      if ( a_side == 0 ) {
         bc_name = m_bdry_name[DIR2_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_name = m_bdry_name[DIR2_UPPER];
      }
      else {
         MayDay::Error("FieldBC::getBdryName(): Invalid side argument");
      }
   }
#endif
   else {
      cout << " a_dir = " << a_dir << " a_side = " << a_dir << endl;
      MayDay::Error("FieldBC::getBdryName(): Invalid direction argument");
   }

   return bc_name;
}

void FieldBC::parseParameters( ParmParse&  a_pp )
{
   
   //a_pp.query( "insulator_conductor_variable", m_insulator_conductor_bc );
   //if( m_insulator_conductor_bc ) { // define insulator/conductor object
   //   m_InsulatorConductorBC = RefCountedPtr<InsulatorConductorBC> (new InsulatorConductorBC(m_species_name,m_variable_name,m_bdry_name));
   //}
   
   for (int i(0); i<NUM_BOUNDARIES; i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string function_name;
      fpp.query( "type", m_bc_type[i] );
   }
   
}

void FieldBC::printParameters() const
{
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "FieldBC =======================================" << std::endl;
      std::cout << "- variable: "  << m_variable_name << std::endl;
      for (int i(0); i<m_domain_bc_type.size(); i++) {
         std::cout << "  " << m_domain_bdry_name[i] << ": " << std::endl;
         std::cout << "  " << m_domain_bc_type[i] << ": " << std::endl;
      }
      std::cout << "===============================================" << std::endl;
   }
}

void FieldBC::applyCellBC( LevelData<FArrayBox>&   a_dst,
                     const Real                    a_time )
{
   CH_TIME("FieldBC::applyCellBC()");

   /*
   // apply BCs to cell variable
   LevelData<FArrayBox>& cell_var( a_species_phys.cell_var( m_variable_name ) );
   FluidBCUtils::setCellBC( cell_var,
                            m_bdry_layout,
                            m_domain_bc_type );
   
   if(m_insulator_conductor_bc) { // apply InsulatorConductorBC if defined
      m_InsulatorConductorBC->applyBC( cell_var, m_bdry_layout, m_domain_bc_type, a_time );
   }
   
//#define NEW_CODIM_CORNER_BC_METHOD
#ifdef NEW_CODIM_CORNER_BC_METHOD
   // interpolate to fill physical corner boundaries
   // JRA, something wrong here for multi-block? (see two_field_neutrals)
   // JRA, Also, code is actually slower in 3D. While time spent here goes
   //      down, for some reason time spent doing exchange in convertToPhysical
   //      goes up by a lot. Don't understand this.
   CodimBC::setCodimCornerValues( cell_var, coord_sys );
#else
   // fills both physical and internal corners with extrapolation
   CodimBC::setCodimBoundaryValues( cell_var, coord_sys );
   
   // copy internal corners with non-corner values via exchange
   geometry.fillCorners(cell_var, cell_var.ghostVect(), SpaceDim); 
#endif
  */
 
}

void FieldBC::applyFluxBC( LevelData<FluxBox>&     a_dst,
                     const Real                    a_time )
{
   CH_TIME("FieldBC::applyFluxBC()");
  
   /*
   // apply BCs to a_dst
   FluidBCUtils::setFluxBC( a_dst,
                            m_bdry_layout,
                            m_domain_bc_type,
                            geometry );
   */
}

void FieldBC::applyEdgeBC( LevelData<EdgeDataBox>&  a_dst,
                     const Real                     a_time )
{
   CH_TIME("FieldBC::applyEdgeBC()");
   
   /*
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   if(dst_ghost_vect > IntVect::Zero) {
      FluidBCUtils::setEdgeBC( a_dst,
                               m_bdry_layout,
                               m_domain_bc_type );
   }
   
   if(m_insulator_conductor_bc) {
      const MagGeom& geometry( a_species_phys.configurationSpaceGeometry() );
      m_InsulatorConductorBC->applyEdgeBC( a_dst, m_bdry_layout, m_domain_bc_type, geometry );
   }
   */

}
/*
void FieldBC::setFluxBC(const FluidSpecies&        a_species_phys,
                           LevelData<FluxBox>&        a_dst,
                           const LevelData<FluxBox>&  a_src,
                           const Real                 a_time )
{
   CH_TIME("FieldBC::setFluxBC()");
      
   // check that src and dst ghost vects are what they need to be
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   const IntVect& src_ghost_vect( a_src.ghostVect() );
   CH_assert(dst_ghost_vect==src_ghost_vect || src_ghost_vect==IntVect::Zero);   

   // copy a_src to a_dst in boundary region
   FluidBCUtils::setFluxBC( a_dst,
                            m_bdry_layout,
                            a_src ) ;
   
}

void FieldBC::setEdgeBC(const FluidSpecies&            a_species_phys,
                           LevelData<EdgeDataBox>&        a_dst,
                           const LevelData<EdgeDataBox>&  a_src,
                           const Real                     a_time )
{
   CH_TIME("FieldBC::setEdgeBC()");
   
   // check that src and dst ghost vects are what they need to be
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   const IntVect& src_ghost_vect( a_src.ghostVect() );
   CH_assert(dst_ghost_vect==src_ghost_vect || src_ghost_vect==IntVect::Zero);   
   
   // copy a_src to a_dst in boundary region
   FluidBCUtils::setEdgeBC( a_dst,
                            m_bdry_layout,
                            a_src );
   
}

void FieldBC::setOnAxisCurlBC( LevelData<EdgeDataBox>&  a_curl_covar,
                            const LevelData<FArrayBox>&    a_By_phys,
                            const FluidSpecies&            a_species_phys )
{
   CH_TIME("FieldBC::setOnAxisCurlBC()");
   
   // Compute the appropriate BC on axis for the curl of a function
   // that is odd wrt the axis for axisymmetric systems where the Jacobian is zero.
   //
   // Assumed a_curl is covariant curl(a_By_phys) on cell edges   
   
   const MagGeom& geometry( a_species_phys.configurationSpaceGeometry() );
   
   FluidBCUtils::setOnAxisCurlBC( a_curl_covar,
                                  a_By_phys,
                                  geometry,
                                  m_bdry_layout,
                                  m_domain_bc_type );

}

void FieldBC::setInsulatorBC( const FluidSpecies&    a_species_comp,
                                 LevelData<FArrayBox>&  a_dst,
                           const LevelData<FArrayBox>&  a_src,
                           const Real                   a_time )
{
   CH_TIME("FieldBC::setInsulatorBC()");
  
   // set a_dst equal to a_src in boundary region where an insulator exists
   
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   
   if(m_insulator_conductor_bc) {
      m_InsulatorConductorBC->setInsulatorBC( a_dst, a_src, m_bdry_layout, m_domain_bc_type, geometry, a_time );
   }
   
}
*/
      
void FieldBC::setDomainTypesAndNames( const BoundaryBoxLayoutPtrVect&  a_bdry_layout )
{
   m_domain_bc_type.resize(a_bdry_layout.size());
   m_domain_bdry_name.resize(a_bdry_layout.size());
   for (int i(0); i<a_bdry_layout.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(a_bdry_layout[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      m_domain_bc_type[i] = getBCType(dir, side);
      m_domain_bdry_name[i] = getBdryName(dir, side);
   }
}

#include "NamespaceFooter.H"
