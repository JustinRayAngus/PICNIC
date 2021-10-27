#include "FieldBC.H"
#include "CodimBC.H"
#include "FieldBCUtils.H"

#include "NamespaceHeader.H"


FieldBC::FieldBC( const DomainGrid&  a_mesh,
                  const int          a_verbosity )
   : m_mesh(a_mesh),
     m_verbosity(a_verbosity)
{
   // parse the input file
   string pp_prefix = "BC.em_fields";
   ParmParse pp(pp_prefix.c_str());
   parseParameters(pp);
   
   if(m_verbosity) printParameters();

}

FieldBC::~FieldBC()
{
   m_bdry_name.resize(0);
   m_bc_type.resize(0);
}

void FieldBC::parseParameters( ParmParse&  a_pp )
{
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   m_bc_type.resize(bdry_layout.size());
   m_bdry_name.resize(bdry_layout.size());
   
   for (int i(0); i<bdry_layout.size(); i++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[i]) );
      m_bdry_name[i] = this_bdry_layout.name();
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      if( fpp.contains("type") ) {
         fpp.get( "type", m_bc_type[i] );
         CH_assert( m_bc_type[i] == "symmetry" ||    
                    m_bc_type[i] == "odd" ||
                    m_bc_type[i] == "even" ||
                    m_bc_type[i] == "zero" ||
                    m_bc_type[i] == "conductor" ||
                    m_bc_type[i] == "extrapolate" || 
                    m_bc_type[i] == "extrapolate_zeroBv" || 
                    m_bc_type[i] == "insulator_conductor" );
         if(m_bc_type[i]=="insulator_conductor") {
            m_InsulatorBC = RefCountedPtr<InsulatorBC> (new InsulatorBC(i,m_bdry_name[i],m_mesh));
         }
      }
   }
   
}

void FieldBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "FieldBC =======================================" << std::endl;
      for (int i(0); i<m_bc_type.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << m_bc_type[i] << ": " << std::endl;
      }
      std::cout << "===============================================" << std::endl << std::endl;
   }
}

void FieldBC::applyCellBC( LevelData<FArrayBox>&   a_dst,
                     const Real                    a_time )
{
   CH_TIME("FieldBC::applyCellBC()");

   // apply BCs to cell variable (virtual B)

   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   FieldBCUtils::setCellBC( a_dst,
                            bdry_layout,
                            m_bc_type,
                            m_InsulatorBC,
                            a_time );
   
   CodimBC::setCodimCornerValues( a_dst, m_mesh );
    
}

void FieldBC::applyFluxBC( LevelData<FluxBox>&     a_dst,
                     const Real                    a_time )
{
   CH_TIME("FieldBC::applyFluxBC()");
   
   // apply BCs to face variable (in-plane B)
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   FieldBCUtils::setFluxBC( a_dst,
                            bdry_layout,
                            m_bc_type,
                            m_InsulatorBC,
                            a_time );
   
}

void FieldBC::applyEdgeBC( LevelData<EdgeDataBox>&  a_dst,
                     const Real                     a_time )
{
   CH_TIME("FieldBC::applyEdgeBC()");
   
   // apply BCs to face variable (in-plane E)
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   FieldBCUtils::setEdgeBC( a_dst,
                            bdry_layout,
                            m_bc_type,
                            m_InsulatorBC,
                            a_time );
   
}

void FieldBC::applyNodeBC( LevelData<NodeFArrayBox>&  a_dst,
                     const Real                       a_time )
{
   CH_TIME("FieldBC::applyNodeBC()");

   // apply BCs to node variable (virtual E)
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   FieldBCUtils::setNodeBC( a_dst,
                            bdry_layout,
                            m_bc_type,
                            m_InsulatorBC,
                            a_time );

}

void FieldBC::setFluxBC( LevelData<FluxBox>&  a_dst,
                   const LevelData<FluxBox>&  a_src,
                   const Real                 a_time )
{
   CH_TIME("FieldBC::setFluxBC()");
      
   // check that src and dst ghost vects are what they need to be
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   const IntVect& src_ghost_vect( a_src.ghostVect() );
   CH_assert(dst_ghost_vect==src_ghost_vect || src_ghost_vect==IntVect::Zero);   
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();

   // copy a_src to a_dst in boundary region
   FieldBCUtils::setFluxBC( a_dst,
                            bdry_layout,
                            a_src ) ;
   
}

void FieldBC::setEdgeBC( LevelData<EdgeDataBox>&  a_dst,
                   const LevelData<EdgeDataBox>&  a_src,
                   const Real                     a_time )
{
   CH_TIME("FieldBC::setEdgeBC()");
   
   // check that src and dst ghost vects are what they need to be
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   const IntVect& src_ghost_vect( a_src.ghostVect() );
   CH_assert(dst_ghost_vect==src_ghost_vect || src_ghost_vect==IntVect::Zero);   
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   
   // copy a_src to a_dst in boundary region
   FieldBCUtils::setEdgeBC( a_dst,
                            bdry_layout,
                            a_src );
   
}

#include "NamespaceFooter.H"
