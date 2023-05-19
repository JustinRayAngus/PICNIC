#include "FieldBC.H"
#include "CodimBC.H"
#include "FieldBCUtils.H"
#include "FieldsF_F.H"

#include "NamespaceHeader.H"


FieldBC::FieldBC( const DomainGrid&  a_mesh,
                  const int          a_verbosity )
   : m_conservative_wall(false),
     m_mesh(a_mesh),
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
   m_InsulatorBC.resize(0);
}

void FieldBC::parseParameters( ParmParse&  a_pp )
{
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   m_bc_type.resize(bdry_layout.size());
   m_bdry_name.resize(bdry_layout.size());
   m_InsulatorBC.resize(bdry_layout.size());

   a_pp.query("conservative_wall",m_conservative_wall);

   for (int b(0); b<bdry_layout.size(); b++) {
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      m_bdry_name[b] = this_bdry_layout.name();
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[b];
      ParmParse fpp( prefix.c_str() );
      if( fpp.contains("type") ) {
         fpp.get( "type", m_bc_type[b] );
         CH_assert( m_bc_type[b] == "symmetry" ||    
                    m_bc_type[b] == "odd" ||
                    m_bc_type[b] == "even" ||
                    m_bc_type[b] == "zero" ||
                    m_bc_type[b] == "conductor" ||
                    m_bc_type[b] == "axis" ||
                    m_bc_type[b] == "extrapolate" || 
                    m_bc_type[b] == "extrapolate_zeroBv" || 
                    m_bc_type[b] == "insulator_conductor" );
         if(m_bc_type[b]=="insulator_conductor") {
            m_InsulatorBC[b] = RefCountedPtr<InsulatorBC> (new InsulatorBC(b,m_bdry_name[b],m_mesh));
         }
         if(m_bc_type[b]=="axis" && !m_mesh.axisymmetric()) m_bc_type[b] = "symmetry";
         if(m_bc_type[b]=="symmetry" && m_mesh.axisymmetric()) {
            int bdry_side(this_bdry_layout.side());
            const RealVect& Xmin(m_mesh.getXmin());
            if(bdry_side==0 && Xmin[0]==0.0) m_bc_type[b] = "axis";
         }
      }
   }
   
}

void FieldBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "FieldBC =======================================" << std::endl;
      std::cout << "conservative wall = " << m_conservative_wall << std::endl;
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
   FieldBCUtils::setCellBC( a_dst,
                            m_mesh,
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
                            m_conservative_wall,
                            a_time );
   
}

void FieldBC::applyNodeBC( LevelData<NodeFArrayBox>&  a_dst,
                     const Real                       a_time )
{
   CH_TIME("FieldBC::applyNodeBC()");

   // apply BCs to node variable (virtual E)
   FieldBCUtils::setNodeBC( a_dst,
                            m_mesh,
                            m_bc_type,
                            m_InsulatorBC,
                            m_conservative_wall,
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

void FieldBC::applyOnAxisCurlBC( LevelData<EdgeDataBox>&  a_curlB,
                           const LevelData<FArrayBox>&    a_B )
{
   CH_TIME("FieldBC::applyOnAxisCurlBC()");
   if(!m_mesh.axisymmetric()) return;      
         
   // set boundry value for curlB assuming B~r ==> 1/r*d(rB)/dr = 4*B(1/2)/dX
   
   const RealVect& dX(m_mesh.getdX());      

   const BoundaryBoxLayoutPtrVect& all_bdry_layouts = m_mesh.getBoundaryLayout();
   for(int b(0); b<all_bdry_layouts.size(); b++) {
      const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
      const int bdry_dir = bdry_layout.dir();
      const int bdry_side(bdry_layout.side());
      const std::string this_bc_type (m_bc_type[b]);
      if(this_bc_type!="axis") continue;      
         
      int curl_dir = 0;
      if(bdry_dir==0) curl_dir = 1;

      const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const FArrayBox& this_B( a_B[bdry_layout.dataIndex(dit)] );
         FArrayBox& this_curlB( a_curlB[bdry_layout.dataIndex(dit)][bdry_dir] );
         
         // create fill box right on axis
         const Box& bdry_box( bdry_grids[dit] );
         Box edge_box = surroundingNodes(bdry_box);
         edge_box.enclosedCells(curl_dir);
               
         // collapse edge_box to 1 cell thick in bdry_dir direction                  
         if(bdry_side==0) edge_box.setSmall(bdry_dir,edge_box.bigEnd(bdry_dir));
         if(bdry_side==1) edge_box.setBig(bdry_dir,edge_box.smallEnd(bdry_dir));

         // set boundry value for curlB
         IntVect ib, shift;
         shift = IntVect::Zero;
         shift[bdry_dir] = -bdry_side;
         Real bdry_val, B0;
         BoxIterator bit(edge_box);
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            B0 = this_B.get(ib+shift,0);
            bdry_val = 4.0*B0/dX[bdry_dir];
            this_curlB.set(ib,0,bdry_val);
         }

      }

   }

}

void FieldBC::applyOnAxisCurlBC( LevelData<NodeFArrayBox>&  a_curlB,
                           const LevelData<FArrayBox>&      a_B )
{
   CH_TIME("FieldBC::applyOnAxisCurlBC()");
   if(!m_mesh.axisymmetric()) return;      
   if(a_curlB.nComp()==1) return;   
      
   // set boundry value for curlBy assuming By~r ==> d(rBy)/dr/r = 2*By(1/2)/r(1/2)
   
   const RealVect& dX(m_mesh.getdX());      

   const BoundaryBoxLayoutPtrVect& all_bdry_layouts = m_mesh.getBoundaryLayout();
   for(int b(0); b<all_bdry_layouts.size(); b++) {
      const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
      const int bdry_dir = bdry_layout.dir();
      const int bdry_side(bdry_layout.side());
      const std::string this_bc_type (m_bc_type[b]);
      if(this_bc_type!="axis") continue;      
         
      const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const FArrayBox& this_B( a_B[bdry_layout.dataIndex(dit)] );
         FArrayBox& this_curlB( a_curlB[bdry_layout.dataIndex(dit)].getFab() );
         
         // create fill box right on axis
         const Box& bdry_box( bdry_grids[dit] );
         Box node_box = surroundingNodes(bdry_box);
               
         // collapse node_box to 1 cell thick in bdry_dir direction                  
         if(bdry_side==0) node_box.setSmall(bdry_dir,node_box.bigEnd(bdry_dir));
         if(bdry_side==1) node_box.setBig(bdry_dir,node_box.smallEnd(bdry_dir));

         // set boundry value for curlB
         IntVect ib, shift;
         shift = IntVect::Zero;
         shift[bdry_dir] = -bdry_side;
         Real bdry_val, B0;
         const int this_comp = 1; // Z-dir is comp 1 for 1D
         BoxIterator bit(node_box);
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            B0 = this_B.get(ib+shift,0);
            bdry_val = 4.0*B0/dX[bdry_dir];
            this_curlB.set(ib,this_comp,bdry_val);
         }

      }

   }

}

void FieldBC::applyOnAxisDivBC( LevelData<NodeFArrayBox>&  a_divE,
                          const LevelData<EdgeDataBox>&    a_E )
{
   CH_TIME("FieldBC::applyOnAxisDivBC()");
   if(!m_mesh.axisymmetric()) return;      
      
   // set boundry value for d(rEr)/dr/r assuming Er~r ==> d(rEr)/dr/r = 2*Er(1/2)/r(1/2)
   
   const RealVect& dX(m_mesh.getdX());      

   const BoundaryBoxLayoutPtrVect& all_bdry_layouts = m_mesh.getBoundaryLayout();
   for(int b(0); b<all_bdry_layouts.size(); b++) {
      const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
      const int bdry_dir = bdry_layout.dir();
      const int bdry_side(bdry_layout.side());
      const std::string this_bc_type (m_bc_type[b]);
      if(this_bc_type!="axis") continue;      
         
      const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const FArrayBox& this_Edir( a_E[bdry_layout.dataIndex(dit)][bdry_dir] );
         FArrayBox& this_divE( a_divE[bdry_layout.dataIndex(dit)].getFab() );
         
         // create fill box right on axis
         const Box& bdry_box( bdry_grids[dit] );
         Box node_box = surroundingNodes(bdry_box);
               
         // collapse node_box to 1 cell thick in bdry_dir direction                  
         if(bdry_side==0) node_box.setSmall(bdry_dir,node_box.bigEnd(bdry_dir));
         if(bdry_side==1) node_box.setBig(bdry_dir,node_box.smallEnd(bdry_dir));

         // set boundry value for d(rEr)/dr/r
         IntVect ib, shift;
         shift = IntVect::Zero;
         shift[bdry_dir] = -bdry_side;
         Real bdry_val, E0;
         BoxIterator bit(node_box);
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            E0 = this_Edir.get(ib+shift,0);
            bdry_val = 4.0*E0/dX[bdry_dir];
            this_divE.set(ib,0,bdry_val);
         }

      }

   }

}
      
void FieldBC::computeIntSdA( RealVect&                  a_intSdA_lo, 
                             RealVect&                  a_intSdA_hi,
                       const LevelData<EdgeDataBox>&    a_E, 
                       const LevelData<FluxBox>&        a_B,
                       const LevelData<NodeFArrayBox>&  a_Ev,
                       const LevelData<FArrayBox>&      a_Bv )
{
   CH_TIME("FieldBC::computIntSdA()");

   const BoundaryBoxLayoutPtrVect& all_bdry_layouts = m_mesh.getBoundaryLayout();

   const RealVect& dA(m_mesh.getMappedFaceArea());      
   const LevelData<NodeFArrayBox>& masked_Jnc = m_mesh.getMaskedJnc();
#if CH_SPACEDIM==1 
   const LevelData<FArrayBox>& Xcc = m_mesh.getXcc();
   const LevelData<NodeFArrayBox>& Xnc = m_mesh.getXnc();
#elif CH_SPACEDIM==2 
   const LevelData<EdgeDataBox>& masked_Jec = m_mesh.getMaskedJec();
#endif
      
   for(int b(0); b<all_bdry_layouts.size(); b++) {
      const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
      const std::string this_bc_type (m_bc_type[b]);
         
      const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const Box& bdry_box( bdry_grids[dit] );
         const int bdry_dir = bdry_layout.dir();
         const int bdry_side(bdry_layout.side());
         const Real bdry_dA = dA[bdry_dir]*2.0; // 2.0 needed because masked Jacobians have 
                                                // factor of 0.5 on surface. For 1D, can use
                                                // normal Jacobians, but for D>1 need masked
                                                // to get shared edges/corners correct
  
         const FArrayBox& this_Bv( a_Bv[bdry_layout.dataIndex(dit)] );
         const FArrayBox& this_Ev( a_Ev[bdry_layout.dataIndex(dit)].getFab() );
         const FArrayBox& this_Jnc( masked_Jnc[bdry_layout.dataIndex(dit)].getFab() );         
     
         // get node box
         Box node_box = surroundingNodes(bdry_box);
               
         // collapse node_box to 1 cell thick in bdry_dir direction                  
         if(bdry_side==0) node_box.setSmall(bdry_dir,node_box.bigEnd(bdry_dir));
         if(bdry_side==1) node_box.setBig(bdry_dir,node_box.smallEnd(bdry_dir));

#if CH_SPACEDIM==1
         const FArrayBox& this_Xcc( Xcc[bdry_layout.dataIndex(dit)] );         
         const FArrayBox& this_Xnc( Xnc[bdry_layout.dataIndex(dit)].getFab() );         
            
         FArrayBox EyBz(node_box,1);
         FArrayBox EzBy(node_box,1);

         // Ey*Bz
         SpaceUtils::interpolateStag(EyBz,node_box,0,this_Bv,1,bdry_dir);
         EyBz.mult(this_Ev,node_box,0,0,1);
         
         // Ez*By
         if(m_mesh.axisymmetric()) {
            if(this_bc_type=="axis") EzBy.setVal(0.0);
            else {
               FArrayBox rBy(this_Bv.box(),1);
               rBy.copy(this_Bv,0,0,1);           
               rBy.mult(this_Xcc,this_Bv.box(),0,0,1);           
               SpaceUtils::interpolateStag(EzBy,node_box,0,rBy,0,bdry_dir);
               EzBy.mult(this_Ev,node_box,1,0,1);
               EzBy.divide(this_Xnc,node_box,0,0,1);
            }
         }
         else {
            SpaceUtils::interpolateStag(EzBy,node_box,0,this_Bv,0,bdry_dir);
            EzBy.mult(this_Ev,node_box,1,0,1);
         }

         //if(bdry_side==0) a_intSdA_lo[bdry_dir] = -(EyBz.sum(0) - EzBy.sum(0))*bdry_dA; 
         //if(bdry_side==1) a_intSdA_hi[bdry_dir] =  (EyBz.sum(0) - EzBy.sum(0))*bdry_dA; 
         if(bdry_side==0) {
            a_intSdA_lo[bdry_dir] = -( EyBz.dotProduct(this_Jnc,node_box) 
                                  -    EzBy.dotProduct(this_Jnc,node_box) )*bdry_dA;
         } 
         if(bdry_side==1) {
            a_intSdA_hi[bdry_dir] =  ( EyBz.dotProduct(this_Jnc,node_box) 
                                  -    EzBy.dotProduct(this_Jnc,node_box) )*bdry_dA; 
         }

#elif CH_SPACEDIM==2 
            
         if(bdry_dir==0) {

            const FArrayBox& this_Ey( a_E[bdry_layout.dataIndex(dit)][1] );
            const FArrayBox& this_By( a_B[bdry_layout.dataIndex(dit)][1] );
            const FArrayBox& this_J( masked_Jec[bdry_layout.dataIndex(dit)][1] );

            Box edge_box = enclosedCells(node_box,1);
            FArrayBox EyBz(edge_box,1);
            FArrayBox EzBy(node_box,1);
               
            // Ey*Bz
            SpaceUtils::interpolateStag(EyBz,edge_box,0,this_Bv,0,bdry_dir);
            EyBz.mult(this_Ey,edge_box,0,0,1);
               
            // Ez*By
            SpaceUtils::interpolateStag(EzBy,node_box,0,this_By,0,bdry_dir);
            EzBy.mult(this_Ev,node_box,0,0,1);
         
            if(bdry_side==0) {
               a_intSdA_lo[bdry_dir] = -( EyBz.dotProduct(this_J,edge_box) 
                                     -    EzBy.dotProduct(this_Jnc,node_box) )*bdry_dA; 
            }
            if(bdry_side==1) {
               a_intSdA_hi[bdry_dir] =  ( EyBz.dotProduct(this_J,edge_box) 
                                     -    EzBy.dotProduct(this_Jnc,node_box) )*bdry_dA;
            } 

         } 
            
         else if(bdry_dir==1) {
            
            const FArrayBox& this_Ex( a_E[bdry_layout.dataIndex(dit)][0] );
            const FArrayBox& this_Bx( a_B[bdry_layout.dataIndex(dit)][0] );
            const FArrayBox& this_J( masked_Jec[bdry_layout.dataIndex(dit)][0] );
               
            Box edge_box = enclosedCells(node_box,0);
            FArrayBox ExBz(edge_box,1);
            FArrayBox EzBx(node_box,1);
               
            // Ex*Bz
            SpaceUtils::interpolateStag(ExBz,edge_box,0,this_Bv,0,bdry_dir);
            ExBz.mult(this_Ex,edge_box,0,0,1);
               
            // Ez*Bx
            SpaceUtils::interpolateStag(EzBx,node_box,0,this_Bx,0,bdry_dir);
            EzBx.mult(this_Ev,node_box,0,0,1);
         
            if(bdry_side==0) {
               a_intSdA_lo[bdry_dir] = -( EzBx.dotProduct(this_Jnc,node_box) 
                                     -    ExBz.dotProduct(this_J,edge_box) )*bdry_dA; 
            }
            if(bdry_side==1) {
               a_intSdA_hi[bdry_dir] =  ( EzBx.dotProduct(this_Jnc,node_box) 
                                     -    ExBz.dotProduct(this_J,edge_box) )*bdry_dA; 
            }

         }
#endif 

      }
 
   }

}

void FieldBC::applyToJ( LevelData<EdgeDataBox>&    a_J_inPlane,
                        LevelData<NodeFArrayBox>&  a_J_virtual ) const
{
   CH_TIME("FieldBC::applyToJ()");
   
   FieldBCUtils::applyToJ_PIC( a_J_inPlane,
                               a_J_virtual, 
                               m_mesh,
                               m_bc_type );

}

#if CH_SPACEDIM==1
void FieldBC::prepJforBC( LevelData<EdgeDataBox>&    a_Jx,
                          LevelData<NodeFArrayBox>&  a_Jv,
                    const LevelData<EdgeDataBox>&    a_Jx0,
                    const LevelData<NodeFArrayBox>&  a_Jv0,
                    const LevelData<EdgeDataBox>&    a_Ex0,
                    const LevelData<NodeFArrayBox>&  a_Ev0,
                    const LevelData<EdgeDataBox>&    a_Ex,
                    const LevelData<NodeFArrayBox>&  a_Ev,
                    const LevelData<EdgeDataBox>&    a_sigma_xx,
                    const LevelData<EdgeDataBox>&    a_sigma_xy,
                    const LevelData<EdgeDataBox>&    a_sigma_xz,
                    const LevelData<NodeFArrayBox>&  a_sigma_yx,
                    const LevelData<NodeFArrayBox>&  a_sigma_yy,
                    const LevelData<NodeFArrayBox>&  a_sigma_yz,
                    const LevelData<NodeFArrayBox>&  a_sigma_zx,
                    const LevelData<NodeFArrayBox>&  a_sigma_zy,
                    const LevelData<NodeFArrayBox>&  a_sigma_zz ) const
{
   CH_TIME("FieldBC::prepJforBC()");
   
   // WARNING.....
   //
   // NOTE THAT CARE HAS NOT BEEN TAKEN TO ENSURE CONSISTENT BEHAVIOR
   // WITH applyToJ() CALL FROM PIC SPECIES CLASS....
   // NEED TO ADDRESS THIS ISSUE AT SOME POINT
 
   // loop over non-periodic boundaries and apply BCs to J
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
         
         const FArrayBox& Ex( a_Ex[interior_dit][0] );
         const FArrayBox& Ev( a_Ev[interior_dit].getFab() );
         const FArrayBox& Jx0( a_Jx0[interior_dit][0] );
         const FArrayBox& Jv0( a_Jv0[interior_dit].getFab() );
         const FArrayBox& Ex0( a_Ex0[interior_dit][0] );
         const FArrayBox& Ev0( a_Ev0[interior_dit].getFab() );
          
         //
         // compute Jx in ghost cells from deposited mass matrices
         //  
         
         Box edge_box = surroundingNodes(bdry_box);
         edge_box.enclosedCells(bdry_dir);

         FArrayBox& Jx( a_Jx[interior_dit][0] );
         const FArrayBox& sigxx = a_sigma_xx[interior_dit][0];
         const FArrayBox& sigxy = a_sigma_xy[interior_dit][0];
         const FArrayBox& sigxz = a_sigma_xz[interior_dit][0];

         FORT_COMPUTE_JX_BDRY_FROM_MASS_MATRIX( CHF_BOX(edge_box),
                                                CHF_CONST_INT(bdry_side),
                                                CHF_CONST_FRA(sigxx),
                                                CHF_CONST_FRA(sigxy),
                                                CHF_CONST_FRA(sigxz),
                                                CHF_CONST_FRA1(Ex0,0),
                                                CHF_CONST_FRA1(Ex,0),
                                                CHF_CONST_FRA1(Ev0,0),
                                                CHF_CONST_FRA1(Ev,0),
                                                CHF_CONST_FRA1(Ev0,1),
                                                CHF_CONST_FRA1(Ev,1),
                                                CHF_CONST_FRA1(Jx0,0),
                                                CHF_FRA1(Jx,0) );

         //
         // compute Jy and Jz in ghost cells from deposited mass matrices
         //
         
         Box node_box = surroundingNodes(bdry_box);
         if(bdry_side==0) node_box.growHi(bdry_dir,-1);
         if(bdry_side==1) node_box.growLo(bdry_dir,-1);

         FArrayBox& Jv( a_Jv[interior_dit].getFab() );
         const FArrayBox& sigyx = a_sigma_yx[interior_dit].getFab();
         const FArrayBox& sigyy = a_sigma_yy[interior_dit].getFab();
         const FArrayBox& sigyz = a_sigma_yz[interior_dit].getFab();
         const FArrayBox& sigzx = a_sigma_zx[interior_dit].getFab();
         const FArrayBox& sigzy = a_sigma_zy[interior_dit].getFab();
         const FArrayBox& sigzz = a_sigma_zz[interior_dit].getFab();

         const int z_comp = 1;
         
         FORT_COMPUTE_JY_BDRY_FROM_MASS_MATRIX( CHF_BOX(node_box),
                                                CHF_CONST_INT(bdry_side),
                                                CHF_CONST_FRA(sigyx),
                                                CHF_CONST_FRA(sigyy),
                                                CHF_CONST_FRA(sigyz),
                                                CHF_CONST_FRA1(Ex0,0),
                                                CHF_CONST_FRA1(Ex,0),
                                                CHF_CONST_FRA1(Ev0,0),
                                                CHF_CONST_FRA1(Ev,0),
                                                CHF_CONST_FRA1(Ev0,z_comp),
                                                CHF_CONST_FRA1(Ev,z_comp),
                                                CHF_CONST_FRA1(Jv0,0),
                                                CHF_FRA1(Jv,0) );
         
         FORT_COMPUTE_JZ_BDRY_FROM_MASS_MATRIX( CHF_BOX(node_box),
                                                CHF_CONST_INT(bdry_side),
                                                CHF_CONST_FRA(sigzx),
                                                CHF_CONST_FRA(sigzy),
                                                CHF_CONST_FRA(sigzz),
                                                CHF_CONST_FRA1(Ex0,0),
                                                CHF_CONST_FRA1(Ex,0),
                                                CHF_CONST_FRA1(Ev0,0),
                                                CHF_CONST_FRA1(Ev,0),
                                                CHF_CONST_FRA1(Ev0,z_comp),
                                                CHF_CONST_FRA1(Ev,z_comp),
                                                CHF_CONST_FRA1(Jv0,z_comp),
                                                CHF_FRA1(Jv,z_comp) );
         
      }
      
   }

}
#endif

#include "NamespaceFooter.H"
