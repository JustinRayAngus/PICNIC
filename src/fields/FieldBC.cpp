#include "FieldBC.H"
#include "CodimBC.H"
#include "FieldBCUtils.H"
#include "FieldsF_F.H"

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
                    m_bc_type[i] == "axis" ||
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
         
   // set boundry value for curlB assuming B~r ==> 1/r*d(rB)/dr = 4*B(1/2)/X(1/2)
   
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
   const LevelData<EdgeDataBox>& masked_Jec = m_mesh.getMaskedJec();
   const LevelData<NodeFArrayBox>& masked_Jnc = m_mesh.getMaskedJnc();
      
   for(int b(0); b<all_bdry_layouts.size(); b++) {
      const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
         
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

         if(SpaceDim==1) {
            
            FArrayBox EyBz(node_box,1);
            FArrayBox EzBy(node_box,1);

            // Ey*Bz
            SpaceUtils::interpolateStag(EyBz,node_box,0,this_Bv,1,bdry_dir);
            EyBz.mult(this_Ev,node_box,0,0,1);
         
            // Ez*By
            SpaceUtils::interpolateStag(EzBy,node_box,0,this_Bv,0,bdry_dir);
            EzBy.mult(this_Ev,node_box,1,0,1);

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

         } 
            
         if(SpaceDim==2 && bdry_dir==0) {

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
            
         if(SpaceDim==2 && bdry_dir==1) {
            
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

      }
 
   }

}

void FieldBC::applyToJ( LevelData<EdgeDataBox>&    a_J_inPlane,
                        LevelData<NodeFArrayBox>&  a_J_virtual )
{
   CH_TIME("PicSpeciesBC::applyToJ()");
 
   // currently this is only used for symmetry boundaries. The contribution of the
   // current density from the mirror particles across the symmetry boundary are
   // not accounted for when J is computed (only particles inside physical domain
   // are used in J calc). For J parallel to a symmetry boundary, we multiply by
   // a factor of 2 on the physical boundary (higher-order shape functions are not
   // considered yet). For J perp to a symmetry boundary, we subtract off J in the
   // corresonding ghost cell across the symmetry boundary.
  
   // loop over non-periodic boundaries and apply BCs to J
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {

      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];

      if(this_bc=="symmetry" || this_bc=="axis") {
      
         for(DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

            const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
            const Box bdry_box( bdry_grids[dit] );
            for (int dir(0); dir<SpaceDim; dir++) {
          
               // convert cell bdry box to a edge bdry box
               Box edge_box = surroundingNodes(bdry_box);
               edge_box.enclosedCells(dir);

               // collapse edge_box to 1 cell thick in bdry_dir direction                  
               if(bdry_side==0) edge_box.setSmall(bdry_dir,edge_box.bigEnd(bdry_dir));
               if(bdry_side==1) edge_box.setBig(bdry_dir,edge_box.smallEnd(bdry_dir));
                  
               // adjust J in cells near boundary appropriately to account for
               // mirror particles across symmetry boundary
               FArrayBox& this_J( a_J_inPlane[interior_dit][dir] );
               if(dir==bdry_dir) {
                  Box dst_box = edge_box;
                  Box src_box = edge_box;
                  const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
                  for (int n=0; n<nG; n++) {
                     if(bdry_side==0) dst_box.shift(bdry_dir,1);
                     if(bdry_side==1) dst_box.shift(bdry_dir,-1);
                     this_J.minus(this_J,src_box,dst_box,0,0,1);
                     if(bdry_side==0) src_box.shift(bdry_dir,-1);
                     if(bdry_side==1) src_box.shift(bdry_dir,1);
                  }

               }
               else { // work exactly on bdry only
                  // would need to modify this if using higher-order particle shape functions
                  if(!m_mesh.axisymmetric()) this_J.mult(2.0,edge_box,0,this_J.nComp());
               }
               
            }
           
            if(SpaceDim<3 && !m_mesh.axisymmetric()) {
                       
               // convert cell bdry box to a node bdry box
               Box node_box = surroundingNodes(bdry_box);

               // collapse node_box to 1 cell thick in bdry_dir direction                  
               if(bdry_side==0) node_box.setSmall(bdry_dir,node_box.bigEnd(bdry_dir));
               if(bdry_side==1) node_box.setBig(bdry_dir,node_box.smallEnd(bdry_dir));

               // would need to modify this if using higher-order 
               // particle shape functions
               FArrayBox& this_J( a_J_virtual[interior_dit].getFab() );
               this_J.mult(2.0,node_box,0,this_J.nComp());
            }
            
         }

      }
   }
   
}

void FieldBC::prepJforBC( LevelData<EdgeDataBox>&  a_Jx,
                    const LevelData<EdgeDataBox>&  a_Jx0,
                    const LevelData<EdgeDataBox>&  a_Ex,
                    const LevelData<NodeFArrayBox>&  a_Ev,
                    const LevelData<EdgeDataBox>&  a_sigma_xx,
                    const LevelData<EdgeDataBox>&  a_sigma_xy,
                    const LevelData<EdgeDataBox>&  a_sigma_xz,
                    const Real                     a_volume_scale )
{
   CH_TIME("FieldBC::prepJforBC()");
 
   // loop over non-periodic boundaries and apply BCs to J
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {

      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
      const std::string this_bc = m_bc_type[b];

      if(this_bc=="symmetry" || this_bc=="axis") {
      
         for(DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

            const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
            const Box bdry_box( bdry_grids[dit] );
            for (int dir(0); dir<SpaceDim; dir++) {
          
               // convert cell bdry box to a edge bdry box
               Box edge_box = surroundingNodes(bdry_box);
               edge_box.enclosedCells(dir);

               // collapse edge_box to 1 cell thick in bdry_dir direction                  
               if(bdry_side==0) edge_box.setSmall(bdry_dir,edge_box.bigEnd(bdry_dir));
               if(bdry_side==1) edge_box.setBig(bdry_dir,edge_box.smallEnd(bdry_dir));

               // adjust normal J boundary cells to account for
               // deposition there
               FArrayBox& Jx( a_Jx[interior_dit][dir] );
               const FArrayBox& Jx0( a_Jx0[interior_dit][dir] );
               const FArrayBox& Ex( a_Ex[interior_dit][dir] );
               const FArrayBox& Ev( a_Ev[interior_dit].getFab() );
               if(dir==bdry_dir) {
                  const FArrayBox& sigxx = a_sigma_xx[interior_dit][dir];
                  const FArrayBox& sigxy = a_sigma_xy[interior_dit][dir];
                  const FArrayBox& sigxz = a_sigma_xz[interior_dit][dir];

                  FORT_COMPUTE_JX_BDRY_FROM_MASS_MATRIX( CHF_BOX(edge_box),
                                        CHF_CONST_INT(bdry_side),
                                        CHF_CONST_FRA(sigxx),
                                        CHF_CONST_FRA(sigxy),
                                        CHF_CONST_FRA(sigxz),
                                        CHF_CONST_FRA1(Ex,0),
                                        CHF_CONST_FRA1(Ev,0),
                                        CHF_CONST_FRA1(Ev,1),
                                        CHF_CONST_FRA1(Jx0,0),
                                        CHF_FRA1(Jx,0) );
                  Jx.mult(1.0/a_volume_scale,edge_box,0,1);
               }
            }
         }
      }
   }

}

#include "NamespaceFooter.H"
