#include "SpaceUtils.H"
#include "SpaceUtilsF_F.H"
//#include "altFaceAveragesF_F.H"
//#include "upwindSchemesF_F.H"

#include "NamespaceHeader.H"

void
SpaceUtils::applyBinomialFilter( LevelData<EdgeDataBox>&  a_var )
{
   const IntVect var_gv = a_var.ghostVect();
   CH_assert( var_gv >= IntVect::Unit );
   CH_TIME("SpaceUtils::applyBinomialFilter() EdgeDataBox");
  
   const DisjointBoxLayout& grids( a_var.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir(0); dir<SpaceDim; dir++) {

         FArrayBox& this_var( a_var[dit][dir] );

         Box edge_box = grids[dit];
         edge_box.surroundingNodes( );  // grow hi end by one in all dirs
         edge_box.enclosedCells( dir ); // shrink hi end by 1 in dir
	 applyBinomialFilter( this_var, edge_box );
    
      }

   }

}

void
SpaceUtils::applyBinomialFilter( LevelData<NodeFArrayBox>&  a_var )
{
   const IntVect var_gv = a_var.ghostVect();
   CH_assert( var_gv >= IntVect::Unit );
   CH_TIME("SpaceUtils::applyBinomialFilter() NodeFArrayBox");
  
   const DisjointBoxLayout& grids( a_var.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FArrayBox& this_var( a_var[dit].getFab() );

      Box node_box = grids[dit];
      node_box.surroundingNodes( );
      applyBinomialFilter( this_var, node_box );

   }

}

void
SpaceUtils::applyBinomialFilter( FArrayBox&  a_Q,
	                   const Box&        a_grid_box	)
{
   const int ncomp = a_Q.nComp();
   FArrayBox Q2( a_Q.box(), ncomp );

   const IntVect ivzero = IntVect::Zero;
   BoxIterator bit(a_grid_box);
   for (bit.begin(); bit.ok(); ++bit) {

      IntVect iv = bit();
      for (int n=0; n<ncomp; n++) {
#if CH_SPACEDIM==1
         IntVect iup0 = ivzero; iup0[0]++;  // [1 2 1]/4
         IntVect idn0 = ivzero; idn0[0]--;
         Q2(iv,n) = a_Q(iv+iup0,n) + a_Q(iv+idn0,n);
#elif CH_SPACEDIM==2
         IntVect iup0 = ivzero; iup0[0]++; // [1 2 1
         IntVect idn0 = ivzero; idn0[0]--; //  2 4 2
         IntVect iup1 = ivzero; iup1[1]++; //  1 2 1]/16
         IntVect idn1 = ivzero; idn1[1]--;
	 Q2(iv,n) = 2.0*(a_Q(iv+iup0,n) + a_Q(iv+idn0,n))
	          + 2.0*(a_Q(iv+iup1,n) + a_Q(iv+idn1,n))
	          + a_Q(iv+iup0+iup1,n) + a_Q(iv+iup0+idn1,n)
	          + a_Q(iv+idn0+iup1,n) + a_Q(iv+idn0+idn1,n);
#elif CH_SPACEDIM==3
         IntVect iup0 = ivzero; iup0[0]++; // [1 2 1 ... 2 4 2 ... 1 2 1
         IntVect idn0 = ivzero; idn0[0]--; //  2 4 2 ... 4 8 4 ... 2 4 2
         IntVect iup1 = ivzero; iup1[1]++; //  1 2 1 ... 2 4 2 ... 1 2 1]/64
         IntVect idn1 = ivzero; idn1[1]--;
         IntVect iup2 = ivzero; iup2[2]++;
         IntVect idn2 = ivzero; idn2[2]--;
	 Q2(iv,n) = 4.0*(a_Q(iv+iup0,n) + a_Q(iv+idn0,n))
	          + 4.0*(a_Q(iv+iup1,n) + a_Q(iv+idn1,n))
	          + 4.0*(a_Q(iv+iup2,n) + a_Q(iv+idn2,n))
	          + 2.0*(a_Q(iv+iup0+iup1,n) + a_Q(iv+iup0+idn1,n))
	                       + 2.0*(a_Q(iv+iup0+iup2,n) + a_Q(iv+iup0+idn2,n))
	                       + 2.0*(a_Q(iv+iup1+iup2,n) + a_Q(iv+iup1+idn2,n))
	                       + 2.0*(a_Q(iv+idn1+iup2,n) + a_Q(iv+idn1+idn2,n))
	                       + 2.0*(a_Q(iv+idn0+iup1,n) + a_Q(iv+idn0+idn1,n))
	                       + 2.0*(a_Q(iv+idn0+iup2,n) + a_Q(iv+idn0+idn2,n))
	                       + a_Q(iv+iup0+iup1+iup2,n) + a_Q(iv+iup0+iup1+idn2,n)
	                       + a_Q(iv+iup0+idn1+iup2,n) + a_Q(iv+iup0+idn1+idn2,n)
	                       + a_Q(iv+idn0+iup1+iup2,n) + a_Q(iv+idn0+iup1+idn2,n)
	                       + a_Q(iv+idn0+idn1+iup2,n) + a_Q(iv+idn0+idn1+idn2,n);
#endif
      } // end loop over components

   } // end loop over grid cells

   const Real factor0 = std::pow(2.0,SpaceDim);
   const Real factor1 = std::pow(2.0,2*SpaceDim);
   a_Q.mult(factor0);
   a_Q.plus(Q2);
   a_Q.divide(factor1);
      
}
      
void
SpaceUtils::upWindToFaces( LevelData<FluxBox>&    a_face_phi,
                     const LevelData<FArrayBox>&  a_cell_phi,
                     const LevelData<FluxBox>&    a_norm_vel,
                     const std::string&           a_method )
{
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FluxBox& this_face_phi( a_face_phi[dit] );
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FluxBox& this_norm_vel( a_norm_vel[dit] );
      const Box& this_dbl_box( grids[dit] ); // this box has no ghost cells
      upWindToFaces( this_face_phi, this_cell_phi, this_norm_vel, this_dbl_box, a_method );
       
   } 

}

void
SpaceUtils::upWindToFaces( FluxBox&      this_face_phi,
                     const FArrayBox&    this_cell_phi,
                     const FluxBox&      a_norm_vel,
                     const Box&          a_dbl_box,
                     const std::string&  a_method )
{
   CH_TIME("SpaceUtils::upWindToFaces()");
   CH_assert( a_method=="c2" ||
              a_method=="TVD1st" || a_method=="TVDvanleer" ||
              a_method=="TVDminmod"  || a_method=="TVDsuperbee" ||
              a_method=="uw1" || a_method=="uw3" || a_method=="uw5" ||
              a_method=="quick" || a_method=="weno5" || a_method=="bweno");

   int this_limiter;
   if(a_method=="TVD1st")      this_limiter=0;
   if(a_method=="TVDvanleer")  this_limiter=1;
   if(a_method=="TVDminmod")   this_limiter=2;
   if(a_method=="TVDsuperbee") this_limiter=3;


   int normVelcomp = 0;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      if(a_norm_vel.nComp()==SpaceDim) normVelcomp = dir; 

      Box face_box( a_dbl_box ); 
    
      // for 4th order, need extra face in mapped-grid,
      // transverse direction to handle 4th-order products  
      //for (int tdir(0); tdir<SpaceDim; tdir++) {
      //   if (tdir!=dir) {
      //      const int TRANSVERSE_GROW(1);
      //      face_box.grow( tdir, TRANSVERSE_GROW );
      //   }
      //}
      
      face_box.surroundingNodes( dir );
      //face_box.grow( dir, 1 );

      //  now compute limited face values
         
      FArrayBox& this_face_phi_dir( this_face_phi[dir] );
      const FArrayBox& this_norm_vel_dir( a_norm_vel[dir] );
      if(a_method=="c2") {
         FORT_C2FACE( CHF_FRA( this_face_phi_dir ),
                      CHF_CONST_FRA( this_cell_phi ),
                      CHF_BOX( face_box ),
                      CHF_CONST_INT( dir ) );
      } else
      if(a_method=="uw1") {
         FORT_UW1FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      } else
      if(a_method=="uw3") {
         FORT_UW3FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      } else
      if(a_method=="uw5") {
         FORT_UW5FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      } else
      if(a_method=="quick") {
         FORT_QUICKFACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
      if(a_method=="weno5") {
         FArrayBox this_unit_smooth(this_cell_phi.box(),SpaceDim);
         this_unit_smooth.setVal(1.0);  
         FORT_WENO5FACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_CONST_FRA( this_unit_smooth ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
      if(a_method=="bweno") {
         FORT_BWENOFACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
      if(a_method=="TVD1st"    || a_method=="TVDvanleer" || 
         a_method=="TVDminmod" || a_method=="TVDsuperbee") {
         FORT_TVDFACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( this_limiter ),
                       CHF_CONST_INT( dir ) );
      }
   } // end loop over directions

}

void
SpaceUtils::interpToFaces( LevelData<FluxBox>&    a_face_phi,
                     const LevelData<FArrayBox>&  a_cell_phi,
                     const LevelData<FArrayBox>&  a_cell_fun,
                     const LevelData<FluxBox>&    a_norm_vel,
                     const std::string&           a_method )
{
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FluxBox& this_face_phi( a_face_phi[dit] );
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FArrayBox& this_cell_fun( a_cell_fun[dit] );
      const FluxBox& this_norm_vel( a_norm_vel[dit] );
      const Box& this_dbl_box( grids[dit] ); // this box has no ghost cells
      interpToFaces( this_face_phi, this_cell_phi, this_cell_fun, this_norm_vel, this_dbl_box, a_method );
       
   } 

}

void
SpaceUtils::interpToFaces( FluxBox&      this_face_phi,
                     const FArrayBox&    this_cell_phi,
                     const FArrayBox&    this_cell_fun,
                     const FluxBox&      a_norm_vel,
                     const Box&          a_dbl_box,
                     const std::string&  a_method )
{
   CH_TIME("SpaceUtils::interpToFaces()");
   CH_assert( a_method=="c2" || a_method == "uw1" || 
              a_method=="uw1c21st" || a_method=="uw1c2vanleer" || 
              a_method== "uw1c2minmod" || a_method == "uw1c2superbee" ||  
              a_method=="uw1c2vanAlbada1" || a_method=="uw1c2vanAlbada2" );
   int this_limiter;
   if(a_method=="uw1c21st")      this_limiter=0;
   if(a_method=="uw1c2vanleer")  this_limiter=1;
   if(a_method=="uw1c2minmod")   this_limiter=2;
   if(a_method=="uw1c2superbee") this_limiter=3;
   if(a_method=="uw1c2vanAlbada1") this_limiter=4;
   if(a_method=="uw1c2vanAlbada2") this_limiter=5;

   int normVelcomp = 0;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      if(a_norm_vel.nComp()==SpaceDim) normVelcomp = dir; 

      Box face_box( a_dbl_box );         // no ghosts 
      face_box.surroundingNodes( dir );  // grow hi end by 1 in dir direction
    
      // compute limited face values
      //
      FArrayBox& this_face_phi_dir( this_face_phi[dir] );
      const FArrayBox& this_norm_vel_dir( a_norm_vel[dir] );
      if(a_method=="c2") {
         FORT_C2FACE( CHF_FRA( this_face_phi_dir ),
                      CHF_CONST_FRA( this_cell_phi ),
                      CHF_BOX( face_box ),
                      CHF_CONST_INT( dir ) );
      } else
      if(a_method=="uw1") {
         FORT_UW1FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      } 
      else {
         FORT_UW1C2FACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA( this_cell_fun ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( this_limiter ),
                         CHF_CONST_INT( dir ) );
      }
   }
}

void
SpaceUtils::interpToFacesWENO( LevelData<FluxBox>&   a_face_phi,
                         const LevelData<FArrayBox>& a_cell_phi,
                         const LevelData<FluxBox>&   a_norm_vel,
                         const LevelData<FArrayBox>& a_smooth,
                         const std::string&          a_method )
{
   CH_TIME("SpaceUtils::interpToFacesWENO()");
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FluxBox& this_face_phi( a_face_phi[dit] );
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FluxBox& this_norm_vel( a_norm_vel[dit] );
      const FArrayBox& this_smooth( a_smooth[dit] );
      const Box& this_dbl_box( grids[dit] ); // this box has no ghost cells
      interpToFacesWENO( this_face_phi, this_cell_phi, this_norm_vel, this_smooth, this_dbl_box, a_method );
      
   } // end loop over grids

}

void
SpaceUtils::interpToFacesWENO( FluxBox&      this_face_phi,
                         const FArrayBox&    this_cell_phi,
                         const FluxBox&      a_norm_vel,
                         const FArrayBox&    this_smooth,
                         const Box&          a_dbl_box,
                         const std::string&  a_method )
{
   CH_TIME("SpaceUtils::interpToFacesWENO()");
   CH_assert( a_method=="weno5" || a_method=="bweno");

   
   //FluxBox normal_vel( a_norm_vel.box(), 1 );
   int normVelcomp = 0;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      if(a_norm_vel.nComp()==SpaceDim) normVelcomp = dir; 
      //normal_vel.copy(a_norm_vel,normVelcomp,0,1);

      Box face_box( a_dbl_box ); 
      /*   
      for (int tdir(0); tdir<SpaceDim; tdir++) {
         if (tdir!=dir) {
            const int TRANSVERSE_GROW(1);
            face_box.grow( tdir, TRANSVERSE_GROW );
         }
      }
      */
      face_box.surroundingNodes( dir );
      face_box.grow( dir, 1 );

      // Interp fluxes from cell centers to cell faces
      //
      FArrayBox& this_face_phi_dir( this_face_phi[dir] );
      //const FArrayBox& this_norm_vel_dir( normal_vel[dir] );
      const FArrayBox& this_norm_vel_dir( a_norm_vel[dir] );
      if(a_method=="weno5") {
         FORT_WENO5FACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_CONST_FRA( this_smooth ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
      if(a_method=="bweno") {
         FORT_BWENOFACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
   } // end loop over directions

}

void
SpaceUtils::interpCellToEdges( LevelData<EdgeDataBox>&   a_edge_phi,
                         const LevelData<FArrayBox>&     a_cell_phi,
                         const LevelData<EdgeDataBox>&   a_norm_vel,
                         const std::string&              a_method )
{
   const IntVect cellgv = a_cell_phi.ghostVect();
   const IntVect edgegv = a_edge_phi.ghostVect();
   CH_assert( cellgv >= edgegv || cellgv >= IntVect::Unit );
   //CH_assert( a_cell_phi.ghostVect() >= a_edge_phi.ghostVect() );
   //CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );  
   CH_assert( a_method=="c2" ); 
  
   // Note that EdgeDataBox has extra value in each of the non-dir directions,
   // opposed to FluxBox which has extra value in dir direction. 

   const DisjointBoxLayout& grids( a_edge_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      EdgeDataBox& this_edge_phi( a_edge_phi[dit] );
      
      for (int dir(0); dir<SpaceDim; dir++) {
         Box edge_box;
         if(cellgv >= edgegv )
            edge_box = this_cell_phi.box();    // has ghost 
         else {
            edge_box = grids[dit];    // no ghost 
         }
         edge_box.surroundingNodes( );  // grow hi end by one in all dirs
         edge_box.enclosedCells( dir ); // shrink hi end by 1 in dir
         
         FArrayBox& this_edge_phi_dir( this_edge_phi[dir] );
         FORT_C2EDGE( CHF_FRA( this_edge_phi_dir ),
                      CHF_CONST_FRA( this_cell_phi ),
                      CHF_BOX( edge_box ),
                      CHF_CONST_INT( dir ) );
      } // end loop over directions
   } // end loop over grids

}

void
SpaceUtils::interpEdgesToCell( LevelData<FArrayBox>&    a_cell_phi,
                         const LevelData<EdgeDataBox>&  a_edge_phi,
                         const std::string&             a_method )
{
   CH_TIME("SpaceUtils::interpEdgesToCell()");
   const int cell_ncomp = a_cell_phi.nComp();
   const int edge_ncomp = a_edge_phi.nComp();
   CH_assert( (cell_ncomp==SpaceDim && edge_ncomp==1 )
            || cell_ncomp==edge_ncomp );
   CH_assert( a_edge_phi.ghostVect()>=a_cell_phi.ghostVect() );
   CH_assert( a_method=="c2" );

   const DisjointBoxLayout& grids( a_edge_phi.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& this_cell_phi( a_cell_phi[dit] );
      Box cell_box( this_cell_phi.box() ); // has ghosts

      if(cell_ncomp==SpaceDim && edge_ncomp==1) { // edge_phi is vector
         for (int dir(0); dir<SpaceDim; dir++) {
            const FArrayBox& this_edge_phi_dir( a_edge_phi[dit][dir] );
            if(a_method=="c2") {
               FORT_C2CELL( CHF_BOX( cell_box ),
                            CHF_CONST_INT( dir ),
                            CHF_CONST_FRA1( this_edge_phi_dir,0 ),
                            CHF_FRA1( this_cell_phi,dir ) );
            }
         }
      }
      else { // edge_phi is scalar
         this_cell_phi.setVal(0.0);
         for (int dir(0); dir<SpaceDim; dir++) {
            const FArrayBox& this_edge_phi_dir( a_edge_phi[dit][dir] );
            if(a_method=="c2") {
               FORT_EDGE_SCALAR_TO_CELL( CHF_BOX( cell_box ),
                                         CHF_CONST_INT( dir ),
                                         CHF_CONST_FRA( this_edge_phi_dir ),
                                         CHF_FRA( this_cell_phi ) );
            }
         }
      }

   }

}

void
SpaceUtils::interpEdgesToEdges( LevelData<EdgeDataBox>&  a_edge_out,
                          const LevelData<EdgeDataBox>&  a_edge_in,
                          const std::string&             a_method )
{
   CH_TIME("SpaceUtils::interpEdgesToEdges()");
   CH_assert( a_edge_in.nComp()==1 );
   CH_assert( a_edge_in.ghostVect() >= 1*IntVect::Unit);
   CH_assert( a_method=="c2" );
   
   // interpolate values that live on one edge to the other edges
   // example: Ein = [Ex_onz, Ez_onx], Eout = [Ex_onx, Ez_onz]
   //
   // Hack for 2D right now... need to think about 3D
   //
   
   int dir0=1;
   const DisjointBoxLayout& grids( a_edge_out.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir(0); dir<SpaceDim; dir++) {
         
         Box edge_box( grids[dit] ); // no ghosts
         edge_box.surroundingNodes( );  // grow hi end by one in all dirs
         edge_box.enclosedCells( dir ); // shrink hi end by 1 in dir
        
         if(dir==1) dir0=0;
         const FArrayBox& this_edge_in(  a_edge_in[dit][dir0] );
               FArrayBox& this_edge_out( a_edge_out[dit][dir] );
        
         FORT_C2_EDGES_TO_EDGES( CHF_BOX( edge_box ),
                                 CHF_CONST_FRA1( this_edge_in,0 ),
                                 CHF_FRA1( this_edge_out,0 ),
                                 CHF_CONST_INT( dir ) );
      }

   }

}

void
SpaceUtils::interpNodesToEdges( LevelData<EdgeDataBox>&    a_edgePhi_out,
                          const LevelData<NodeFArrayBox>&  a_nodePhi_in,
                          const std::string&               a_method )
{
   CH_TIME("SpaceUtils::interpNodesToEdges()");
   CH_assert( a_edgePhi_out.nComp()==a_nodePhi_in.nComp() );
   CH_assert( a_method=="c2" );
   
   // interpolate values that live cell nodes to cell edges
   //
   const DisjointBoxLayout& grids( a_edgePhi_out.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir(0); dir<SpaceDim; dir++) {
         
         Box edgebox( grids[dit] ); // no ghosts
         edgebox.surroundingNodes( );  // grow hi end by one in all dirs
         edgebox.enclosedCells( dir ); // shrink hi end by 1 in dir
        
         const FArrayBox& this_node_in(  a_nodePhi_in[dit].getFab() );
               FArrayBox& this_edge_out( a_edgePhi_out[dit][dir] );
        
         FORT_C2_NODES_TO_EDGES( CHF_BOX( edgebox ),
                                 CHF_CONST_INT( dir ),
                                 CHF_CONST_FRA( this_node_in ),
                                 CHF_FRA( this_edge_out ) );
      }

   }
   
}

void
SpaceUtils::interpCellsToNodes( LevelData<NodeFArrayBox>&  a_node_phi,
                          const LevelData<FArrayBox>&      a_cell_phi,
                          const std::string&               a_method )
{
   const IntVect cellgv = a_cell_phi.ghostVect();
   const IntVect nodegv = a_node_phi.ghostVect();
   CH_assert( cellgv > nodegv || cellgv >= IntVect::Unit );
   CH_assert( a_cell_phi.nComp() == a_node_phi.nComp() );
   CH_assert( a_method=="c2" );
   
   
   const DisjointBoxLayout& grids( a_node_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
            FArrayBox& this_node_phi( a_node_phi[dit].getFab() );
      
      Box node_box;
      if(cellgv > nodegv )
         node_box = this_cell_phi.box();    // has ghost 
      else {  
         node_box = grids[dit];    // no ghost 
      }
      node_box.surroundingNodes( );  // grow hi end by one in all dirs

      FORT_C2_CELLS_TO_NODES( CHF_BOX( node_box ),
                             CHF_CONST_FRA( this_cell_phi ),
                             CHF_FRA( this_node_phi ) );
   
   }

}

void
SpaceUtils::interpNodesToCells( LevelData<FArrayBox>&      a_cellPhi_out,
                          const LevelData<NodeFArrayBox>&  a_nodePhi_in,
                          const std::string&               a_method )
{
   CH_TIME("SpaceUtils::interpNodesToCells()");
   CH_assert( a_cellPhi_out.nComp()==a_nodePhi_in.nComp() );
   CH_assert( a_nodePhi_in.ghostVect() >= a_cellPhi_out.ghostVect() );
   CH_assert( a_method=="c2" );
   
   // interpolate values that live cell nodes to cell centers
   //
   const DisjointBoxLayout& grids( a_cellPhi_out.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& this_node_in(  a_nodePhi_in[dit].getFab() );
            FArrayBox& this_cell_out( a_cellPhi_out[dit] );
      
      //Box cellbox( grids[dit] ); // no ghosts
      //cellbox.grow( 1 );         // add one layer of ghost cells
      const Box& cellbox( this_cell_out.box() );        

      FORT_C2_NODES_TO_CELLS( CHF_BOX( cellbox ),
                              CHF_CONST_FRA( this_node_in ),
                              CHF_FRA( this_cell_out ) );
   }

}

void
SpaceUtils::PerpEdgeGradientAtCells( FArrayBox&    a_gradPhi,
                               const EdgeDataBox&  a_edgePhi,
                               const RealVect&     a_dx,
                               const Box&          a_box,
                               const std::string&  a_method )
{
   CH_TIME("SpaceUtils::PerpEdgeGradientAtCells()");
   CH_assert( a_gradPhi.nComp()==SpaceDim );
   CH_assert( a_method=="c2" );
   
   // calculate cross gradient comps at cells from edge data
   // i.e.
   // gradPhi(comp=0) = d(edgePhi0)/dX1
   // gradPhi(comp=1) = d(edgePhi1)/dX0
   //
          
   int dir0 = 1;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      const FArrayBox& this_edgePhi( a_edgePhi[dir] );

      if(dir==1) dir0 = 0;
      double thisdx = a_dx[dir0];

      FORT_EDGE_GRAD_AT_CELLS( CHF_BOX( a_box ),
                               CHF_CONST_INT( dir0 ),
                               CHF_CONST_REAL( thisdx ),
                               CHF_CONST_FRA1( this_edgePhi,0 ),
                               CHF_FRA1( a_gradPhi,dir ) );
   }
      
}

void
SpaceUtils::interpFacesToCell( LevelData<FArrayBox>&  a_cell_phi,
                         const LevelData<FluxBox>&    a_face_phi,
                         const std::string&           a_method )
{
   CH_TIME("SpaceUtils::interpFacesToCell()");
   CH_assert( (a_cell_phi.nComp() == SpaceDim && a_face_phi.nComp() == 1) ||
              (a_cell_phi.nComp() == a_face_phi.nComp()) );
   CH_assert( a_method=="c2" );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& this_cell_phi( a_cell_phi[dit] );
      Box cell_box( grids[dit] );
      if(a_cell_phi.nComp() == SpaceDim && a_face_phi.nComp() == 1) {

         for (int dir(0); dir<SpaceDim; dir++) {
            const FArrayBox& this_face_phi_dir( a_face_phi[dit][dir] );
            FORT_C2FACETOCELL( CHF_BOX( cell_box ),
                               CHF_CONST_INT( dir ),
                               CHF_CONST_FRA1( this_face_phi_dir,0 ),
                               CHF_FRA1( this_cell_phi,dir ) );
         }
      }
      else { // face_phi is scalar
      
         this_cell_phi.setVal(0.0);
         for (int dir(0); dir<SpaceDim; dir++) {
            const FArrayBox& this_face_phi_dir( a_face_phi[dit][dir] );
            FORT_FACE_SCALAR_TO_CELL( CHF_BOX( cell_box ),
                                      CHF_CONST_INT( dir ),
                                      CHF_CONST_FRA( this_face_phi_dir ),
                                      CHF_FRA( this_cell_phi ) );
         }                            
      }  
      
   } 

}

void
SpaceUtils::computeLaxSplitting( LevelData<FArrayBox>&  a_fluxR,
                                 LevelData<FArrayBox>&  a_fluxL,
                           const LevelData<FArrayBox>&  a_flux,
                           const LevelData<FArrayBox>&  a_Cspeed,
                           const LevelData<FArrayBox>&  a_fun )
{
   CH_assert( a_flux.nComp()==a_Cspeed.nComp() );

   const DisjointBoxLayout& grids( a_fun.getBoxes() );

   //   get left and right going flux at cell-center
   //   fluxR = 0.5*(flux + Cspeed*fun),
   //   fluxL = 0.5*(flux - Cspeed*fun),
   //   Cspeed = abs(max(eigenValue of Flux Jacobian))
   //
   a_fluxR.define(grids, a_flux.nComp(), a_flux.ghostVect());
   a_fluxL.define(grids, a_flux.nComp(), a_flux.ghostVect());
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_fluxR[dit].copy(a_Cspeed[dit]);
      a_fluxL[dit].copy(a_Cspeed[dit]);
      a_fluxL[dit].mult(-1.0);
      for (int n(0); n < a_flux[dit].nComp(); n++) {
         a_fluxR[dit].mult(a_fun[dit],0,n,1);
         a_fluxL[dit].mult(a_fun[dit],0,n,1);
      }
      a_fluxR[dit].plus(a_flux[dit]);
      a_fluxL[dit].plus(a_flux[dit]);
      a_fluxR[dit].mult(0.5);
      a_fluxL[dit].mult(0.5);
   }

}

void
SpaceUtils::computeLaxSplitting( FArrayBox&  a_fluxR,
                                 FArrayBox&  a_fluxL,
                           const FArrayBox&  a_flux,
                           const FArrayBox&  a_Cspeed,
                           const FArrayBox&  a_fun )
{
   CH_TIME("SpaceUtils::computeLaxSplitting()");
   CH_assert( a_flux.nComp()==a_Cspeed.nComp() );

   const Box& thisbox = a_fun.box();

   //   get left and right going flux at cell-center
   //   fluxR = 0.5*(flux + Cspeed*fun),
   //   fluxL = 0.5*(flux - Cspeed*fun),
   //   Cspeed = abs(max(eigenValue of Flux Jacobian))
   //
   a_fluxR.define(thisbox, a_flux.nComp());
   a_fluxL.define(thisbox, a_flux.nComp());

   a_fluxR.copy(a_Cspeed, a_fluxR.box());
   a_fluxL.copy(a_Cspeed, a_fluxL.box());
   //a_fluxL.mult(-1.0);
   a_fluxL.negate();
   for (int n(0); n < a_flux.nComp(); n++) {
      a_fluxR.mult(a_fun,0,n,1);
      a_fluxL.mult(a_fun,0,n,1);
   }
   a_fluxR.plus(a_flux);
   a_fluxL.plus(a_flux);
   a_fluxR.mult(0.5);
   a_fluxL.mult(0.5);

}

void
SpaceUtils::computeLaxSplitting( LevelData<FArrayBox>&  a_fluxR,
                                 LevelData<FArrayBox>&  a_fluxL,
                           const LevelData<FArrayBox>&  a_flux,
                           const LevelData<FArrayBox>&  a_Cspeed,
                           const LevelData<FArrayBox>&  a_fun,
                           const int  a_fun_comp )
{

   const DisjointBoxLayout& grids( a_fun.getBoxes() );
   LevelData<FArrayBox> this_a_fun(grids, 1, a_fun.ghostVect());
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      this_a_fun[dit].copy(a_fun[dit],a_fun_comp,0,1);
   }
   computeLaxSplitting(a_fluxR,a_fluxL,a_flux,a_Cspeed,this_a_fun);

}

void
SpaceUtils::faceInterpolate(  const int         a_dir,
                              const Box&        a_fc_box,
                              const int         a_order,
                              const FArrayBox&  a_var,
                              FArrayBox&        a_face_var )
{
  int ncomp = a_var.nComp();
  CH_assert(a_face_var.nComp() == ncomp);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert (a_dir >= 0 && a_dir < SpaceDim);

  for (int n=0; n<ncomp; n++) {
    FORT_FACE_INTERPOLATE(  CHF_CONST_INT(a_dir),
                            CHF_BOX(a_fc_box),
                            CHF_CONST_INT(a_order),
                            CHF_CONST_FRA1(a_var, n),
                            CHF_FRA1(a_face_var, n) );
  }

  return;
}


void
SpaceUtils::faceInterpolate(  const int         a_dir,
                              const Box&        a_fc_box,
                              const Box&        a_cc_box,
                              const int         a_order,
                              const FArrayBox&  a_var,
                              FArrayBox&        a_face_var )
{
  int ncomp = a_var.nComp();
  CH_assert(a_face_var.nComp() == ncomp);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert (a_dir >= 0 && a_dir < SpaceDim);

  int imin = a_cc_box.smallEnd(a_dir),
      imax = a_cc_box.bigEnd(a_dir);

  if (imin == imax) {

    /* box is flat in along this dimension,
     * so copy instead of interpolate */
    BoxIterator bit(a_cc_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      IntVect iw = iv; iw[a_dir]++;

      for (int n=0; n<ncomp; n++) {
        a_face_var(iv,n) = a_var(iv,n);
        a_face_var(iw,n) = a_var(iv,n);
      }
    }

  } else {

    faceInterpolate(a_dir, a_fc_box, a_order, a_var, a_face_var);

  }

  return;
}

void
SpaceUtils::interpolateStag( FArrayBox&  a_dst, 
                       const Box&        a_dst_box,
                       const int         a_dst_comp,
                       const FArrayBox&  a_src,
                       const int         a_src_comp,
                       const int         a_dir )
{
   // c2 interpolate from cell to node or node to cell 

   const IntVect dst_box_type = a_dst.box().type();
   const IntVect src_box_type = a_src.box().type();
   const int src_is_stag = src_box_type[a_dir];
   CH_assert(dst_box_type[a_dir]!=src_is_stag);  
   
   FORT_STAG_INTERPOLATE( CHF_BOX(a_dst_box),
                          CHF_CONST_INT(a_dir),
                          CHF_CONST_INT(src_is_stag),
                          CHF_CONST_FRA1(a_src, a_src_comp),
                          CHF_FRA1(a_dst, a_dst_comp) );

}

void 
SpaceUtils::simpleStagGradComp( FArrayBox&  a_dst,
                          const Box&        a_grid_box,
                          const int         a_dst_comp,
                          const FArrayBox&  a_src,
                          const int         a_src_comp,
                          const Real        a_dX,
                          const int         a_dir_dX,
                          const int         a_additive ) 
{
   CH_TIME("SpaceUtils::simpleDifference()");

   // compute simple gradient in a_dir_dX direction for staggered
   // variables - cell to node or node to cell

   const IntVect dst_box_type = a_dst.box().type();
   const IntVect src_box_type = a_src.box().type();
   const int src_is_stag = src_box_type[a_dir_dX];
   CH_assert(dst_box_type[a_dir_dX]!=src_is_stag);  
   
   // create the appropriate box
   Box dst_box = a_grid_box;
   for (int dir=0; dir<SpaceDim; dir++) {
      if(dst_box_type[dir]) dst_box.surroundingNodes(dir);  // grow dir hi end by one
   }

   //const Box& dst_box = a_dst.box();
   FORT_STAG_GRAD_COMPONENT( CHF_BOX(dst_box),
                             CHF_CONST_INT(a_dir_dX),
                             CHF_CONST_REAL(a_dX),
                             CHF_CONST_FRA1(a_src, a_src_comp),
                             CHF_CONST_INT(a_additive),
                             CHF_CONST_INT(src_is_stag),
                             CHF_FRA1(a_dst, a_dst_comp) );

} 


void
SpaceUtils::cellCenteredGradientComponent(  const Box&        a_box,
                                            const int         a_dir,
                                            const FArrayBox&  a_var,
                                            const RealVect&   a_dx,
                                            const int         a_order,
                                            FArrayBox&        a_grad_var )
{
  CH_assert(a_var.nComp() == 1);
  CH_assert(a_grad_var.nComp() > a_dir);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  FORT_CELL_CENTERED_GRAD_COMPONENT(  CHF_BOX(a_box),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_FRA1(a_var,0),
                                      CHF_CONST_REALVECT(a_dx),
                                      CHF_CONST_INT(a_order),
                                      CHF_CONST_FRA1(a_grad_var, a_dir) );

  return;
}

void
SpaceUtils::faceCenteredGradientComponent(  const Box&        a_box,
                                            const int         a_dir,
                                            const FArrayBox&  a_var,
                                            const RealVect&   a_dx,
                                            const int         a_order,
                                            FArrayBox&        a_grad_var )
{
  CH_assert(a_var.nComp() == 1);
  CH_assert(a_grad_var.nComp() > a_dir);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  FORT_FACE_CENTERED_GRAD_COMPONENT(  CHF_BOX(a_box),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_FRA1(a_var, 0),
                                      CHF_CONST_REALVECT(a_dx),
                                      CHF_CONST_INT(a_order),
                                      CHF_FRA1(a_grad_var, a_dir) );

  return;
}

void
SpaceUtils::extrapBoundaryGhostsForCC( FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const int                               a_order,
                                       const int                               a_side)
{
   // If a_order = 2, this function second-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   // If a_order = 4, this function fourth-order extrapolates to fill three layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   int imin = a_interiorbox.smallEnd(a_dir);
   int imax = a_interiorbox.bigEnd(a_dir);

   int depth = (a_order==2)? 2: 3;
   
   Box dstbox;
   bool isbc = false;
      
   switch(a_side) 
     {
     case -1:
       // Handle the low side
       isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
       if (isbc) {
	 dstbox = adjCellLo(a_interiorbox, a_dir, depth);
       }
       break;
     case 1:
       // Handle high side
       isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
       if (isbc) {
	 dstbox = adjCellHi(a_interiorbox, a_dir, depth);
       }
       break;
     }

   if (isbc) {

     CH_assert(bx.contains(dstbox));
     
     if (imin == imax) {

       BoxIterator bit(dstbox);
       for (bit.begin(); bit.ok(); ++bit) {
	 IntVect i1 = bit();
	 IntVect i2 = i1; i2[a_dir] = imin;
	 for (int n=0; n<a_data.nComp(); n++) {
	   a_data(i1,n) = a_data(i2,n);
	 }
       }
       
     } else {

       CH_assert( (a_order==4 && a_interiorbox.size(a_dir)>=5)
		  || (a_order==2 && a_interiorbox.size(a_dir)>=3));
       FORT_EXTRAP_FOR_CC_OPS(CHF_CONST_INT(a_dir),
			      CHF_CONST_INT(a_side),
			      CHF_CONST_INT(a_order),
			      CHF_BOX(dstbox),
			      CHF_BOX(a_interiorbox),
			      CHF_FRA(a_data));

     }
   }

   return;
} 

void
SpaceUtils::extrapBoundaryGhostsForFC( FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const int                               a_order,
				       const int                               a_side)
{
   // This function fourth-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   int imin = a_interiorbox.smallEnd(a_dir);
   int imax = a_interiorbox.bigEnd(a_dir);

   Box dstbox;
   bool isbc = false;
      
   switch(a_side) 
     {
     case -1:
       // Handle the low side
       isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
       if (isbc) {
	 dstbox = adjCellLo(a_interiorbox, a_dir, 2);
       }
       break;
     case 1:
       // Handle high side
       isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
       if (isbc) {
	 dstbox = adjCellHi(a_interiorbox, a_dir, 2);
       }
       break;
     }
   
   if (isbc) {
     
     CH_assert(bx.contains(dstbox));

     if (imin == imax) {

       BoxIterator bit(dstbox);
       for (bit.begin(); bit.ok(); ++bit) {
	 IntVect i1 = bit();
	 IntVect i2 = i1; i2[a_dir] = imin;
	 for (int n=0; n<a_data.nComp(); n++) {
	   a_data(i1,n) = a_data(i2,n);
	 }
       }
       
     } else {
       
       CH_assert( (a_order==4 && a_interiorbox.size(a_dir)>=5)
		  || (a_order==2 && a_interiorbox.size(a_dir)>=3));
       FORT_EXTRAP_FOR_FC_OPS(CHF_CONST_INT(a_dir),
			      CHF_CONST_INT(a_side),
			      CHF_CONST_INT(a_order),
			      CHF_BOX(dstbox),
			      CHF_BOX(a_interiorbox),
			      CHF_FRA(a_data));
       
       
     }
   }

   return;
} 


void
SpaceUtils::secondOrderTransExtrapAtDomainBdry(FArrayBox&           a_data,
                                               const int            a_dir,
                                               const Box&           a_interiorbox,
                                               const ProblemDomain& a_domain,
                                               const int            a_maxdim )
{
   CH_TIME("secondOrderTransExtrapAtDomainBdry");

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));
   const Box& domBox = a_domain.domainBox();

   for (int tdir=0; tdir<a_maxdim; ++tdir)
   {
      if (tdir != a_dir && !a_domain.isPeriodic(tdir))
      {
         CH_assert(a_interiorbox.size(tdir)>=3);

         for (int side=-1; side<2; side+=2)
         {
            Box dstbox;
            bool isbc = false;

            switch(side)
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(tdir) == domBox.smallEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellLo(a_interiorbox,tdir,1);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(tdir) == domBox.bigEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellHi(a_interiorbox,tdir,1);
               }
               break;
            }

            if (isbc)
            {
               CH_assert(bx.contains(dstbox));
               FORT_SECOND_ORDER_EXTRAPOLATION( CHF_CONST_INT(tdir),
                                                CHF_CONST_INT(side),
                                                CHF_BOX(a_interiorbox),
                                                CHF_BOX(dstbox),
                                                CHF_FRA(a_data) );
            }
         }
      }
   }

   return;
}



void 
SpaceUtils::fillGhostCellsSimple( FArrayBox& a_phi,
                                  const Box& a_bx_int,
                                  const int  a_dir )
{
  CH_assert((a_dir >=0) && (a_dir < SpaceDim));
  const Box& phi_bx = a_phi.box();

  for (int side=-1; side<2; side+=2) {

    int ii;
    Box bdrybox;

    int n_gpt;
    if (side == -1) {
      n_gpt = a_bx_int.smallEnd(a_dir) - phi_bx.smallEnd(a_dir);
      if (n_gpt > 0) {
        ii = a_bx_int.smallEnd(a_dir);
        bdrybox = adjCellLo(a_bx_int, a_dir, n_gpt);
      }
    } else {
      n_gpt = phi_bx.bigEnd(a_dir) - a_bx_int.bigEnd(a_dir);
      if (n_gpt > 0) {
        ii = a_bx_int.bigEnd(a_dir);
        bdrybox = adjCellHi(a_bx_int, a_dir, n_gpt);
      }
    }

    if (n_gpt > 0) {
      BoxIterator bit(bdrybox);
      for (bit.begin(); bit.ok(); ++bit) {
  
        IntVect iv(bit());
        int j = (Real) (side < 0 ? ii-iv[a_dir] : iv[a_dir]-ii );
  
        IntVect i_int_0(iv); i_int_0[a_dir] = ii;
        IntVect i_int_1(iv); i_int_1[a_dir] = ii - side;
        IntVect i_int_2(iv); i_int_2[a_dir] = ii - 2*side;
        IntVect i_int_3(iv); i_int_3[a_dir] = ii - 3*side;
  
        Real c0 = ((1.0+j)*(2.0+j)*(3.0+j))/6.0;
        Real c1 = -(j*(2.0+j)*(3.0+j))/2.0;
        Real c2 = (j*(1.0+j)*(3.0+j))/2.0;
        Real c3 = -(j*(1.0+j)*(2.0+j))/6.0;
  
        for (int n=0; n<a_phi.nComp(); n++) {
          a_phi(iv,n) =   c0 * a_phi(i_int_0,n)
                        + c1 * a_phi(i_int_1,n)
                        + c2 * a_phi(i_int_2,n)
                        + c3 * a_phi(i_int_3,n);
        }
  
      }
    }
  }

  return;
}

void
SpaceUtils::copyAndFillGhostCellsSimple(  LevelData<FArrayBox>&       a_var_wg,
                                          const LevelData<FArrayBox>& a_var )
{
  const DisjointBoxLayout& var_dbl = a_var.disjointBoxLayout();

  /* interior */
  for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var_wg[dit].setVal(0.0);
    a_var_wg[dit].copy(a_var[dit], var_dbl[dit]);
  }
  a_var_wg.exchange();

  /* codim-1 boundaries */
  for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
    for (int dir=0; dir<SpaceDim; dir++) {
      fillGhostCellsSimple(a_var_wg[dit], var_dbl[dit], dir);
    }
  }
  a_var_wg.exchange();

  /* higher codim boundaries */
  for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
    Box bx_int(var_dbl[dit]);
    const Box& bx_var(a_var_wg[dit].box());
    for (int dir=0; dir<SpaceDim-1; dir++) {
      bx_int.growLo(dir, (bx_int.smallEnd(dir)-bx_var.smallEnd(dir)) );
      bx_int.growHi(dir, (bx_var.bigEnd(dir)-bx_int.bigEnd(dir)) );
      fillGhostCellsSimple(a_var_wg[dit], bx_int, dir+1);
    }
  }
  a_var_wg.exchange();

  /* done - hopefully! */
  return;
}

void 
SpaceUtils::inspectFArrayBox( const FArrayBox&  a_F0,
                              const Box&        a_box,
                              const int         a_comp )
{
   CH_assert(a_comp<a_F0.nComp());
   FORT_INSPECT_FARRAYBOX( CHF_BOX(a_box), 
                           CHF_CONST_FRA1(a_F0,a_comp) );
}

void 
SpaceUtils::inspectFArrayBox(const LevelData<FArrayBox>&  a_F0,
                             const int                    a_comp, 
                             const int                    a_ghosts )
{
   CH_assert(a_comp<a_F0.nComp());
   CH_assert(a_F0.ghostVect() >= a_ghosts*IntVect::Unit);
   const DisjointBoxLayout& grids( a_F0.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& F0_on_patch = a_F0[dit];   
      Box cellbox( grids[dit] ); // no ghosts
      if(a_ghosts) cellbox.grow( a_ghosts );

      FORT_INSPECT_FARRAYBOX( CHF_BOX(cellbox), 
                              CHF_CONST_FRA1(F0_on_patch,a_comp) );
   }
}

void 
SpaceUtils::inspectFluxBox(const LevelData<FluxBox>&  a_Flux,
                           const int                  a_dir)
{
   const DisjointBoxLayout& grids( a_Flux.getBoxes() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& Flux_on_dir = a_Flux[dit][a_dir]; 
      Box facebox = Flux_on_dir.box();
      //Box facebox( grids[dit] );
      //facebox.surroundingNodes( a_dir );

      FORT_INSPECT_FLUXBOX( CHF_BOX(facebox), 
                            CHF_CONST_FRA(Flux_on_dir),
                            CHF_CONST_INT(a_dir) );
   }
}

void 
SpaceUtils::inspectEdgeDataBox(const LevelData<EdgeDataBox>&  a_Edge,
                               const int                      a_dir)
{
   const DisjointBoxLayout& grids( a_Edge.getBoxes() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const EdgeDataBox& Edge_on_patch = a_Edge[dit]; 
      const FArrayBox& Edge_on_dir = Edge_on_patch[a_dir]; 

      const Box& edgebox = Edge_on_dir.box();
      FORT_INSPECT_FLUXBOX( CHF_BOX(edgebox), 
                            CHF_CONST_FRA(Edge_on_dir),
                            CHF_CONST_INT(a_dir) );
   }
}

void 
SpaceUtils::inspectNodeFArrayBox(const LevelData<NodeFArrayBox>&  a_F0,
                                 const int                        a_comp )
{
   CH_assert(a_comp<a_F0.nComp());
   const DisjointBoxLayout& grids( a_F0.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FArrayBox& F0_on_patch = a_F0[dit].getFab();   
      const Box node_box = F0_on_patch.box();
      FORT_INSPECT_FARRAYBOX( CHF_BOX(node_box), 
                              CHF_CONST_FRA1(F0_on_patch,a_comp) );
   }

}

#if 1
void
SpaceUtils::copyEdgeDataBox( EdgeDataBox& a_dst,
		       const EdgeDataBox& a_src )
{
  
  CH_assert(a_dst.nComp() == a_src.nComp());
  for (int dir=0; dir<SpaceDim; dir++) {
     FArrayBox& dst_dir = a_dst[dir];
     const FArrayBox& src_dir = a_src[dir];
   
     Box thisbox(dst_dir.box()); // use the smaller box
     if(!src_dir.box().contains(dst_dir.box())) thisbox = src_dir.box();

     copy( dst_dir, src_dir, thisbox );
  }

}

void
SpaceUtils::copy(  FArrayBox& a_dst,
		   const FArrayBox& a_src )
{
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert(a_src.box().contains(a_dst.box()));
  FORT_COPY(CHF_BOX(a_dst.box()), 
	    CHF_FRA( a_dst ),
	    CHF_CONST_FRA( a_src ));

}

void
SpaceUtils::copy(  FArrayBox& a_dst,
             const FArrayBox& a_src,
             const Box&       a_box )
{
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert(a_dst.box().contains(a_box));
  CH_assert(a_src.box().contains(a_box));
  FORT_COPY(CHF_BOX( a_box ),
            CHF_FRA( a_dst ),
            CHF_CONST_FRA( a_src ));

}

void
SpaceUtils::copyNodeToCell(  FArrayBox&  a_dst,
                       const FArrayBox&  a_src,
                       const Box&        a_box,
                       const int         a_dir )
{ 
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert(a_dst.box().type(a_dir) == 0);
  CH_assert(a_src.box().type(a_dir) == 1);
  //cout << "copyNodeToCell: src.box().type() = " << a_src.box().type() << endl;
  FORT_COPY(CHF_BOX( a_box ),
            CHF_FRA( a_dst ),
            CHF_CONST_FRA( a_src ));

}

void
SpaceUtils::setVal( FArrayBox&  a_dst,
              const Real        a_val,
              const int         a_comp )
{
  CH_TIME("SpaceUtils::setVal()");

  CH_assert(a_comp<a_dst.nComp());

  FORT_SETVAL( CHF_BOX(a_dst.box()), 
               CHF_CONST_REAL(a_val),
               CHF_FRA1(a_dst,a_comp) );

}

void
SpaceUtils::localVectorNorm(  FArrayBox& a_dst,
                        const FArrayBox& a_src )
{ 
  CH_assert(a_dst.nComp() == 1);
  CH_assert(a_src.box().contains(a_dst.box()));
  FORT_VECTOR_NORM( CHF_BOX(a_dst.box()),
                    CHF_FRA1( a_dst,0 ), 
                    CHF_CONST_FRA( a_src ) );

}

void
SpaceUtils::exchangeFluxBox( LevelData<FluxBox>& a_Face )
{
   CH_TIME("SpaceUtils::exchangeFluxBox()");
   
   const DisjointBoxLayout& grids( a_Face.getBoxes() );
   DataIterator dit = grids.dataIterator();

   // copy to a temporary
   LevelData<FluxBox> tmp_Face(grids, a_Face.nComp());
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) { 
         const FArrayBox& Face_on_dir     = a_Face[dit][dir]; 
               FArrayBox& tmp_Face_on_dir = tmp_Face[dit][dir]; 
         copy(tmp_Face_on_dir,Face_on_dir);
      }
   }

   // call exchange on passed FluxBox
   a_Face.exchange();

   // now use the original data stored in temp to modify
   // the exchanged Face data. In this example, I'm copying data
   // on the low-side grid edge from temp->Face, which has the
   // effect of declaring the high-side value on the shared edge
   // to be the "correct" value

   for (dit.begin(); dit.ok(); ++dit) {
      FluxBox& thisTemp = tmp_Face[dit];
      FluxBox& thisFace = a_Face[dit];
      for (int dir=0; dir<SpaceDim; dir++) {
         Box loFaceBox=thisTemp[dir].box();
         loFaceBox.setBig(dir,loFaceBox.smallEnd(dir));
         copy(thisFace[dir],thisTemp[dir],loFaceBox);
      }
   }
 
   // exchange again is for the purpose of actual ghost cells (not shared cells)
   a_Face.exchange();

}

void
SpaceUtils::exchangeEdgeDataBox( LevelData<EdgeDataBox>& a_Edge )
{
   CH_TIME("SpaceUtils::exchangeEdgeDataBox()");
   
#if CH_SPACEDIM==1
   a_Edge.exchange();
#else
   const DisjointBoxLayout& grids( a_Edge.getBoxes() );
   DataIterator dit = grids.dataIterator();

   // copy to a temporary
   LevelData<EdgeDataBox> tmp_Edge(grids, a_Edge.nComp());
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) { 
         const FArrayBox& Edge_on_dir     = a_Edge[dit][dir]; 
               FArrayBox& tmp_Edge_on_dir = tmp_Edge[dit][dir]; 
         copy(tmp_Edge_on_dir,Edge_on_dir);
      }
   }

   // call exchange on passed EdgeDataBox
   a_Edge.exchange();

   // now use the original data stored in temp to modify
   // the exchanged Edge data. In this example, I'm copying data
   // on the low-side grid edge from temp->Edge, which has the
   // effect of declaring the high-side value on the shared edge
   // to be the "correct" value

   for (dit.begin(); dit.ok(); ++dit) {
      EdgeDataBox& thisTemp = tmp_Edge[dit];
      EdgeDataBox& thisEdge = a_Edge[dit];
      for (int dir=0; dir<SpaceDim; dir++) {
         for (int dir0=0; dir0<SpaceDim; dir0++) {
             if(dir0!=dir) {
                Box loEdgeBox=thisTemp[dir].box();
                loEdgeBox.setBig(dir0,loEdgeBox.smallEnd(dir0));
                copy(thisEdge[dir],thisTemp[dir],loEdgeBox);
             }
         }
      }
   }
   
   // exchange again is for the purpose of actual ghost cells (not shared cells)
   a_Edge.exchange();
#endif

}

void
SpaceUtils::exchangeNodeFArrayBox( LevelData<NodeFArrayBox>& a_Node )
{
   CH_TIME("SpaceUtils::exchangeNodeFArrayBox()");
  
   // tested/verified for 2D with doubly-periodic BCs, periodic only in X,
   // periodic only in Y, and not periodic in either direction

   CH_assert(SpaceDim<3); // not generalized for 3D
                          // not tested in 1D
      
   const DisjointBoxLayout& grids( a_Node.getBoxes() );
   DataIterator dit = grids.dataIterator();

   // copy to a temporary
   LevelData<NodeFArrayBox> tmp_Node(grids, a_Node.nComp());
   for (dit.begin(); dit.ok(); ++dit) {
      const FArrayBox& this_Node = a_Node[dit].getFab(); 
            FArrayBox& this_tmp_Node = tmp_Node[dit].getFab(); 
       Box NodeBox = grids[dit];
       NodeBox.surroundingNodes();
       this_tmp_Node.copy(this_Node,NodeBox);
   }

   // call exchange on passed NodeFArrayBox
   a_Node.exchange();

   // now use the original data stored in temp to modify
   // the exchanged Node data such that all shared values are identical

   const ProblemDomain& domain( grids.physDomain() ); 
   Box domainNodeBox = domain.domainBox();
   domainNodeBox.surroundingNodes();
   const IntVect domainNodeHi = domainNodeBox.bigEnd();
   const IntVect domainNodeLo = domainNodeBox.smallEnd();
   
   bool is_periodic_any_dir = false;
   for (int dir=0; dir<SpaceDim; dir++) {
      if(domain.isPeriodic(dir)) {is_periodic_any_dir = true;}
   }

   for (dit.begin(); dit.ok(); ++dit) {
      const Box& gridBox = grids[dit];
      FArrayBox& this_Temp = tmp_Node[dit].getFab();
      FArrayBox& this_Node = a_Node[dit].getFab();

      // first copy all non-corner shared cells on low ends
      for(int dir=0; dir<SpaceDim; dir++) {
         Box loNodeBox = gridBox;
         loNodeBox.surroundingNodes();                  // convert to node type
         loNodeBox.setBig(dir,loNodeBox.smallEnd(dir)); // collapse to dir small end
         for(int dir0=0; dir0<SpaceDim; dir0++) {       // remove corners
            if(dir0!=dir) loNodeBox.growHi(dir0,-1);
            if(dir0!=dir) loNodeBox.growLo(dir0,-1);
         }
         copy(this_Node,this_Temp,loNodeBox);
      }
   
      //
      // now focus on the corners
      //
      
      Box NodeBox = gridBox;
      NodeBox.surroundingNodes();

      // check if low-low corner touchs any physical
      // boundaries. If not, then do copy
      int copyCorner00 = 1;
      for(int dir=0; dir<SpaceDim; dir++) {
         const int thisBoxLo = NodeBox.smallEnd(dir);
         if(thisBoxLo==domainNodeLo[dir] && domain.isPeriodic(dir)) {
            const int thisBoxHi = NodeBox.bigEnd(dir);
	    if(thisBoxHi!=domainNodeHi[dir]) copyCorner00 = copyCorner00*0;
         }
      }
      if(copyCorner00) {
         Box corner00Box = NodeBox;
         for(int dir=0; dir<SpaceDim; dir++) {
            corner00Box.setBig(dir,NodeBox.smallEnd(dir));
         }
         copy(this_Node,this_Temp,corner00Box);
      }
      
      if(SpaceDim==2) { // this would need to be generalized for SpaceDim>2
      
         // check if high end touchs a physical in dir0, but 
         // does not touch low end in dir1
         const int thisBoxHi0 = NodeBox.bigEnd(0);
         const int thisBoxLo1 = NodeBox.smallEnd(1);
	 bool copy_corner10 = false;
         if(domainNodeHi[0]==thisBoxHi0) { // sufficient for periodic in X and not in Y
	    if(domain.isPeriodic(0) && !domain.isPeriodic(1)) { copy_corner10 = true; }
	    else if(domainNodeLo[1]!=thisBoxLo1) { copy_corner10 = true; }
	    if(copy_corner10) {
            Box corner10Box = NodeBox;
            corner10Box.setSmall(0,NodeBox.bigEnd(0));
            corner10Box.setBig(1,NodeBox.smallEnd(1));
            copy(this_Node,this_Temp,corner10Box);
	    }
         }
      
         // check if high end touchs a physical in dir1, but 
         // does not touch low end in dir0
         const int thisBoxHi1 = NodeBox.bigEnd(1);
         const int thisBoxLo0 = NodeBox.smallEnd(0);
	 bool copy_corner01 = false;
         if(domainNodeHi[1]==thisBoxHi1) { // sufficient for periodic in Y and not in X
            if(domain.isPeriodic(1) && !domain.isPeriodic(0)) { copy_corner01 = true; }
	    else if(domainNodeLo[0]!=thisBoxLo0) {copy_corner01 = true; }
	    if(copy_corner01) {
            Box corner01Box = NodeBox;
            corner01Box.setSmall(1,NodeBox.bigEnd(1));
            corner01Box.setBig(0,NodeBox.smallEnd(0));
            copy(this_Node,this_Temp,corner01Box);
	    }
         }

      }

      if(is_periodic_any_dir) {

         // copy over upper domain corner
         bool copy_corner11 = true;
         Box corner11Box = NodeBox;
         for(int dir=0; dir<SpaceDim; dir++) {
            const int thisBoxHi = NodeBox.bigEnd(dir);
            if(domainNodeHi[dir]==thisBoxHi) {
               corner11Box.setSmall(dir,NodeBox.bigEnd(dir));
            }
            else { copy_corner11 = false; }
         }
         if(copy_corner11) {
	    copy(this_Node,this_Temp,corner11Box);
	 }

      }
      
   }
   
   // exchange again is for the purpose of actual ghost cells (not shared cells)
   if(is_periodic_any_dir) { a_Node.exchange(); }
   
}

void
SpaceUtils::checkForNAN( const LevelData<FArrayBox>&  a_var,
		         const std::string&           a_string )
{
  const DisjointBoxLayout& grids( a_var.getBoxes() );
  for (DataIterator dit(grids); dit.ok(); ++dit) {
    const FArrayBox& this_var = a_var[dit];
    checkForNAN(this_var,a_string);
  }
}

void
SpaceUtils::checkForNAN( const LevelData<FluxBox>&  a_var,
		         const std::string&         a_string )
{
  const DisjointBoxLayout& grids( a_var.getBoxes() );
  for (DataIterator dit(grids); dit.ok(); ++dit) {
    for (int dir = 0; dir<SpaceDim; dir++) {
      const FArrayBox& this_var = a_var[dit][dir];
      checkForNAN(this_var,a_string);
    }
  }

}

void
SpaceUtils::checkForNAN( const LevelData<EdgeDataBox>&  a_var,
		         const std::string&             a_string )
{
  const DisjointBoxLayout& grids( a_var.getBoxes() );
  for (DataIterator dit(grids); dit.ok(); ++dit) {
    for (int dir = 0; dir<SpaceDim; dir++) {
      const FArrayBox& this_var = a_var[dit][dir];
      const std::string dir_string = a_string + " at dir = " + std::to_string(dir);
      if(procID()==0) checkForNAN(this_var,dir_string);
    }
  }

}

void
SpaceUtils::checkForNAN( const LevelData<NodeFArrayBox>&  a_var,
		         const std::string&               a_string )
{
  const DisjointBoxLayout& grids( a_var.getBoxes() );
  for (DataIterator dit(grids); dit.ok(); ++dit) {
    const FArrayBox& this_var = a_var[dit].getFab();
    checkForNAN(this_var,a_string);
  }
}

void
SpaceUtils::checkForNAN( const FArrayBox&    a_var,
		         const std::string&  a_string )
{
  int ncomp = a_var.nComp();
  Real this_val;

  BoxIterator bit(a_var.box());
  for (bit.begin(); bit.ok(); ++bit) {
    IntVect iv = bit();

    for (int n=0; n<ncomp; n++) {
      this_val = a_var(iv,n);
      if(this_val!=this_val) {
        cout << "procID() = " << procID() << endl;
	cout << a_string << endl;
        cout << "NAN found at iv = " << iv << " and comp = " << n << endl;
        cout << "this_val = " << this_val << endl;
      }
    }

  }

}


#endif

#include "NamespaceFooter.H"
