#ifndef _FIELD_BC_UTILS_H_
#define _FIELD_BC_UTILS_H_

#include "Box.H"
#include "LoHiSide.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "IntVect.H"
#include "BoxIterator.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"

#include "EMFields.H"
#include "DomainGrid.H"
#include "BoundaryBoxLayout.H"
#include "BoundaryConditions.H"
#include "InsulatorBC.H"
//#include "SpaceUtils.H"

#include "NamespaceHeader.H"

namespace FieldBCUtils {

   inline void getInsulatorBox( Box&        a_box,
                          const FArrayBox&  a_binary,
                          const int&        a_bdry_dir )
   {
   
      // modify passed boundary box to only cover insulator region
    
      IntVect this_index = a_box.smallEnd();
      Real this_binary;
      for (int dir=0; dir<SpaceDim; dir++) {

         if(dir==a_bdry_dir) continue;
         
         const int box_imin = a_box.smallEnd(dir);
         const int box_imax = a_box.bigEnd(dir);

         // set new small index of box for this dir
         for (int i(box_imin); i<box_imax; i++) {
            this_index[dir] = i;
            this_binary = a_binary.get(this_index,0);
            if(this_binary>0) {
               a_box.setSmall(dir,i);
            break;
            }
         }
   
         // set new big index of box for this dir
         for (int i(box_imax); i>box_imin; i--) {
            this_index[dir] = i;
            this_binary = a_binary.get(this_index,0);
            if(this_binary>0) {
               a_box.setBig(dir,i);
               break;
            }
         }
  
      }

   }
   
   inline bool getInsulatorBox_PC( Box&        a_box,
                             const FArrayBox&  a_binary,
                             const int&        a_bdry_dir )
   {
      bool flag = false;
      // modify passed boundary box to only cover insulator region
    
      IntVect this_index = a_box.smallEnd();
      Real this_binary;
      for (int dir=0; dir<SpaceDim; dir++) {

         if(dir==a_bdry_dir) continue;
         
         const int box_imin = a_box.smallEnd(dir);
         const int box_imax = a_box.bigEnd(dir);


         // set new small index of box for this dir
         for (int i(box_imin); i<=box_imax; i++) {
            this_index[dir] = i;
            this_binary = a_binary.get(this_index,0);
            if(this_binary>0) {
               a_box.setSmall(dir,i);
               flag = true;
               break;
            }
         }
   
         // set new big index of box for this dir
         for (int i(box_imax); i>=box_imin; i--) {
            this_index[dir] = i;
            this_binary = a_binary.get(this_index,0);
            if(this_binary>0) {
               a_box.setBig(dir,i);
               flag = true;
               break;
            }
         }

      }

      return flag;
   }
   
   inline
   void setInsulatorConductorBC_Bv( LevelData<FArrayBox>&  a_dst,
                              const DomainGrid&            a_mesh,
                              const BoundaryBoxLayout&     a_bdry_layout,
                              const InsulatorBC&           a_insulator_bc,
                              const Real                   a_time )
   {
      CH_TIME("FieldBCUtils::setInsulatorConductorBC_Bv()");
         
      const LevelData<FArrayBox>& Xphys = a_mesh.getXcc();
      
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<FArrayBox>& ICbinary = a_insulator_bc.getICbinary_cc();
      const RealVect& Xmin_insulator = a_insulator_bc.getXmin();
      const int bdry_dir = a_bdry_layout.dir();

      Real bdry_val = a_insulator_bc.getTimeFunctionValue(a_time);

      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         const FArrayBox& this_ICbinary( ICbinary[interior_dit] );
         const FArrayBox& this_X( Xphys[interior_dit] );
         FArrayBox& this_dst( a_dst[interior_dit] );
        
         // grow fill box to include tranverse ghosts 
         const Box fill_box( bdry_grids[dit] );
         IntVect grow_vect = a_dst.ghostVect();
         grow_vect[bdry_dir] = 0;
         Box insulator_box = grow(fill_box,grow_vect);
           
         // reduce grown fill_box to insulator box
         getInsulatorBox(insulator_box,this_ICbinary,bdry_dir);     
               
         // set the value in the boundary for the appropriate comp
         const int TEM_comp = 0; // this should come from insulator_bc object
         IntVect ib;         // cell index
         Real local_X, local_val;
         BoxIterator bit(insulator_box); // grid indices for boundary box    
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            if(a_mesh.axisymmetric()) {
               local_X = this_X.get(ib,bdry_dir);
               local_val = bdry_val*Xmin_insulator[bdry_dir]/local_X;
            }
            else { local_val = bdry_val; }
	    for (int n=0; n<this_dst.nComp(); n++) {
               if(n==TEM_comp) { this_dst.set(ib,n,local_val); }
	       else { this_dst.set(ib,n,0.0); }
	    }
         }

      }
   }
   
   inline
   void setInsulatorConductorBC_B( LevelData<FluxBox>&  a_dst,
                             const BoundaryBoxLayout&   a_bdry_layout,
                             const InsulatorBC&         a_insulator_bc,
                             const Real                 a_time )
   {
      CH_TIME("FieldBCUtils::setInsulatorConductorBC_B()");
         
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<FluxBox>& ICbinary = a_insulator_bc.getICbinary_fc();
      
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         const int bdry_dir = a_bdry_layout.dir();
            
         for (int dir(0); dir<SpaceDim; dir++) {
            
            FArrayBox& this_dst( a_dst[a_bdry_layout.dataIndex(dit)][dir] );
            const FArrayBox& this_ICbinary( ICbinary[interior_dit][dir] );
               
            // convert box to node type in dir direction and grow box to include tranverse ghosts 
            Box fill_box_grown = surroundingNodes(fill_box,dir);
            IntVect grow_vect = a_dst.ghostVect();
            grow_vect[bdry_dir] = 0;
            fill_box_grown.grow(grow_vect);
                  
            // take care of conductor part of bcs
            string B_bc_type = "even";
            if(dir==bdry_dir) { // multiply by ICbinary to set zero on conductor
               B_bc_type = "odd";
               this_dst.mult(this_ICbinary,fill_box_grown,0,0,this_dst.nComp());

               // after mult, adjust box to not include boundary
               const int ISIDE(a_bdry_layout.side()); 
               if(ISIDE==0) fill_box_grown.growHi(bdry_dir,-1);
               if(ISIDE==1) fill_box_grown.growLo(bdry_dir,-1);
            }

            // apply bc for conductor
            BoundaryConditions::setBC( this_dst,
                                       fill_box_grown,
                                       0,
                                       B_bc_type,
                                       bdry_dir,
                                       a_bdry_layout.side() );
            
            // reduce grown fill_box to insulator box and apply extrap bc
            Box insulator_box = fill_box_grown;
            getInsulatorBox(insulator_box,this_ICbinary,bdry_dir);     
                  
            if(dir==bdry_dir) { B_bc_type = "extrapolate"; }
	    else { B_bc_type = "zero"; }
            BoundaryConditions::setBC( this_dst,
                                       insulator_box,
                                       0,
                                       B_bc_type,
                                       bdry_dir,
                                       a_bdry_layout.side() );

         } // end loop over dirs

      } // end loop over boxes

   } // end function
   
   inline
   void setInsulatorConductorBC_E( LevelData<EdgeDataBox>&  a_dst,
                             const BoundaryBoxLayout&       a_bdry_layout,
                             const InsulatorBC&             a_insulator_bc,
                             const bool                     a_conservative_wall,
                             const Real                     a_time )
   {
      CH_TIME("FieldBCUtils::setInsulatorConductorBC_E()");
         
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<EdgeDataBox>& ICbinary = a_insulator_bc.getICbinary_ec();
   
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         const int bdry_dir = a_bdry_layout.dir();
            
         for (int dir(0); dir<SpaceDim; dir++) {
            
            FArrayBox& this_dst( a_dst[a_bdry_layout.dataIndex(dit)][dir] );
            const FArrayBox& this_ICbinary( ICbinary[interior_dit][dir] );
               
            // convert box to node type in dir direction and grow box to include tranverse ghosts 
            Box fill_box_grown = surroundingNodes(fill_box);
            fill_box_grown.enclosedCells(dir);
            IntVect grow_vect = a_dst.ghostVect();
            grow_vect[bdry_dir] = 0;
            fill_box_grown.grow(grow_vect);
            
            // apply bc for conductor
         
            std::string E_bc_type;
            if(dir!=bdry_dir) {
               if(a_conservative_wall) { E_bc_type = "odd"; } 
               else { E_bc_type = "odd"; }
               this_dst.mult(this_ICbinary,fill_box_grown,0,0,this_dst.nComp());
   
               const int ISIDE(a_bdry_layout.side());
               if(ISIDE==0) { fill_box_grown.growHi(bdry_dir,-1); }
               if(ISIDE==1) { fill_box_grown.growLo(bdry_dir,-1); }
            }
            else {
               if(a_conservative_wall) { E_bc_type = "even"; }
               else { E_bc_type = "extrapolate"; }
            }
            
            BoundaryConditions::setBC( this_dst,
                                       fill_box_grown,
                                       0,
                                       E_bc_type,
                                       bdry_dir,
                                       a_bdry_layout.side() );
            
            // overwrite bc in insulator region
            
            if(dir!=bdry_dir) {
               // reduce grown fill_box to insulator box and setBC in insulator region
               Box insulator_box = fill_box_grown;
               getInsulatorBox(insulator_box,this_ICbinary,bdry_dir);     

               if(a_conservative_wall) { E_bc_type = "extrapolate"; }
               else { E_bc_type = "extrapolate"; }
               BoundaryConditions::setBC( this_dst,
                                          insulator_box,
                                          0,
                                          E_bc_type,
                                          bdry_dir,
                                          a_bdry_layout.side() );
            }

         } // end loop over dirs

      } // end loop over boxes

   } // end function
   
   inline
   void setInsulatorConductorBC_Ev( LevelData<NodeFArrayBox>&  a_dst,
                              const BoundaryBoxLayout&         a_bdry_layout,
                              const InsulatorBC&               a_insulator_bc,
                              const bool                       a_conservative_wall,
                              const Real                       a_time )
   {
      CH_TIME("FieldBCUtils::setInsulatorConductorBC_Ev()");
         
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<NodeFArrayBox>& ICbinary = a_insulator_bc.getICbinary_nc();
      
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         const FArrayBox& this_ICbinary( ICbinary[interior_dit].getFab() );
         FArrayBox& this_dst( a_dst[interior_dit].getFab() );
        
         // grow fill box to include tranverse ghosts and convert to node type
         const Box fill_box( bdry_grids[dit] );
         IntVect grow_vect = a_dst.ghostVect();
         grow_vect[a_bdry_layout.dir()] = 0;
         Box fill_box_grown = grow(fill_box,grow_vect);
         fill_box_grown.surroundingNodes();           
            
         // multiply by ICbinary to set zero on conductor
         for (auto n=0; n<this_dst.nComp(); ++n) {
            this_dst.mult(this_ICbinary,fill_box_grown,0,n,1);
         }

         // dont change value on boundary in setBC()
         const int ISIDE(a_bdry_layout.side());
         if(ISIDE==0) { fill_box_grown.growHi(a_bdry_layout.dir(),-1); }
         if(ISIDE==1) { fill_box_grown.growLo(a_bdry_layout.dir(),-1); }
         
         // apply bc for conductor
         std::string Ev_bc_type;
        
         // JRA, energy conservation requires zero BC for Epar,
         // since Jpar at ghosts does not enter field solve
         if(a_conservative_wall) { Ev_bc_type = "zero"; } // JRA, zero correct for virtual E
         else { Ev_bc_type = "odd"; }
         for (int n(0); n<this_dst.nComp(); n++) {
            BoundaryConditions::setBC( this_dst,
                                       fill_box_grown,
                                       n,
                                       Ev_bc_type,
                                       a_bdry_layout.dir(),
                                       a_bdry_layout.side() );
         }
         
         // overwrite bc in insulator region

         // reduce grown fill_box to insulator box and setBC in insulator region
         Box insulator_box = fill_box_grown;
         getInsulatorBox(insulator_box,this_ICbinary,a_bdry_layout.dir());     

         if(a_conservative_wall) { Ev_bc_type = "zero"; }
         else { Ev_bc_type = "extrapolate"; }
         for (int n(0); n<this_dst.nComp(); n++) {
            BoundaryConditions::setBC( this_dst,
                                       insulator_box,
                                       n,
                                       Ev_bc_type,
                                       a_bdry_layout.dir(),
                                       a_bdry_layout.side() );
         }

      } // end loop over boxes

   } // end function
   
   //////////////////////////////////////////////////////////////////////
   //
   //
   //
   //////////////////////////////////////////////////////////////////////

   inline
   void setCellBC( LevelData<FArrayBox>&        a_dst,
             const DomainGrid&                  a_mesh,
             const Vector<std::string>&         a_all_bc_type,
             const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc,
             const Real                         a_time )
   {
      CH_TIME("FieldBCUtils::setCellBC (virtual magnetic field)");
        
      // 1D CAR: Bv(n=0) = By,  Bv(n=1) = Bz     2D CAR: Bv(n=0) = Bz
      // 1D CYL: Bv(n=0) = Bth, Bv(n=1) = Bz     2D CYL: Bv(n=0) = Bth
      // 1D SPH: Bv(n=0) = Bth, Bv(n=1) = Bphi

      //const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
      const LevelData<FArrayBox>& Xcc = a_mesh.getXcc();
      const string& geom_type = a_mesh.geomType();

      const BoundaryBoxLayoutPtrVect& all_bdry_layouts = a_mesh.getBoundaryLayout();
      for (int b(0); b<all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const std::string this_bc_type (a_all_bc_type[b]);
         string Bv_bc_type = this_bc_type;

         if(this_bc_type=="symmetry") { Bv_bc_type = "odd"; }
         if(this_bc_type=="axis") { Bv_bc_type = "odd"; }
         if(this_bc_type=="conductor") { Bv_bc_type = "even"; }
         if(this_bc_type=="insulator_conductor") { Bv_bc_type = "even"; }
         if(this_bc_type=="extrapolate_zeroBv") { Bv_bc_type = "zero"; }

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            FArrayBox& this_dst( a_dst[interior_dit] );
        
            // grow fill box to include tranverse ghosts 
            const Box fill_box( bdry_grids[dit] );
            IntVect grow_vect = a_dst.ghostVect();
            grow_vect[bdry_layout.dir()] = 0;
            Box fill_box_grown = grow(fill_box,grow_vect);
         
            for (int n(0); n<this_dst.nComp(); n++) {
	       bool scale_by_r = false; 
               if(this_bc_type=="axis" && geom_type=="cyl_R") {
	          if(n==0) { Bv_bc_type = "odd"; }
		  else if(n==1) { Bv_bc_type = "even"; }
	       }
               if( (this_bc_type=="conductor" || this_bc_type=="insulator_conductor") 
	          && a_mesh.axisymmetric() ) {
		  if(geom_type=="cyl_R" && n==1) { scale_by_r = false; }
		  else { scale_by_r = true; }
	       }
	       if(scale_by_r) {
                  const FArrayBox& this_Xcc( Xcc[interior_dit] );
                  FArrayBox rBn;
                  rBn.define(this_dst.box(),1);
                  rBn.copy(this_dst,n,0,1);
                  rBn.mult(this_Xcc,0,0,1);
                  BoundaryConditions::setBC( rBn,
                                             fill_box_grown,
                                             0,
                                             Bv_bc_type,
                                             bdry_layout.dir(),
                                             bdry_layout.side() );
                  this_dst.copy(rBn,fill_box_grown,0,fill_box_grown,n,1);
                  this_dst.divide(this_Xcc,fill_box_grown,0,n,1);
               }
               else { 
                  BoundaryConditions::setBC( this_dst,
                                             fill_box_grown,
                                             n,
                                             Bv_bc_type,
                                             bdry_layout.dir(),
                                             bdry_layout.side() );
               }
            }

         }

         if(this_bc_type=="insulator_conductor") { // re-do bcs in insulator region
            setInsulatorConductorBC_Bv(a_dst,a_mesh,bdry_layout,*a_all_insulator_bc[b],a_time);
         }
      
      } // end loop over boundaries

   } // end function
   
   inline
   void setFluxBC( LevelData<FluxBox>&        a_dst,
             const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
             const Vector<std::string>&       a_all_bc_type,
             const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc,
             const Real                       a_time )
   {
      CH_TIME("FieldBCUtils::setFluxBC() (magnetic field)");
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         
         if(this_bc_type=="insulator_conductor") {
            setInsulatorConductorBC_B(a_dst,bdry_layout,*a_all_insulator_bc[b],a_time);
         }
         else {
      
            const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
            for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
               const Box fill_box( bdry_grids[dit] );
               
               for (int dir(0); dir<SpaceDim; dir++) {
            
                  FArrayBox& this_dst( a_dst[bdry_layout.dataIndex(dit)][dir] );

                  // convert box to node type in dir direction and grow box to include tranverse ghosts 
                  Box fill_box_grown = surroundingNodes(fill_box,dir);
                  const int bdry_dir = bdry_layout.dir();
                  IntVect grow_vect = a_dst.ghostVect();
                  grow_vect[bdry_dir] = 0;
                  fill_box_grown.grow(grow_vect);
                     
                  std::string B_bc_type = this_bc_type;
                  if(this_bc_type=="axis") {
                     if(dir==bdry_layout.dir()) { B_bc_type = "odd"; }
		     else { B_bc_type = "even"; }
                  }
                  if(this_bc_type=="symmetry") {
                     if(dir==bdry_layout.dir()) { B_bc_type = "even"; }
		     else { B_bc_type = "odd"; }
                  }
                  if(this_bc_type=="conductor") {
                     if(dir==bdry_layout.dir()) {
                        B_bc_type = "odd";
                        this_dst.setVal(0.0,fill_box_grown,0,this_dst.nComp());
                     }
		     else { B_bc_type = "even"; }
                  }
                  if(this_bc_type=="extrapolate_zeroBv") { B_bc_type = "extrapolate"; }

                  if(dir==bdry_layout.dir()) {
                     const int ISIDE(bdry_layout.side());
                     if(ISIDE==0) { fill_box_grown.growHi(bdry_layout.dir(),-1); }
                     if(ISIDE==1) { fill_box_grown.growLo(bdry_layout.dir(),-1); }
                     if(this_bc_type=="axis") { // enforce zero on axis for Br
	                Box axis_box = fill_box_grown;
                        axis_box.setBig(dir,0);
                        axis_box.setSmall(dir,0);
			this_dst.setVal(0.0,axis_box,0,1);
		     }
                  } 
                  BoundaryConditions::setBC( this_dst,
                                             fill_box_grown,
                                             0,
                                             B_bc_type,
                                             bdry_layout.dir(),
                                             bdry_layout.side() );
               
               } // end loop over dirs
         
            } // end loop over boxes

         } // end if/else insulator_conductor

      } // end loop over boundaries

   } // end function
   
   inline
   void setEdgeBC( LevelData<EdgeDataBox>&      a_dst,
             const BoundaryBoxLayoutPtrVect&    a_all_bdry_layouts,
             const Vector<std::string>&         a_all_bc_type,
             const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc,
             const bool                         a_conservative_wall,
             const Real                         a_time )
   {
      CH_TIME("FieldBCUtils::setEdgeBC() (electric field)");
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         
         if(this_bc_type=="insulator_conductor") {
            setInsulatorConductorBC_E( a_dst,bdry_layout, *a_all_insulator_bc[b],
                                       a_conservative_wall, a_time );
         }
         else {
         
            const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
            for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
               const Box fill_box( bdry_grids[dit] );
               
               for (int dir(0); dir<SpaceDim; dir++) {
            
                  FArrayBox& this_dst( a_dst[bdry_layout.dataIndex(dit)][dir] );

                  // convert box to node type in !dir direction and grow box to include tranverse ghosts 
                  Box fill_box_grown = surroundingNodes(fill_box);
                  fill_box_grown.enclosedCells(dir);
                  const int bdry_dir = bdry_layout.dir();
                  IntVect grow_vect = a_dst.ghostVect();
                  grow_vect[bdry_dir] = 0;
                  fill_box_grown.grow(grow_vect);
                     
                  std::string E_bc_type = this_bc_type;
                  if(this_bc_type=="axis" || this_bc_type=="symmetry") {
                     if(dir==bdry_dir) { E_bc_type = "odd"; }
		     else { E_bc_type = "even"; }
                  }
                  if(this_bc_type=="conductor") {
                     if(dir!=bdry_layout.dir()) {
                        if(a_conservative_wall) { E_bc_type = "odd"; }
                        else { E_bc_type = "odd"; }
                        this_dst.setVal(0.0,fill_box_grown,0,this_dst.nComp());
                     }
                     else {
                        if(a_conservative_wall) { E_bc_type = "even"; }
                        else { E_bc_type = "extrapolate"; }
                     }
                  }
                  if(this_bc_type=="extrapolate_zeroBv") E_bc_type = "extrapolate";

                  if(dir!=bdry_layout.dir()) { // dont change value on edge in setBC()
                     const int ISIDE(bdry_layout.side());
                     if(ISIDE==0) fill_box_grown.growHi(bdry_layout.dir(),-1);
                     if(ISIDE==1) fill_box_grown.growLo(bdry_layout.dir(),-1);
                  } 
                  BoundaryConditions::setBC( this_dst,
                                             fill_box_grown,
                                             0,
                                             E_bc_type,
                                             bdry_layout.dir(),
                                             bdry_layout.side() );

               } // end loop over dirs
         
            } // end loop over boxes

         } // end if/else insulator_conduct

      } // end loop over boundaries
      
   } // end function
   
   inline
   void setNodeBC( LevelData<NodeFArrayBox>&    a_dst,
             const DomainGrid&                  a_mesh,
             const Vector<std::string>&         a_all_bc_type,
             const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc,
             const bool                         a_conservative_wall,
             const Real                         a_time )
   {
      CH_TIME("FieldBCUtils::setNodeBC (virtual electric field)");

      const string& geom_type = a_mesh.geomType();
      const BoundaryBoxLayoutPtrVect& all_bdry_layouts = a_mesh.getBoundaryLayout();
      for (int b(0); b<all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         
         const std::string this_bc_type(a_all_bc_type[b]);
         if(this_bc_type=="insulator_conductor") {
            setInsulatorConductorBC_Ev( a_dst,bdry_layout, *a_all_insulator_bc[b], 
                                        a_conservative_wall, a_time );
         }
         else {

            for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         
               const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
               FArrayBox& this_dst( a_dst[interior_dit].getFab() );
        
               // convert to node type and grow box to include tranverse ghosts 
               const Box fill_box( bdry_grids[dit] );
               Box fill_box_grown = surroundingNodes(fill_box);
               IntVect grow_vect = a_dst.ghostVect();
               grow_vect[bdry_layout.dir()] = 0;
               fill_box_grown.grow(grow_vect);
      
               std::string Ev_bc_type;
               if(this_bc_type=="symmetry") { Ev_bc_type = "even"; }
               if(this_bc_type=="conductor") {
                  if(a_conservative_wall) { Ev_bc_type = "zero"; }
                  else { Ev_bc_type = "odd"; }
                  this_dst.setVal(0.0,fill_box_grown,0,this_dst.nComp());
               }
               if(this_bc_type=="extrapolate_zeroBv") Ev_bc_type = "extrapolate";
            
               // dont change value on face in setBC()
               const int ISIDE(bdry_layout.side());
               if(ISIDE==0) fill_box_grown.growHi(bdry_layout.dir(),-1);
               if(ISIDE==1) fill_box_grown.growLo(bdry_layout.dir(),-1);
            
               for (int n(0); n<this_dst.nComp(); n++) {
                  if(this_bc_type=="axis" && a_mesh.axisymmetric()) {
                     if(geom_type=="cyl_R" && n==1) { 
			Ev_bc_type = "even"; // Ez
		     }
		     else {
		        Ev_bc_type = "odd";  // Eth or Ephi
	                Box axis_box = fill_box_grown;
                        axis_box.setBig(0,0);
                        axis_box.setSmall(0,0);
		        this_dst.setVal(0.0,axis_box,n,1);
		     }
		  }
                  BoundaryConditions::setBC( this_dst,
                                             fill_box_grown,
                                             n,
                                             Ev_bc_type,
                                             bdry_layout.dir(),
                                             bdry_layout.side() );
               }

            } // end loop over boxes
 
         } // end if/else insulator_conduct
      
      } // end loop over boundaries

   } // end function
  
   //////////////////////////////////////////////////////////////////////
   //
   //
   //
   //////////////////////////////////////////////////////////////////////
 
   inline
   void setFluxBC( LevelData<FluxBox>& a_dst,
             const BoundaryBoxLayoutPtrVect& a_all_bdry_layouts,
             const LevelData<FluxBox>& a_src )
   {
      CH_TIME("FieldBCUtils::setFluxBC() from setFluxBC()");
      CH_assert(a_dst.nComp() == a_src.nComp());

      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const int bdry_dir = bdry_layout.dir();
         const int bdry_side = bdry_layout.side();

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
    
            FluxBox& this_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const FluxBox& this_src( a_src[bdry_layout.dataIndex(dit)] );
            
            if(a_src.ghostVect()==IntVect::Zero) {
               for (int dir=0; dir<SpaceDim; dir++) {
                  FArrayBox& dst_dir(this_dst[dir]);
                  const FArrayBox& src_dir(this_src[dir]);
                  if(dir==bdry_dir) {
                     Box fill_box_mod = src_dir.box();
                     if(bdry_side==0) fill_box_mod.setBig(bdry_dir,fill_box_mod.smallEnd(bdry_dir));
                     if(bdry_side==1) fill_box_mod.setSmall(bdry_dir,fill_box_mod.bigEnd(bdry_dir));
                     dst_dir.copy(src_dir, fill_box_mod);
                  }
               }
            } 
            else {
               // grow fill box to include tranverse ghosts 
               const Box& fill_box( bdry_grids[dit] );
               IntVect grow_vect = a_dst.ghostVect();
               grow_vect[bdry_dir] = 0;
               Box fill_box_grown = grow(fill_box,grow_vect);

               const Box& src_box = this_src.box();
               if(bdry_side==1) fill_box_grown.setBig(bdry_dir,src_box.bigEnd(bdry_dir));
               if(bdry_side==0) fill_box_grown.setSmall(bdry_dir,src_box.smallEnd(bdry_dir));
               this_dst.copy(this_src, fill_box_grown);
            }
         }

      }
   }

   inline
   void setEdgeBC( LevelData<EdgeDataBox>& a_dst,
             const BoundaryBoxLayoutPtrVect& a_all_bdry_layouts,
             const LevelData<EdgeDataBox>& a_src )
   {
      CH_TIME("FieldBCUtils::setEdgeBC() from setEdgeBC()");
      CH_assert(a_dst.nComp() == a_src.nComp());

      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const int bdry_dir = bdry_layout.dir();
         const int bdry_side = bdry_layout.side();

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
    
            EdgeDataBox& this_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const EdgeDataBox& this_src( a_src[bdry_layout.dataIndex(dit)] );
            
            if(a_src.ghostVect()==IntVect::Zero) {
               for (int dir=0; dir<SpaceDim; dir++) {
                  FArrayBox& dst_dir(this_dst[dir]);
                  const FArrayBox& src_dir(this_src[dir]);
                  if(dir!=bdry_dir) {
                     Box fill_box_mod = src_dir.box();
                     if(bdry_side==0) fill_box_mod.setBig(bdry_dir,fill_box_mod.smallEnd(bdry_dir));
                     if(bdry_side==1) fill_box_mod.setSmall(bdry_dir,fill_box_mod.bigEnd(bdry_dir));
                     dst_dir.copy(src_dir, fill_box_mod);
                  }
               }
            } 
            else {
               // grow fill box to include tranverse ghosts 
               const Box& fill_box( bdry_grids[dit] );
               IntVect grow_vect = a_dst.ghostVect();
               grow_vect[bdry_dir] = 0;
               Box fill_box_grown = grow(fill_box,grow_vect);

               const Box& src_box = this_src.box();
               if(bdry_side==1) fill_box_grown.setBig(bdry_dir,src_box.bigEnd(bdry_dir));
               if(bdry_side==0) fill_box_grown.setSmall(bdry_dir,src_box.smallEnd(bdry_dir));
               this_dst.copy(fill_box_grown, a_dst.interval(), fill_box_grown, this_src, a_src.interval());
            }
         }

      }
   }

   //////////////////////////////////////////////////////////////////////
   //
   //
   //
   //////////////////////////////////////////////////////////////////////

   inline
   void applyToJ_PIC( LevelData<EdgeDataBox>&    a_J_inPlane,
                      LevelData<NodeFArrayBox>&  a_J_virtual,
                const DomainGrid&                a_mesh,
                const Vector<std::string>&       a_all_bc_type,
                const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc )
   {
      CH_TIME("FieldBCUtils::applyToJ_PIC()");

      // Comment on J normal to boundary. For Symmetry, we subtract ghost Jnorm
      // from mirror location across boundary, which is the same as depositing
      // the mirror current to the interior. For other boundaries, we add ghost
      // Jnorm, which amount so using NGP for Jnorm close to boundary and is
      // needed to maintain exact energy conservation (only verified for insulator
      // BC so far)
      
      const string& geom_type = a_mesh.geomType();

      // loop over non-periodic boundaries and apply BCs to J
      const BoundaryBoxLayoutPtrVect& bdry_layout = a_mesh.getBoundaryLayout();
      for (int b(0); b<bdry_layout.size(); b++) {

         const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
         const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
         const int bdry_dir = this_bdry_layout.dir();
         const int bdry_side(this_bdry_layout.side());
         std::string this_bc = a_all_bc_type[b];
         
         for(DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

            const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
            const Box bdry_box( bdry_grids[dit] );
            for (int dir(0); dir<SpaceDim; dir++) {
          
               // convert cell bdry box to a edge bdry box
               Box edge_box = surroundingNodes(bdry_box);
               edge_box.enclosedCells(dir);

	       // grow ghost cells in transverse direction
               IntVect grow_vect = a_J_inPlane.ghostVect();
               grow_vect[bdry_dir] = 0;
               edge_box.grow(grow_vect);

               // collapse edge_box to 1 cell thick in bdry_dir direction                  
               if(bdry_side==0) { edge_box.setSmall(bdry_dir,edge_box.bigEnd(bdry_dir)); }
               if(bdry_side==1) { edge_box.setBig(bdry_dir,edge_box.smallEnd(bdry_dir)); }
                  
               FArrayBox& this_J( a_J_inPlane[interior_dit][dir] );
               if(dir==bdry_dir) {
                  Box dst_box = edge_box;
                  Box src_box = edge_box;
                  const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
                  for (int n=0; n<nG; n++) {
                     if(bdry_side==0) { dst_box.shift(bdry_dir,1); }
                     if(bdry_side==1) { dst_box.shift(bdry_dir,-1); }
                     if(this_bc=="symmetry" || this_bc=="axis") {
                        this_J.minus(this_J,src_box,dst_box,0,0,1);
                     }
                     else {
                        this_J.plus(this_J,src_box,dst_box,0,0,1);
                     }
                     if(bdry_side==0) { src_box.shift(bdry_dir,-1); }
                     if(bdry_side==1) { src_box.shift(bdry_dir,1); }
                  }

               }
               else {
                  if(this_bc=="symmetry" || this_bc=="axis") { 
                     Box dst_box = edge_box;
                     Box src_box = edge_box;
                     //this_J.plus(this_J,src_box,dst_box,0,0,1); // 2X on bdry (done via corrected Jacobian)
                     const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
                     for (int n=0; n<nG; n++) {
                        dst_box.shift(bdry_dir,1-2*bdry_side);
                        src_box.shift(bdry_dir,2*bdry_side-1);
                        this_J.plus(this_J,src_box,dst_box,0,0,1);
                     }
                  }
                  else if(this_bc=="conductor") {
                     Box dst_box = edge_box;
                     Box src_box = edge_box;
                     this_J.minus(this_J,src_box,dst_box,0,0,1); // set zero right on conductor
                     const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
                     for (int n=0; n<nG; n++) {
                        dst_box.shift(bdry_dir,1-2*bdry_side);
                        src_box.shift(bdry_dir,2*bdry_side-1);
                        this_J.minus(this_J,src_box,dst_box,0,0,1);
                     }
		  }
                  else if(this_bc=="insulator_conductor") {
            
                     // For CC1, only have J in first ghost cell past boundary
	             // Add 2*Jpara(-1) to Jpara(0) and subtract Jpara(-1) from
		     // Jpara(1). This gives exact charge conservation for insulator.
	             // Have not yet thought about conductor
                     Box dst_box = edge_box;
                     Box src_box = edge_box;
                     if(bdry_side==0) { src_box.shift(bdry_dir,-1); }
                     if(bdry_side==1) { src_box.shift(bdry_dir,1); }
		     const Real scale = 2.0;
                     this_J.plus(this_J,src_box,dst_box,scale,0,0,1);
                     if(bdry_side==0) { dst_box.shift(bdry_dir,1); }
                     if(bdry_side==1) { dst_box.shift(bdry_dir,-1); }
                     this_J.minus(this_J,src_box,dst_box,0,0,1);
		     
		     const InsulatorBC& this_insulator_bc = *a_all_insulator_bc[b];
                     const LevelData<EdgeDataBox>& ICbinary = this_insulator_bc.getICbinary_ec();
                     const FArrayBox& this_ICbinary( ICbinary[interior_dit][dir] );

	             // set zero right on conductor boundary
                     this_J.mult(this_ICbinary,edge_box,0,0,1);

                  }
               }

            }
           
#if CH_SPACEDIM<3
            
            // collapse node_box to 1 cell thick in bdry_dir direction (on boundary)
            Box node_box = surroundingNodes(bdry_box);
	    
	    // grow ghost cells in transverse direction
            IntVect grow_vect = a_J_virtual.ghostVect();
            grow_vect[bdry_dir] = 0;
            node_box.grow(grow_vect);

            if(bdry_side==0) node_box.setSmall(bdry_dir,node_box.bigEnd(bdry_dir));
            if(bdry_side==1) node_box.setBig(bdry_dir,node_box.smallEnd(bdry_dir));

            FArrayBox& this_J( a_J_virtual[interior_dit].getFab() );
            if(this_bc=="symmetry" || this_bc=="axis") { 
     	       if(this_bc=="axis") { // enforce Jth = 0 right on axis
		                     // 1 + 1 = 2 factor for Jz (done via corrected Jacobian)
				     // 1 - 1 = 0 factor for Jth (enforced here)
	          Box axis_box = node_box;
                  axis_box.setBig(0,0);
                  axis_box.setSmall(0,0);
	          this_J.setVal(0.0,axis_box,0,1); // Jth
                  if(geom_type=="sph_R") { this_J.setVal(0.0,axis_box,1,1); } // Jphi 
	       }
     
	       Box dst_box = node_box;
               Box src_box = node_box;
               const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
               for (int n=0; n<nG; n++) {
                  dst_box.shift(bdry_dir,1-2*bdry_side);
                  src_box.shift(bdry_dir,2*bdry_side-1);
                  if(this_bc=="axis") {
#if CH_SPACEDIM==1
                     if(geom_type=="sph_R") { 
                        this_J.minus(this_J,src_box,dst_box,0,0,1); // Jth
			this_J.minus(this_J,src_box,dst_box,1,1,1); // Jphi
		     }
		     else { 
                        this_J.minus(this_J,src_box,dst_box,0,0,1); // Jth
			this_J.plus(this_J,src_box,dst_box,1,1,1);  // Jz
		     }
#endif
                  }
                  else this_J.plus(this_J,src_box,dst_box,0,0,this_J.nComp());
               }
	    }
#endif
            
         } // end loop over boxes

      } // end loop over boundaries

   } // end applyToJ_PIC()
   
   //////////////////////////////////////////////////////////////////////
   //
   //
   //
   //////////////////////////////////////////////////////////////////////

   inline
   void setCellPCMask( LevelData<FArrayBox>&                a_dst,
                 const DomainGrid&                          a_mesh,
                 const Vector<std::string>&                 a_all_bc_type,
                 const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc,
                 const Real                                 a_time )
   {
      CH_TIME("FieldBCUtils::setCellPCMask (virtual magnetic field)");
         
      const BoundaryBoxLayoutPtrVect& all_bdry_layouts = a_mesh.getBoundaryLayout();
      for (int b(0); b<all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const std::string this_bc_type (a_all_bc_type[b]);

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const int bdry_dir = bdry_layout.dir();
            const int ISIDE(bdry_layout.side());
            FArrayBox& this_dst( a_dst[bdry_layout.dataIndex(dit)] );
        
            if(this_bc_type=="symmetry") {
              // get the boundary layer of faces in dir direction
              Box bdry_cells = bdry_grids[dit];
              // shift layer by 1 cell to be the near-boundary interior cells
              bdry_cells.shift(bdry_dir,1-2*ISIDE);
              this_dst.setVal(2.0, bdry_cells, 2*bdry_dir+(1-ISIDE));
            }

            if(this_bc_type=="axis") {
              // yet to implement
            }

            if(this_bc_type=="conductor") {
              // get the boundary layer of faces in dir direction
              Box bdry_cells = bdry_grids[dit];
              // shift layer by 1 cell to be the near-boundary interior cells
              bdry_cells.shift(bdry_dir,1-2*ISIDE);
              this_dst.setVal(0.0, bdry_cells, 2*bdry_dir+(1-ISIDE));
            }

            if(this_bc_type=="insulator_conductor") {
              {
                // get the boundary layer of faces in dir direction
                Box bdry_cells = bdry_grids[dit];
                // shift layer by 1 cell to be the near-boundary interior cells
                bdry_cells.shift(bdry_dir,1-2*ISIDE);
                this_dst.setVal(0.0, bdry_cells, 2*bdry_dir+(1-ISIDE));
              }
              {
                // now, get the insulator portion of the ghost layer
                Box bdry_cells = bdry_grids[dit];
                getInsulatorBox_PC( bdry_cells,
                                    a_all_insulator_bc[b]
                                      ->getICbinary_cc()[bdry_layout.dataIndex(dit)],
                                    bdry_dir );
                // shift layer by 1 cell to be the near-boundary interior cells
                bdry_cells.shift(bdry_dir,1-2*ISIDE);
                this_dst.setVal(1.0, bdry_cells, 2*bdry_dir+(1-ISIDE));
              }
            }

            if(this_bc_type=="extrapolate_zeroBv") {
              // yet to implement
            }

         }

         if(this_bc_type=="insulator_conductor") { // re-do bcs in insulator region
            //setInsulatorConductorBC_Bv(a_dst,a_mesh,bdry_layout,*a_all_insulator_bc[b],a_time);
         }
      } // end loop over boundaries

   } // end function
   
   inline
   void setFluxPCMask( LevelData<FluxBox>&                  a_dst,
                 const BoundaryBoxLayoutPtrVect&            a_all_bdry_layouts,
                 const Vector<std::string>&                 a_all_bc_type,
                 const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc,
                 const Real                                 a_time )
   {
      CH_TIME("FieldBCUtils::setFluxPCMask() (magnetic field)");
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         
            for (int dir(0); dir<SpaceDim; dir++) {
         
               const int bdry_dir = bdry_layout.dir();
               const int ISIDE(bdry_layout.side());
               FArrayBox& this_dst( a_dst[bdry_layout.dataIndex(dit)][dir] );

               if(this_bc_type=="axis") {
                 // yet to implement
               }

               if(this_bc_type=="symmetry") {

                 if(dir == bdry_layout.dir()) {

                   // doesn't participate in field solve;
                   // nothing to do!

                 } else {
                   // get the boundary layer of faces in dir direction
                   Box bdry_cells = surroundingNodes(bdry_grids[dit],dir);
                   // shift layer by 1 cell to be the near-boundary interior cells
                   bdry_cells.shift(bdry_dir,1-2*ISIDE);
                   this_dst.setVal(2.0, bdry_cells, (1-ISIDE));
                 }

               }

               if(this_bc_type=="conductor" || this_bc_type=="axis") {
                 if (dir != bdry_dir) {
                   // get the boundary layer of faces in dir direction
                   Box bdry_cells = surroundingNodes(bdry_grids[dit],dir);
                   // shift layer by 1 cell to be the near-boundary interior cells
                   bdry_cells.shift(bdry_dir,1-2*ISIDE);
                   this_dst.setVal(0.0, bdry_cells, (1-ISIDE));
                 } else {
                   // get the boundary layer of faces in dir direction
                   Box bdry_cells = bdry_grids[dit];
                   // shift layer by 1 cell to be the near-boundary interior cells
                   bdry_cells.shiftHalf(bdry_dir,1-2*ISIDE);
                   // include one ghost face along dir to avoid corner issues
                   if (SpaceDim == 2) {
                     bdry_cells.grow(1-dir,1);
                   }
                   this_dst.setVal(0.0, bdry_cells, 0);
                   this_dst.setVal(0.0, bdry_cells, 1);
                 }
               }

               if(this_bc_type=="insulator_conductor") {
                 if (dir != bdry_dir) {
                   // get the boundary layer of faces in dir direction
                   Box bdry_cells = surroundingNodes(bdry_grids[dit],dir);
                   // shift layer by 1 cell to be the near-boundary interior cells
                   bdry_cells.shift(bdry_dir,1-2*ISIDE);
                   this_dst.setVal(0.0, bdry_cells, (1-ISIDE));
                 } else {
                   // get the boundary layer of faces in dir direction
                   Box bdry_cells = bdry_grids[dit];
                   // shift layer by 1 cell to be the near-boundary interior cells
                   bdry_cells.shiftHalf(bdry_dir,1-2*ISIDE);
                   // include one ghost face along dir to avoid corner issues
                   if (SpaceDim == 2) {
                     bdry_cells.grow(1-dir,1);
                   }
                   this_dst.setVal(0.0, bdry_cells, 0);
                   this_dst.setVal(0.0, bdry_cells, 1);
                 }
                 if (dir != bdry_dir) {
                   // now, get the insulator portion of the ghost layer
                   Box bdry_cells = surroundingNodes(bdry_grids[dit],dir);
                   bool flag = 
                      getInsulatorBox_PC( bdry_cells,
                                          a_all_insulator_bc[b]
                                            ->getICbinary_fc()[bdry_layout.dataIndex(dit)][dir],
                                          bdry_dir );
                   if (flag) {
                     this_dst.setVal(1*(1-2*ISIDE), bdry_cells, this_dst.nComp()-1);
                     // shift layer by 1 cell to be the near-boundary interior cells
                     bdry_cells.shift(bdry_dir,1-2*ISIDE);
                     this_dst.setVal(1.0, bdry_cells, (1-ISIDE));
                     this_dst.setVal(1*(1-2*ISIDE), bdry_cells, this_dst.nComp()-2);
                   }
                 } else {
                   // now, get the insulator portion of the ghost layer
                   Box bdry_cells = bdry_grids[dit];
                   bool flag = 
                     getInsulatorBox_PC( bdry_cells,
                                         a_all_insulator_bc[b]
                                           ->getICbinary_fc()[bdry_layout.dataIndex(dit)][dir],
                                         bdry_dir );
                   if (flag) {
                     // shift layer by 1 cell to be the near-boundary interior cells
                     bdry_cells.shiftHalf(bdry_dir,1-2*ISIDE);
                     this_dst.setVal(1.0, bdry_cells, 0);
                     this_dst.setVal(1.0, bdry_cells, 1);
                   }
                 }
               }

               if(this_bc_type=="extrapolate_zeroBv") {
                 // yet to implement
               }

            } // end loop over dirs
         } // end loop over boxes
      } // end loop over boundaries

   } // end function
   
   inline
   void applyToEforDiv( LevelData<EdgeDataBox>&    a_dst,
                  const DomainGrid&                a_mesh,
                  const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                  const Vector<std::string>&       a_all_bc_type )
   {
      CH_TIME("FieldBCUtils::applyToEforDiv()");
      
      const LevelData<EdgeDataBox>& Jec = a_mesh.getJec();
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {

         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string bc_type (a_all_bc_type[b]);
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
            
	 const int bdry_dir = bdry_layout.dir();

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            FArrayBox& this_dst( a_dst[interior_dit][bdry_dir] );

            if(bc_type=="conductor" || bc_type=="insulator_conductor") {

	       // set ghost values in bdry_dir such that divergence derivative
	       // term in bdry_dir reduces to one-sided derivative.
	       // e.g., upper dir X bdry: (Ex(N)-Ex(N-1))/dX ==> (0-Ex(N-1))/(dX/2)
	       const string dst_bc_type = "odd";
               Box bdry_cells = surroundingNodes(bdry_grids[dit]);
               bdry_cells.enclosedCells(bdry_dir);
	       
	       bool scale_by_Ja = false;
               if( a_mesh.axisymmetric() && bdry_dir==0 ) { scale_by_Ja = true; }
	         
	       if(scale_by_Ja) {
                  const FArrayBox& this_Ja( Jec[interior_dit][bdry_dir] );
                  FArrayBox JaEr;
                  JaEr.define(this_dst.box(),1);
                  JaEr.copy(this_dst,0,0,1);
                  JaEr.mult(this_Ja,0,0,1);
                  BoundaryConditions::setBC( JaEr,
                                             bdry_cells,
                                             0,
                                             dst_bc_type,
                                             bdry_dir,
                                             bdry_layout.side() );
                  this_dst.copy(JaEr,bdry_cells,0,bdry_cells,0,1);
                  this_dst.divide(this_Ja,bdry_cells,0,0,1);
               }
               else { 
                  BoundaryConditions::setBC( this_dst,
                                             bdry_cells,
                                             0,
                                             dst_bc_type,
                                             bdry_dir,
                                             bdry_layout.side() );
               }

            }

	 } // end loop over boxes

      } // end loop over boundaries

   }
   
   inline
   void setEdgePCMask( LevelData<EdgeDataBox>&              a_dst,
                 const BoundaryBoxLayoutPtrVect&            a_all_bdry_layouts,
                 const Vector<std::string>&                 a_all_bc_type,
                 const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc,
                 const bool                                 a_conservative_wall,
                 const Real                                 a_time )
   {
      CH_TIME("FieldBCUtils::setEdgePCMask() (electric field)");
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {

         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         
            const int bdry_dir = bdry_layout.dir();
            const int ISIDE(bdry_layout.side());
            
            for (int dir(0); dir<SpaceDim; dir++) {
         
               FArrayBox& this_dst( a_dst[bdry_layout.dataIndex(dit)][dir] );

               if(this_bc_type=="axis") {
                 // yet to implement
               }

               if(this_bc_type=="symmetry") {
                 // doesn't participate in field solve;
                 // nothing to do!
               }

               if(this_bc_type=="conductor") {

                  if(dir != bdry_dir) {
                    // get the boundary layer of nodes in the dir direction
                    Box bdry_cells = surroundingNodes(bdry_grids[dit]).enclosedCells(dir);
                    // remove the ghost layer of nodes
                    bdry_cells.growDir(bdry_dir, (ISIDE?Side::Hi:Side::Lo), -1);
                    this_dst.setVal(0.0, bdry_cells, ISIDE);
                  } else {
                    // don't do anything?
                  }

               }

               if (this_bc_type=="insulator_conductor") {
                 if(dir != bdry_dir) {
                   // get the boundary layer of nodes in the dir direction
                   Box bdry_cells = surroundingNodes(bdry_grids[dit]).enclosedCells(dir);
                   // remove the ghost layer of nodes
                   bdry_cells.growDir(bdry_dir, (ISIDE?Side::Hi:Side::Lo), -1);
                   this_dst.setVal(0.0, bdry_cells, ISIDE);
                 } else {
                   // don't do anything?
                 }
                 if(dir != bdry_dir) {
                   // now, get the insulator portion of the ghost layer
                   Box bdry_cells = surroundingNodes(bdry_grids[dit]).enclosedCells(dir);
                   bool flag = 
                     getInsulatorBox_PC(  bdry_cells,
                                          a_all_insulator_bc[b]
                                            ->getICbinary_ec()[bdry_layout.dataIndex(dit)][dir],
                                          bdry_dir );
                   if (flag) {
                     // remove the ghost layer of nodes
                     bdry_cells.growDir(bdry_dir, (ISIDE?Side::Hi:Side::Lo), -1);
                     this_dst.setVal(1.0, bdry_cells, ISIDE);
                   }
                 } else {
                   // don't do anything?
                 }
               }

               if(this_bc_type=="extrapolate_zeroBv") {
                 // yet to implement
               }

            } // end loop over dirs
         } // end loop over boxes
      } // end loop over boundaries

   } // end function
   
   inline
   void setNodePCMask( LevelData<NodeFArrayBox>&            a_dst,
                 const DomainGrid&                          a_mesh,
                 const Vector<std::string>&                 a_all_bc_type,
                 const Vector<RefCountedPtr<InsulatorBC>>&  a_all_insulator_bc,
                 const bool                                 a_conservative_wall,
                 const Real                                 a_time )
   {
      CH_TIME("FieldBCUtils::setNodePCMask (virtual electric field)");

      const BoundaryBoxLayoutPtrVect& all_bdry_layouts = a_mesh.getBoundaryLayout();
      for (int b(0); b<all_bdry_layouts.size(); b++) {

         const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const std::string this_bc_type(a_all_bc_type[b]);

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         
            const int bdry_dir = bdry_layout.dir();
            const int ISIDE(bdry_layout.side());
            FArrayBox& this_dst( a_dst[bdry_layout.dataIndex(dit)].getFab() );

            if(this_bc_type=="symmetry") {
              // doesn't participate in field solve;
              // nothing to do!
            }

            if(this_bc_type=="conductor" || this_bc_type=="axis") {
              // get the boundary layer of nodes in the dir direction
              Box bdry_cells = surroundingNodes(bdry_grids[dit]);
              // remove the ghost layer of nodes
              bdry_cells.growDir(bdry_dir, (ISIDE?Side::Hi:Side::Lo), -1);
              for (int dir(0); dir<SpaceDim; dir++) {
                if (dir == bdry_dir) this_dst.setVal(0.0, bdry_cells, 2*dir+ISIDE);
                else {
                  this_dst.setVal(0.0, bdry_cells, 2*dir);
                  this_dst.setVal(0.0, bdry_cells, 2*dir+1);
                }
              }
            }

            if (this_bc_type=="insulator_conductor") {
              {
                // get the boundary layer of nodes in the dir direction
                Box bdry_cells = surroundingNodes(bdry_grids[dit]);
                // remove the ghost layer of nodes
                bdry_cells.growDir(bdry_dir, (ISIDE?Side::Hi:Side::Lo), -1);
                for (int dir(0); dir<SpaceDim; dir++) {
                  if (dir == bdry_dir) this_dst.setVal(0.0, bdry_cells, 2*dir+ISIDE);
                  else {
                    this_dst.setVal(0.0, bdry_cells, 2*dir);
                    this_dst.setVal(0.0, bdry_cells, 2*dir+1);
                  }
                }
              }
              {
                // now, get the insulator portion of the ghost layer
                Box bdry_cells = surroundingNodes(bdry_grids[dit]);
                bool flag = 
                  getInsulatorBox_PC(  bdry_cells,
                                       a_all_insulator_bc[b]
                                         ->getICbinary_nc()[bdry_layout.dataIndex(dit)].getFab(),
                                       bdry_dir );
                if (flag) {
                  // remove the ghost layer of nodes
                  bdry_cells.growDir(bdry_dir, (ISIDE?Side::Hi:Side::Lo), -1);
                  for (int dir(0); dir<SpaceDim; dir++) {
                    if (dir == bdry_dir) this_dst.setVal(1.0, bdry_cells, 2*dir+ISIDE);
                    else {
                      this_dst.setVal(1.0, bdry_cells, 2*dir);
                      this_dst.setVal(1.0, bdry_cells, 2*dir+1);
                    }
                  }
                }
              }
            }

            if(this_bc_type=="extrapolate_zeroBv") {
              // yet to implement
            }
         
         }
      } // end loop over boundaries

   } // end function
   
#if CH_SPACEDIM==1
   inline
   void setPhiBC( LevelData<NodeFArrayBox>&  a_phi,
	    const LevelData<NodeFArrayBox>&  a_rho,
            const DomainGrid&                a_mesh,
            const Real                       a_rhoNorm,
            const Vector<std::string>&       a_all_phi_bc_type )
   {
      CH_TIME("FieldBCUtils::setPhiBC");

      // this is a kluge for getting the potential correction field at the boundaries
      // correct for dirichlet and neumann BCs. Chombo's node-centered potential solver
      // does not properly set the ghost cells when using dirichlet BC such that the 
      // nabla^2(phi) = rho on the boundary. Furthemore, a neumann BC is applied at cell
      // center between the cell on the boundary and the next interior cell. 
      // Note that this kluge only works for 1D. For 2D, can only do periodic BCs.
      // Will eventually use a different solver for the potential.
      CH_assert(SpaceDim==1);

      const RealVect& dX(a_mesh.getdX());
      const BoundaryBoxLayoutPtrVect& all_bdry_layouts = a_mesh.getBoundaryLayout();
      for (int b(0); b<all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const int bdry_side = bdry_layout.side();
         const int bdry_dir = bdry_layout.dir();

         const std::string this_bc_type(a_all_phi_bc_type[b]);

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            FArrayBox& this_phi( a_phi[interior_dit].getFab() );
            const FArrayBox& this_rho( a_rho[interior_dit].getFab() );
        
            // convert to node type and grow box to include tranverse ghosts 
            const Box fill_box( bdry_grids[dit] );
            Box node_box = surroundingNodes(fill_box);
            IntVect grow_vect = a_phi.ghostVect();
            grow_vect[bdry_layout.dir()] = 0;
            node_box.grow(grow_vect);
               
	    // collapse node_box to 1 cell thick in bdry_dir direction                  
            if(bdry_side==0) node_box.setSmall(bdry_dir,node_box.bigEnd(bdry_dir));
            if(bdry_side==1) node_box.setBig(bdry_dir,node_box.smallEnd(bdry_dir));

            IntVect ib;         // cell index
            IntVect ibm = IntVect::Zero; ibm[bdry_dir]--;
            IntVect ibp = IntVect::Zero; ibp[bdry_dir]++;
            Real rho0, phi0, phi0p, phi0m;
            BoxIterator bit(node_box);
            if( this_bc_type == "dirichlet" ) { 
               for(bit.begin(); bit.ok(); ++bit) {
                  ib = bit();
                  rho0 = this_rho.get(ib,0)*a_rhoNorm;
                  phi0 = this_phi.get(ib,0);
		  if(bdry_side==0) {
                     phi0p = this_phi.get(ib+ibp,0);
		     phi0m = 2.0*phi0 - phi0p - rho0*dX[bdry_dir]*dX[bdry_dir];
                     this_phi.set(ib+ibm,0,phi0m);
		  }
		  else {
                     phi0m = this_phi.get(ib+ibm,0);
		     phi0p = 2.0*phi0 - phi0m - rho0*dX[bdry_dir]*dX[bdry_dir];
                     this_phi.set(ib+ibp,0,phi0p);
		  }
               }
            }
	    else if( this_bc_type == "neumann" ) {
               for(bit.begin(); bit.ok(); ++bit) {
                  ib = bit();
                  rho0 = this_rho.get(ib,0)*a_rhoNorm;
		  if(bdry_side==0) {
                     phi0p = this_phi.get(ib+ibp,0);
		     phi0m = phi0p;
		     this_phi.set(ibm,0,phi0m);
		     phi0 = phi0p + rho0*dX[bdry_dir]*dX[bdry_dir]/2.0;
                     this_phi.set(ib,0,phi0);
		  }
		  else {
                     phi0m = this_phi.get(ib+ibm,0);
		     phi0p = phi0m;
		     this_phi.set(ib+ibp,0,phi0p);
		     phi0 = phi0m + rho0*dX[bdry_dir]*dX[bdry_dir]/2.0;
                     this_phi.set(ib,0,phi0);
		  }
	       }
	    }

	    // properly fill ghost cells if more than 1
	    if(a_phi.ghostVect()[bdry_dir]>1) {
               Box fill_box2 = surroundingNodes(fill_box);
               fill_box2.grow(grow_vect);
               if(bdry_side==0) { fill_box2.growHi(bdry_dir,-1); }
               if(bdry_side==1) { fill_box2.growLo(bdry_dir,-1); }
               string dummy_bc_type = "extrapolate";
	       if(this_bc_type=="neumann") dummy_bc_type = "even";
	       else { // shrink by 1 again to keep nabla^2(phi) = -rho at bdry
                 if(bdry_side==0) { fill_box2.growHi(bdry_dir,-1); }
                 if(bdry_side==1) { fill_box2.growLo(bdry_dir,-1); }
	       }
               BoundaryConditions::setBC( this_phi,
                                          fill_box2,
                                          0,
                                          dummy_bc_type,
                                          bdry_dir,
                                          bdry_layout.side() );
	    } 
            
         } // end loop over boxes
 
      } // end loop over boundaries

   } // end function
#endif

}
#include "NamespaceFooter.H"

#endif
