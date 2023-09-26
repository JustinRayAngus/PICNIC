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

#include "ElectroMagneticFields.H"
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
            else local_val = bdry_val;
	    for (int n=0; n<this_dst.nComp(); n++) {
               if(n==TEM_comp) this_dst.set(ib,n,local_val);
	       else this_dst.set(ib,n,0.0);
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
            const Box fill_box( bdry_grids[dit] );
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
                  
            B_bc_type = "extrapolate";
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
            const Box fill_box( bdry_grids[dit] );
            Box fill_box_grown = surroundingNodes(fill_box);
            fill_box_grown.enclosedCells(dir);
            IntVect grow_vect = a_dst.ghostVect();
            grow_vect[bdry_dir] = 0;
            fill_box_grown.grow(grow_vect);
            
            // apply bc for conductor
         
            std::string E_bc_type;
            if(dir!=bdry_dir) {
               // JRA, energy conservation requires zero BC for Epar,
               // since Jpar at ghosts does not enter field solve
               if(a_conservative_wall) E_bc_type = "zero";
               else E_bc_type = "odd";
               this_dst.mult(this_ICbinary,fill_box_grown,0,0,this_dst.nComp());
   
               const int ISIDE(a_bdry_layout.side());
               if(ISIDE==0) fill_box_grown.growHi(bdry_dir,-1);
               if(ISIDE==1) fill_box_grown.growLo(bdry_dir,-1);
            }
            else {
               // JRA, energy conservation requires even BC for Eperp, 
               // which is like using NGP for Eperp near the bdry
               if(a_conservative_wall) E_bc_type = "even";
               else E_bc_type = "extrapolate";
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

               // JRA, energy conservation requires zero BC for Epar,
               // since Jpar at ghosts does not enter field solve
               if(a_conservative_wall) E_bc_type = "zero";
               else E_bc_type = "extrapolate";
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
         if(ISIDE==0) fill_box_grown.growHi(a_bdry_layout.dir(),-1);
         if(ISIDE==1) fill_box_grown.growLo(a_bdry_layout.dir(),-1);
         
         // apply bc for conductor
         std::string Ev_bc_type;
        
         // JRA, energy conservation requires zero BC for Epar,
         // since Jpar at ghosts does not enter field solve
         if(a_conservative_wall) Ev_bc_type = "zero";
         else Ev_bc_type = "odd";
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

         if(a_conservative_wall) Ev_bc_type = "zero";
         else Ev_bc_type = "extrapolate";
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
         
      const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
      const string& geom_type = a_mesh.geomType();

      const BoundaryBoxLayoutPtrVect& all_bdry_layouts = a_mesh.getBoundaryLayout();
      for (int b(0); b<all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const std::string this_bc_type (a_all_bc_type[b]);
         string Bv_bc_type = this_bc_type;

         if(this_bc_type=="symmetry") Bv_bc_type = "odd";
         if(this_bc_type=="axis") Bv_bc_type = "odd";
         if(this_bc_type=="conductor") Bv_bc_type = "even";
         if(this_bc_type=="insulator_conductor") Bv_bc_type = "even";
         if(this_bc_type=="extrapolate_zeroBv") Bv_bc_type = "zero";

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            FArrayBox& this_dst( a_dst[interior_dit] );
        
            // grow fill box to include tranverse ghosts 
            const Box fill_box( bdry_grids[dit] );
            IntVect grow_vect = a_dst.ghostVect();
            grow_vect[bdry_layout.dir()] = 0;
            Box fill_box_grown = grow(fill_box,grow_vect);
         
            for (int n(0); n<this_dst.nComp(); n++) {   
               if(this_bc_type=="axis" && n==1 && geom_type=="cyl_R") Bv_bc_type = "even";
               if(this_bc_type=="conductor" && n==0 && a_mesh.axisymmetric()) {
                  const FArrayBox& this_Ja( Jacobian[interior_dit] );
                  FArrayBox JaBy;
                  JaBy.define(this_dst.box(),0);
                  JaBy.copy(this_dst,n,0,1);
                  JaBy.mult(this_Ja,0,0,1);
                  BoundaryConditions::setBC( JaBy,
                                             fill_box_grown,
                                             n,
                                             Bv_bc_type,
                                             bdry_layout.dir(),
                                             bdry_layout.side() );
                  this_dst.copy(JaBy,fill_box_grown,0,fill_box_grown,n,1);
                  this_dst.divide(this_Ja,fill_box_grown,0,n,1);
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
                     B_bc_type = "even";
                     if(dir==bdry_layout.dir()) B_bc_type = "odd";
                  }
                  if(this_bc_type=="symmetry") {
                     B_bc_type = "odd";
                     if(dir==bdry_layout.dir()) B_bc_type = "even";
                  }
                  if(this_bc_type=="conductor") {
                     B_bc_type = "even";
                     if(dir==bdry_layout.dir()) {
                        B_bc_type = "odd";
                        this_dst.setVal(0.0,fill_box_grown,0,this_dst.nComp());
                     }
                  }
                  if(this_bc_type=="extrapolate_zeroBv") B_bc_type = "extrapolate";

                  if(dir==bdry_layout.dir()) { // dont change value on face in setBC()
                     const int ISIDE(bdry_layout.side());
                     if(ISIDE==0) fill_box_grown.growHi(bdry_layout.dir(),-1);
                     if(ISIDE==1) fill_box_grown.growLo(bdry_layout.dir(),-1);
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
                  if(this_bc_type=="axis") {
                     E_bc_type = "even";
                     if(dir==bdry_layout.dir()) E_bc_type = "odd";
                  }
                  if(this_bc_type=="symmetry") {
                     E_bc_type = "odd";
                     if(dir!=bdry_layout.dir()) E_bc_type = "even";
                  }
                  if(this_bc_type=="conductor") {
                     if(dir!=bdry_layout.dir()) {
                        // JRA, energy conservation requires zero BC for Epar,
                        // since Jpar at ghosts does not enter field solve
                        if(a_conservative_wall) E_bc_type = "zero";
                        else E_bc_type = "odd";
                        this_dst.setVal(0.0,fill_box_grown,0,this_dst.nComp());
                     }
                     else {
                        // JRA, energy conservation requires even BC for Eperp, 
                        // which is like using NGP for Eperp near the bdry
                        if(a_conservative_wall) E_bc_type = "even";
                        else E_bc_type = "extrapolate";
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
               if(this_bc_type=="symmetry") Ev_bc_type = "even";
               if(this_bc_type=="conductor") {
                  // JRA, energy conservation requires zero BC for Epar,
                  // since Jpar at ghosts does not enter field solve
                  if(a_conservative_wall) Ev_bc_type = "zero";
                  else Ev_bc_type = "odd";
                  this_dst.setVal(0.0,fill_box_grown,0,this_dst.nComp());
               }
               if(this_bc_type=="extrapolate_zeroBv") Ev_bc_type = "extrapolate";
            
               // dont change value on face in setBC()
               const int ISIDE(bdry_layout.side());
               if(ISIDE==0) fill_box_grown.growHi(bdry_layout.dir(),-1);
               if(ISIDE==1) fill_box_grown.growLo(bdry_layout.dir(),-1);
            
               for (int n(0); n<this_dst.nComp(); n++) {
                  if(this_bc_type=="axis" && n==1) Ev_bc_type = "even"; // Ez
                  if(this_bc_type=="axis" && n==0) {
		     Ev_bc_type = "odd";  // Eth
	             Box axis_box = fill_box_grown;
                     axis_box.setBig(0,0);
                     axis_box.setSmall(0,0);
		     this_dst.setVal(0.0,axis_box,n,1);
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
                const Vector<std::string>&       a_all_bc_type )
   {
      CH_TIME("FieldBCUtils::applyToJ_PIC()");

      // Note that this is called from pic species and from emfields when
      // computing J from mass matrices. Thus, a_all_bc_type could be that from
      // species or from fields, which are not necessarily the same. This leaves
      // room for a potential future bug. Need to fix this.
      //
      // Comment on J normal to boundary. For Symmetry, we subtract ghost Jnorm
      // from mirror location across boundary, which is the same as depositing
      // the mirror current to the interior. For other boundaries, we add ghost
      // Jnorm, which amount so using NGP for Jnorm close to boundary and is
      // needed to maintain exact energy conservation (only verified for insulator
      // BC so far)
      //

      // loop over non-periodic boundaries and apply BCs to J
      const BoundaryBoxLayoutPtrVect& bdry_layout = a_mesh.getBoundaryLayout();
      for (int b(0); b<bdry_layout.size(); b++) {

         const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
         const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
         const int bdry_dir = this_bdry_layout.dir();
         const int bdry_side(this_bdry_layout.side());
         const std::string this_bc = a_all_bc_type[b];

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
               if(bdry_side==0) edge_box.setSmall(bdry_dir,edge_box.bigEnd(bdry_dir));
               if(bdry_side==1) edge_box.setBig(bdry_dir,edge_box.smallEnd(bdry_dir));
                  
               FArrayBox& this_J( a_J_inPlane[interior_dit][dir] );
               if(dir==bdry_dir) {
                  Box dst_box = edge_box;
                  Box src_box = edge_box;
                  const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
                  for (int n=0; n<nG; n++) {
                     if(bdry_side==0) dst_box.shift(bdry_dir,1);
                     if(bdry_side==1) dst_box.shift(bdry_dir,-1);
                     if(this_bc=="symmetry" || this_bc=="axis") {
                        this_J.minus(this_J,src_box,dst_box,0,0,1);
                     }
                     else {
                        this_J.plus(this_J,src_box,dst_box,0,0,1);
                     }
                     if(bdry_side==0) src_box.shift(bdry_dir,-1);
                     if(bdry_side==1) src_box.shift(bdry_dir,1);
                  }

               }
               else {
                  if(this_bc=="symmetry" || this_bc=="axis") { 
                                            // add ghost J|| back to interior. Needed
                                            // for shape functions higher order than CIC and
                                            // for getting J|| correct for particles with xpbar
                                            // beyond a symmetry boundary during the intermediate 
                                            // stage of the advance, where there is no assumption 
                                            // that particle BCs have been applied.
                     Box dst_box = edge_box;
                     Box src_box = edge_box;
                     const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
                     for (int n=0; n<nG; n++) {
                        dst_box.shift(bdry_dir,1-2*bdry_side);
                        src_box.shift(bdry_dir,2*bdry_side-1);
                        this_J.plus(this_J,src_box,dst_box,0,0,1);
                     }
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
                                      // add ghost J|| back to interior. Needed
                                      // for shape functions higher order than CIC and
                                      // for getting J|| correct for particles with xpbar
                                      // beyond a symmetry boundary during the intermediate 
                                      // stage of the advance, where there is no assumption 
                                      // that particle BCs have been applied.

     	       if(this_bc=="axis") { // enforce Jth = 0 right on axis
		                     // 1 + 1 = 2 factor for Jz (done using corrected Jacobian)
				     // 1 - 1 = 0 factor for Jth (enforced here)
	          Box axis_box = node_box;
                  axis_box.setBig(0,0);
                  axis_box.setSmall(0,0);
	          this_J.setVal(0.0,axis_box,0,1);
	       }
     
	       Box dst_box = node_box;
               Box src_box = node_box;
               const int nG = bdry_box.bigEnd(bdry_dir)-bdry_box.smallEnd(bdry_dir)+1;
               for (int n=0; n<nG; n++) {
                  dst_box.shift(bdry_dir,1-2*bdry_side);
                  src_box.shift(bdry_dir,2*bdry_side-1);
                  if(this_bc=="axis") {
                     this_J.minus(this_J,src_box,dst_box,0,0,1); // Jth
#if CH_SPACEDIM==1
                     this_J.plus(this_J,src_box,dst_box,1,1,1); // Jz
#endif
                  }
                  else this_J.plus(this_J,src_box,dst_box,0,0,this_J.nComp());
               }
	    }
#endif
            
         } // end loop over boxes

      } // end loop over boundaries

   } // end applyToJ_PIC()
   
}
#include "NamespaceFooter.H"

#endif
