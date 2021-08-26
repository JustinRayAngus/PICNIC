#ifndef _INSULATOR_BC_UTILS_H_
#define _INSULATOR_BC_UTILS_H_

#include "Box.H"
#include "LoHiSide.H"
#include "DisjointBoxLayout.H"
#include "IntVect.H"
#include "BoxIterator.H"
#include "EdgeDataBox.H"
#include "SpaceUtils.H"
#include "LevelData.H"

#include "BoundaryBoxLayout.H"
#include "BoundaryConditions.H"

#include "NamespaceHeader.H"

namespace InsulatorBCUtils {

   inline
   void defineInsulatorBinary( LevelData<FArrayBox>&  a_ICbinary,
                         const BoundaryBoxLayout&     a_bdry_layout,
                         const RealVect&              a_Xmin_insulator,
                         const RealVect&              a_Xmax_insulator,
                         const DomainGrid&            a_mesh )
   {
      CH_TIME("InsulatorBCUtils::defineInsulatorBinary()");
 
      // Define a LevelData<FArrayBox> on bdry_grids that contains zeros 
      // at coords that touch conductors and ones at physical locations that 
      // touch insulator. (What about axis?)

      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<FArrayBox>& Xphys_cc = a_mesh.getXcc();

      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         FArrayBox& this_ICbinary( a_ICbinary[interior_dit] );

         // grow fill box to include tranverse ghosts 
         IntVect grow_vect = a_ICbinary.ghostVect();
         grow_vect[a_bdry_layout.dir()] = 0;
         Box fill_box_grown = grow(fill_box,grow_vect);
         
         BoxIterator bit(fill_box_grown); // grid indices for boundary box    
         IntVect ib;         // cell index
         RealVect local_Xcc; // cell center coords at cell index
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            Real this_val = 1.0;
            for(int dir=0; dir<SpaceDim; dir++) {
               local_Xcc[dir] = Xphys_cc[interior_dit].get(ib,dir);
               if(local_Xcc[dir]>a_Xmax_insulator[dir] || 
                  local_Xcc[dir]<a_Xmin_insulator[dir]) {
                  this_val = 0.0;
               } 
            } 
            this_ICbinary.set(ib,0,this_val); // conductor
         } 

      }
      
   }
   
   inline
   void defineInsulatorBinary( LevelData<FluxBox>&  a_ICbinary,
                         const BoundaryBoxLayout&   a_bdry_layout,
                         const RealVect&            a_Xmin_insulator,
                         const RealVect&            a_Xmax_insulator,
                         const DomainGrid&          a_mesh )
   {
      CH_TIME("InsulatorBCUtils::defineInsulatorBinary() Fluxes");
 
      // Define LevelData<FluxBox> on bdry_grids that contains zeros 
      // at coords that touch conductors and ones at physical locations that 
      // touch insulator. (What about axis?)

      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<FluxBox>& Xphys_fc = a_mesh.getXfc();

      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         
         for (int dir=0; dir<SpaceDim; dir++) {         
            FArrayBox& this_ICbinary( a_ICbinary[interior_dit][dir] );

            // grow fill box to include tranverse ghosts 
            IntVect grow_vect = a_ICbinary.ghostVect();
            grow_vect[a_bdry_layout.dir()] = 0;
            Box fill_box_grown = grow(fill_box,grow_vect);
            fill_box_grown.surroundingNodes(dir);       
  
            BoxIterator bit(fill_box_grown); // grid indices for boundary box    
            IntVect ib;         // cell index
            RealVect local_Xfc; // cell center coords at cell index
            for(bit.begin(); bit.ok(); ++bit) {
               ib = bit();
               Real this_val = 1.0;
               for(int dir0=0; dir0<SpaceDim; dir0++) {
                  local_Xfc[dir0] = Xphys_fc[interior_dit][dir].get(ib,dir0);
                  if(local_Xfc[dir0]>a_Xmax_insulator[dir0] || 
                     local_Xfc[dir0]<a_Xmin_insulator[dir0]) {
                     this_val = 0.0;
                  } 
               } 
               this_ICbinary.set(ib,0,this_val); // conductor
            } 
         }
      }     
   }
   
   inline
   void defineInsulatorBinary( LevelData<EdgeDataBox>&  a_ICbinary,
                         const BoundaryBoxLayout&       a_bdry_layout,
                         const RealVect&                a_Xmin_insulator,
                         const RealVect&                a_Xmax_insulator,
                         const DomainGrid&              a_mesh )
   {
      CH_TIME("InsulatorBCUtils::defineInsulatorBinary() Edges");
 
      // Define LevelData<EdgeDataBox> on bdry_grids that contains zeros 
      // at coords that touch conductors and ones at physical locations that 
      // touch insulator. (What about axis?)

      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<EdgeDataBox>& Xphys_ec = a_mesh.getXec();

      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         
         for (int dir=0; dir<SpaceDim; dir++) {         
            FArrayBox& this_ICbinary( a_ICbinary[interior_dit][dir] );

            // grow fill box to include tranverse ghosts 
            IntVect grow_vect = a_ICbinary.ghostVect();
            grow_vect[a_bdry_layout.dir()] = 0;
            Box fill_box_grown = grow(fill_box,grow_vect);
            fill_box_grown.surroundingNodes();
            fill_box_grown.enclosedCells(dir);        
 
            BoxIterator bit(fill_box_grown); // grid indices for boundary box    
            IntVect ib;         // cell index
            RealVect local_Xec; // cell center coords at cell index
            for(bit.begin(); bit.ok(); ++bit) {
               ib = bit();
               Real this_val = 1.0;
               for(int dir0=0; dir0<SpaceDim; dir0++) {
                  local_Xec[dir0] = Xphys_ec[interior_dit][dir].get(ib,dir0);
                  if(local_Xec[dir0]>a_Xmax_insulator[dir0] || 
                     local_Xec[dir0]<a_Xmin_insulator[dir0]) {
                     this_val = 0.0;
                  } 
               }
               this_ICbinary.set(ib,0,this_val); // conductor
            } 
         }
      }     
   }
   
   inline
   void defineInsulatorBinary( LevelData<NodeFArrayBox>&  a_ICbinary,
                         const BoundaryBoxLayout&         a_bdry_layout,
                         const RealVect&                  a_Xmin_insulator,
                         const RealVect&                  a_Xmax_insulator,
                         const DomainGrid&                a_mesh )
   {
      CH_TIME("InsulatorBCUtils::defineInsulatorBinary() Nodes");
 
      // Define a LevelData<FArrayBox> on bdry_grids that contains zeros 
      // at coords that touch conductors and ones at physical locations that 
      // touch insulator. (What about axis?)

      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<NodeFArrayBox>& Xphys_nc = a_mesh.getXnc();

      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         FArrayBox& this_ICbinary( a_ICbinary[interior_dit].getFab() );

         // grow fill box to include tranverse ghosts 
         IntVect grow_vect = a_ICbinary.ghostVect();
         grow_vect[a_bdry_layout.dir()] = 0;
         Box fill_box_grown = grow(fill_box,grow_vect);
         fill_box_grown.surroundingNodes();         

         BoxIterator bit(fill_box_grown); // grid indices for boundary box    
         IntVect ib;         // cell index
         RealVect local_Xnc; // cell center coords at cell index
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            Real this_val = 1.0;
            for(int dir=0; dir<SpaceDim; dir++) {
               local_Xnc[dir] = Xphys_nc[interior_dit].getFab().get(ib,dir);
               if(local_Xnc[dir]>a_Xmax_insulator[dir] || 
                  local_Xnc[dir]<a_Xmin_insulator[dir]) {
                  this_val = 0.0;
               } 
            } 
            this_ICbinary.set(ib,0,this_val); // conductor
         } 

      }
      
   }
   
}
#include "NamespaceFooter.H"

#endif
