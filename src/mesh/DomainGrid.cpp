
#include "DomainGrid.H"
#include "DomainGridF_F.H"
#include <array>
#include <cmath>
#include "BoxIterator.H"
#include "SpaceUtils.H"
#include "Constants.H"

#include "NamespaceHeader.H"

DomainGrid::DomainGrid( ParmParse&          a_ppgrid,
                  const int                 a_numGhosts,
                  const ProblemDomain&      a_domain,
                  const DisjointBoxLayout&  a_grids )
   : m_axisymmetric(false),
     m_anticyclic(false),
     m_write_jacobians(false),
     m_write_corrected_jacobians(false),
     m_ghosts(a_numGhosts),
     m_mapped_cell_volume(1.0)
{
   m_grids  = a_grids;
   m_domain = a_domain;
   IntVect dimensions = a_domain.size(); 
  
   // standard ghost exchange copier
   m_forwardCopier.define( m_grids, m_grids, 
                           m_domain, m_ghosts*IntVect::Unit, true );

   // a reversed version of the above
   m_reverseCopier.define( m_grids, m_grids, 
                           m_domain, m_ghosts*IntVect::Unit, true );
   m_reverseCopier.reverse();
   
   a_ppgrid.query( "write_jacobians", m_write_jacobians );
   a_ppgrid.query( "write_corrected_jacobians", m_write_corrected_jacobians );
   
   a_ppgrid.get( "geometry", m_geom_type );
   if(m_geom_type=="cartesian") {

      a_ppgrid.get("X_min", m_Xmin[0]);
      a_ppgrid.get("X_max", m_Xmax[0]);
      m_dX[0] = (m_Xmax[0] - m_Xmin[0])/(double)dimensions[0];    

      if(SpaceDim==2) {
         a_ppgrid.get("Z_min", m_Xmin[1]);
         a_ppgrid.get("Z_max", m_Xmax[1]);
         m_dX[1] = (m_Xmax[1] - m_Xmin[1])/(double)dimensions[1];    
      }
      if(SpaceDim==3) {
         a_ppgrid.get("Y_min", m_Xmin[1]);
         a_ppgrid.get("Y_max", m_Xmax[1]);
         m_dX[1] = (m_Xmax[1] - m_Xmin[1])/(double)dimensions[1];    
         a_ppgrid.get("Z_min", m_Xmin[2]);
         a_ppgrid.get("Z_max", m_Xmax[2]);
         m_dX[2] = (m_Xmax[2] - m_Xmin[2])/(double)dimensions[2];    
      }
   
   }
   else if(m_geom_type=="cyl_R") {
 
      CH_assert(SpaceDim==1);
      m_axisymmetric = true;

      a_ppgrid.get("R_min", m_Xmin[0]);
      a_ppgrid.get("R_max", m_Xmax[0]);
      m_dX[0] = (m_Xmax[0] - m_Xmin[0])/(double)dimensions[0];    

   }
   else if(m_geom_type=="cyl_RTH") {
 
      CH_assert(SpaceDim==2);

      a_ppgrid.get("R_min", m_Xmin[0]);
      a_ppgrid.get("R_max", m_Xmax[0]);
      m_dX[0] = (m_Xmax[0] - m_Xmin[0])/(double)dimensions[0];    
      
      a_ppgrid.get("TH_min", m_Xmin[1]);
      a_ppgrid.get("TH_max", m_Xmax[1]);
      m_dX[1] = (m_Xmax[1] - m_Xmin[1])/(double)dimensions[1];

   }
   else if(m_geom_type=="cyl_RZ") {
 
      CH_assert(SpaceDim==2);
      m_axisymmetric = true;
      m_anticyclic = true;

      if(a_ppgrid.contains("X_min") && a_ppgrid.contains("X_max")) {
         a_ppgrid.get("X_min", m_Xmin[0]);
         a_ppgrid.get("X_max", m_Xmax[0]);
      }
      else {
         a_ppgrid.get("R_min", m_Xmin[0]);
         a_ppgrid.get("R_max", m_Xmax[0]);
      }
      CH_assert(m_Xmin[0]>=0.0 && m_Xmax[0]>m_Xmin[0]);
      m_dX[0] = (m_Xmax[0] - m_Xmin[0])/(double)dimensions[0];    
      
      a_ppgrid.get("Z_min", m_Xmin[1]);
      a_ppgrid.get("Z_max", m_Xmax[1]);
      CH_assert(m_Xmax[1]>m_Xmin[1]);
      m_dX[1] = (m_Xmax[1] - m_Xmin[1])/(double)dimensions[1];

   }
   else {
      if(!procID()) cout << "m_geom_type " << m_geom_type << " not supported " << endl;
      exit(EXIT_FAILURE);
   }
   
   // compute mapped cell volume and face areas
   for (int dir=0; dir<SpaceDim; ++dir) {
      m_mapped_cell_volume *= m_dX[dir];
      m_mapped_face_area[dir] = 1.0;
      for(int tdir=0; tdir<SpaceDim; ++tdir) {
         if (tdir != dir) m_mapped_face_area[dir] *= m_dX[tdir];
      }
   }

   if(a_ppgrid.contains("axisymmetric")) { // can manually turn this on/off
      a_ppgrid.get("axisymmetric", m_axisymmetric);
      if(m_axisymmetric && SpaceDim==3) {
         MayDay::Error("DomainGrid(): Cannot specify axisymmetry in 3D");
      }
   }
   
   int grid_verbosity;
   a_ppgrid.query( "verbosity", grid_verbosity );
   if(!procID() && grid_verbosity) {
      //cout << "====================== Spatial Grid Parameters =====================" << endl;
      cout << " geometry = " << m_geom_type << endl;
      cout << " axisymmetric = " << m_axisymmetric << endl;
      cout << " anticyclic   = " << m_anticyclic << endl;
      cout << "  X_min, X_max = " << m_Xmin[0] << ", " << m_Xmax[0] << endl;
      cout << "  dX = " << m_dX[0] << endl;
      if(SpaceDim==2) {
          cout << "  Z_min, Z_max = " << m_Xmin[1] << ", " << m_Xmax[1] << endl;
          cout << "  dZ = " << m_dX[1] << endl;
      }
      if(SpaceDim==3) {
         cout << "  Y_min, Y_max = " << m_Xmin[1] << ", " << m_Xmax[1] << endl;
         cout << "  dY = " << m_dX[1] << endl;
         cout << "  Z_min, Z_max = " << m_Xmin[2] << ", " << m_Xmax[2] << endl;
         cout << "  dZ = " << m_dX[2] << endl;
      }
      cout << "====================================================================" << endl;
      cout << endl;
   }

   // set the physical coordinates
   setRealCoords();
   
   // set the Jacobian
   setJacobian();

   // define the vector of boundary box layouts for BCs
   defineBoundaryBoxLayout();
}

void DomainGrid::setRealCoords()
{
   IntVect ghostVect = m_ghosts*IntVect::Unit;
   m_Xcc.define(m_grids,SpaceDim,ghostVect);
   m_Xfc.define(m_grids,SpaceDim,ghostVect);
   m_Xec.define(m_grids,SpaceDim,ghostVect);
   m_Xnc.define(m_grids,SpaceDim,ghostVect);

   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      
      // set the coords at cell center
      FORT_GET_CC_MAPPED_COORDS( CHF_BOX(m_Xcc[dit].box()),
                                 CHF_CONST_REALVECT(m_dX),
                                 CHF_FRA(m_Xcc[dit]) );

      for (int dir=0; dir<SpaceDim; ++dir) {
         m_Xcc[dit].plus(m_Xmin[dir],dir,1);
      }
   
      // set the coords at cell faces
      for (int dir=0; dir<SpaceDim; ++dir) {
         FORT_GET_FC_MAPPED_COORDS( CHF_BOX(m_Xfc[dit][dir].box()),
                                    CHF_CONST_INT(dir),
                                    CHF_CONST_REALVECT(m_dX),
                                    CHF_FRA(m_Xfc[dit][dir]) );
      
         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            m_Xfc[dit][dir].plus(m_Xmin[tdir],tdir,1);
         }
      }
       
      // set the coords at cell edges
      for (int dir=0; dir<SpaceDim; ++dir) {
         FORT_GET_EC_MAPPED_COORDS( CHF_BOX(m_Xec[dit][dir].box()),
                                    CHF_CONST_INT(dir),
                                    CHF_CONST_REALVECT(m_dX),
                                    CHF_FRA(m_Xec[dit][dir]) );
      
         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            m_Xec[dit][dir].plus(m_Xmin[tdir],tdir,1);
         }
      }
      //if(!procID()) cout << "JRA: m_Xcc.box() = " << m_Xcc[dit].box() << endl;      
      //if(!procID()) cout << "JRA: m_Xnc.box() = " << m_Xnc[dit].box() << endl;      
      //if(!procID()) cout << "JRA: m_Xnc.getFab().box() = " << m_Xnc[dit].getFab().box() << endl;      
      // set the coords at cell nodes
      FORT_GET_NC_MAPPED_COORDS( CHF_BOX(surroundingNodes(m_Xnc[dit].box())),
                                 CHF_CONST_REALVECT(m_dX),
                                 CHF_FRA(m_Xnc[dit]) );
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_Xnc[dit].plus(m_Xmin[dir],dir,1);
      }

   }
   
}

void DomainGrid::setJacobian()
{
   IntVect ghostVect = m_ghosts*IntVect::Unit;
   m_Jcc.define(m_grids,1,ghostVect);
   m_Jfc.define(m_grids,1,ghostVect);
   m_Jec.define(m_grids,1,ghostVect);
   m_Jnc.define(m_grids,1,ghostVect);
   
   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      m_Jcc[dit].setVal(1.0); 
      for(int dir=0; dir<SpaceDim; dir++) {
         m_Jfc[dit][dir].setVal(1.0);
         m_Jec[dit][dir].setVal(1.0);
      }
      m_Jnc[dit].getFab().setVal(1.0); 
   }

   // modify Jacobians corresponding to geometry
   if(m_geom_type=="cyl_RTH") {

      for(DataIterator dit(m_grids); dit.ok(); ++dit) {
         m_Jcc[dit].mult(m_Xcc[dit],0,0,1); 
         for(int dir=0; dir<SpaceDim; dir++) {
            m_Jfc[dit][dir].mult(m_Xfc[dit][dir],0,0,1);
            m_Jec[dit][dir].mult(m_Xec[dit][dir],0,0,1);
         }
         m_Jnc[dit].getFab().mult(m_Xnc[dit].getFab(),0,0,1); 
      }

   }
   
   if( (m_geom_type=="cyl_R" || m_geom_type=="cyl_RZ") &&
        m_axisymmetric ) {

      const Real twoPi = Constants::TWOPI; 
      for(DataIterator dit(m_grids); dit.ok(); ++dit) {
         m_Jcc[dit].mult(m_Xcc[dit],0,0,1); 
         m_Jcc[dit].mult(twoPi,0,1);
         for(int dir=0; dir<SpaceDim; dir++) {
            m_Jfc[dit][dir].mult(m_Xfc[dit][dir],0,0,1);
            m_Jfc[dit][dir].mult(twoPi,0,1);
            m_Jec[dit][dir].mult(m_Xec[dit][dir],0,0,1);
            m_Jec[dit][dir].mult(twoPi,0,1);
         }
         m_Jnc[dit].getFab().mult(m_Xnc[dit].getFab(),0,0,1); 
         m_Jnc[dit].getFab().mult(twoPi,0,1);
      }

   }
   
   //
   // define corrected nodal Jacobians used for charge/current
   // density deposit (See J.P. Verboncoeur JCP 2001)
   //

   m_corrected_Jfc.define(m_grids,1,IntVect::Zero);
   m_corrected_Jec.define(m_grids,1,IntVect::Zero);
   m_corrected_Jnc.define(m_grids,1,IntVect::Zero);

   // first define via copy
   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      for(int dir=0; dir<SpaceDim; dir++) {
         m_corrected_Jfc[dit][dir].copy(m_Jfc[dit][dir]);
         m_corrected_Jec[dit][dir].copy(m_Jec[dit][dir]);
      }
      m_corrected_Jnc[dit].getFab().copy(m_Jnc[dit].getFab()); 
   }

   // redefine J on r-nodes if using cyl coords
   const Real Pi = Constants::PI; 
   auto phys_domain( m_grids.physDomain() );
   
   if( (m_geom_type=="cyl_R" || m_geom_type=="cyl_RZ") &&
        m_axisymmetric ) {
      
      for(DataIterator dit(m_grids); dit.ok(); ++dit) {

         // get the domain box on nodes
         Box domain_node_box = phys_domain.domainBox();
         domain_node_box.surroundingNodes();
         int dir0 = 0;
         int dir0_bdry_hi = domain_node_box.bigEnd(dir0);
         int dir0_bdry_lo = domain_node_box.smallEnd(dir0);

         // get Jacobian on nodes and define local node box
         FArrayBox& this_Jnc(m_corrected_Jnc[dit].getFab());
         Box node_box = m_grids[dit];
         node_box.surroundingNodes();

         Real local_J;
         IntVect ig;
         IntVect shift_vect = IntVect::Zero;
         shift_vect[dir0] = 1;
         
         BoxIterator gbit(node_box);
         for(gbit.begin(); gbit.ok(); ++gbit) {
            ig = gbit(); // grid index
            Real rj = m_Xnc[dit].getFab().get(ig,0);
            Real rjp1 = m_Xnc[dit].getFab().get(ig+shift_vect,0);
            Real rjm1 = m_Xnc[dit].getFab().get(ig-shift_vect,0);
            if(ig[dir0]==dir0_bdry_lo) {
                local_J = Pi/3.0*(rjp1 - rj)*(2*rj + rjp1)/m_dX[dir0];
            }
            else if(ig[dir0]==dir0_bdry_hi) {
                local_J = Pi/3.0*(rj - rjm1)*(2*rjm1 + 2*rj)/m_dX[dir0];
            }
            else {
                local_J = Pi/3.0*(rjp1*(rj + rjp1) - rjm1*(rjm1 + rj))/m_dX[dir0];
            }
            this_Jnc.set(ig,0,local_J);
         }

      }

   }
   
   // redefine J on r-edges if using cyl coords
   if( m_geom_type=="cyl_RZ" && m_axisymmetric ) {
      
      for(DataIterator dit(m_grids); dit.ok(); ++dit) {
         const int dir0=0;
         const int dir1=1;

         // get the domain box on edges
         Box domain_edge_box = phys_domain.domainBox();
         domain_edge_box.surroundingNodes();
         domain_edge_box.enclosedCells(dir1);
         int dir0_bdry_hi = domain_edge_box.bigEnd(dir0);
         int dir0_bdry_lo = domain_edge_box.smallEnd(dir0);

         // get Jacobian on edges and define local edge box
         FArrayBox& this_Jec(m_corrected_Jec[dit][dir1]);
         Box edge_box = m_grids[dit];
         edge_box.surroundingNodes();
         edge_box.enclosedCells(dir1);

         Real local_J;
         IntVect ig;
         IntVect shift_vect = IntVect::Zero;
         shift_vect[dir0] = 1;
         BoxIterator gbit(edge_box);
         for(gbit.begin(); gbit.ok(); ++gbit) {
            ig = gbit(); // grid index
            Real rj = m_Xec[dit][dir1].get(ig,0);
            Real rjp1 = m_Xec[dit][dir1].get(ig+shift_vect,0);
            Real rjm1 = m_Xec[dit][dir1].get(ig-shift_vect,0);
            if(ig[dir0]==dir0_bdry_lo) {
                local_J = Pi/3.0*(rjp1 - rj)*(2*rj + rjp1)/m_dX[dir0];
            }
            else if(ig[dir0]==dir0_bdry_hi) {
                local_J = Pi/3.0*(rj - rjm1)*(2*rjm1 + 2*rj)/m_dX[dir0];
            }
            else {
                local_J = Pi/3.0*(rjp1*(rj + rjp1) - rjm1*(rjm1 + rj))/m_dX[dir0];
            }
            this_Jec.set(ig,0,local_J);
         }
      }

   }
   
   // define corected J on r-face if using cyl coords
   if( (m_geom_type=="cyl_R" || m_geom_type=="cyl_RZ")
     && m_axisymmetric ) {
      
      for(DataIterator dit(m_grids); dit.ok(); ++dit) {
         const int dir0=0;

         // get the domain box on face
         Box domain_face_box = phys_domain.domainBox();
         domain_face_box.surroundingNodes(dir0);
         int dir0_bdry_hi = domain_face_box.bigEnd(dir0);
         int dir0_bdry_lo = domain_face_box.smallEnd(dir0);

         // get Jacobian on faces and define local face box
         FArrayBox& this_Jfc(m_corrected_Jfc[dit][dir0]);
         Box face_box = m_grids[dit];
         face_box.surroundingNodes(dir0);

         Real local_J;
         IntVect ig;
         IntVect shift_vect = IntVect::Zero;
         shift_vect[dir0] = 1;
         BoxIterator gbit(face_box);
         for(gbit.begin(); gbit.ok(); ++gbit) {
            ig = gbit(); // grid index
            Real rj = m_Xfc[dit][dir0].get(ig,0);
            Real rjp1 = m_Xfc[dit][dir0].get(ig+shift_vect,0);
            Real rjm1 = m_Xfc[dit][dir0].get(ig-shift_vect,0);
            if(ig[dir0]==dir0_bdry_lo) {
                local_J = Pi/3.0*(rjp1 - rj)*(2*rj + rjp1)/m_dX[dir0];
            }
            else if(ig[dir0]==dir0_bdry_hi) {
                local_J = Pi/3.0*(rj - rjm1)*(2*rjm1 + 2*rj)/m_dX[dir0];
            }
            else {
                local_J = Pi/3.0*(rjp1*(rj + rjp1) - rjm1*(rjm1 + rj))/m_dX[dir0];
            }
            this_Jfc.set(ig,0,local_J);
         }
      }

   }

   //
   // define masked nodal Jacobians used for area/volume integrals
   // Will eventually need to define these using corrected Jacobians
   // Note that exact energy conservation requires dV used in current
   // density to match dV used in field energy integral. Need to
   // think about BCs/geometry when deriving masked from corrected Ja. 
   // 

   m_masked_Jfc.define(m_grids,1,IntVect::Zero);
   m_masked_Jec.define(m_grids,1,IntVect::Zero);
   m_masked_Jnc.define(m_grids,1,IntVect::Zero);

   // first define via copy
   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      for(int dir=0; dir<SpaceDim; dir++) {
         m_masked_Jfc[dit][dir].copy(m_corrected_Jfc[dit][dir]);
         m_masked_Jec[dit][dir].copy(m_corrected_Jec[dit][dir]);
      }
      m_masked_Jnc[dit].getFab().copy(m_corrected_Jnc[dit].getFab()); 
   }
   
   // mask the face centered Jacobian
   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      for(int dir=0; dir<SpaceDim; dir++) {
         Box face_box = m_grids[dit];
         face_box.surroundingNodes(dir);
         Box internal_box = face_box;
         internal_box.grow(dir,-1);

         // correct internal_box on physical boundaries for axisymm
         if( (m_geom_type=="cyl_R" || m_geom_type=="cyl_RZ") && m_axisymmetric ) {
             Box domain_face_box = phys_domain.domainBox();
             domain_face_box.surroundingNodes(dir);
             if(face_box.bigEnd(dir)==domain_face_box.bigEnd(dir)) {
                internal_box.growHi(dir,1);
             }
             if(face_box.smallEnd(dir)==domain_face_box.smallEnd(dir)) {
                internal_box.growLo(dir,1);
             }
         }
         
         FArrayBox mask(face_box,1);
         mask.setVal(0.5);
         mask.plus(0.5,internal_box,0,1); 

         m_masked_Jfc[dit][dir].mult(mask);
      }
   }
   
   // mask the edge centered Jacobian
   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      for(int dir=0; dir<SpaceDim; dir++) {
         Box edge_box = m_grids[dit];
         edge_box.surroundingNodes();
         edge_box.enclosedCells(dir);
         Box internal_box = edge_box;
         for (int adir=0; adir<SpaceDim; adir++) {
            if(adir==dir) continue;
            internal_box.grow(adir,-1);
         }
         
         // correct internal_box on physical boundaries for axisymm
         if( (m_geom_type=="cyl_R" || m_geom_type=="cyl_RZ") && m_axisymmetric ) {
             Box domain_edge_box = phys_domain.domainBox();
             domain_edge_box.surroundingNodes();
             domain_edge_box.enclosedCells(dir);
             for (int adir=0; adir<SpaceDim; adir++) {
                if (adir==dir) continue;
                if(edge_box.bigEnd(adir)==domain_edge_box.bigEnd(adir)) {
                   internal_box.growHi(adir,1);
                }
                if(edge_box.smallEnd(adir)==domain_edge_box.smallEnd(adir)) {
                   internal_box.growLo(adir,1);
                }
             }
         }

         FArrayBox mask(edge_box,1);
         if(SpaceDim==1) mask.setVal(1.0);
         if(SpaceDim==2) {
            mask.setVal(0.5);
            mask.plus(0.5,internal_box,0,1); 
         }
         if(SpaceDim==3) {
            mask.setVal(0.25);
            mask.plus(0.75,internal_box,0,1); 
         }

         m_masked_Jec[dit][dir].mult(mask);
      }
   }
   
   // mask the node centered Jacobian
   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      Box node_box = m_grids[dit];
      node_box.surroundingNodes();
      Box internal_box = node_box;
      for (int dir=0; dir<SpaceDim; dir++) internal_box.grow(dir,-1);
         
      // correct internal_box on physical boundaries for axisymm
      // Is the correct for SpaceDim > 1 ?...
      if( (m_geom_type=="cyl_R" || m_geom_type=="cyl_RZ") && m_axisymmetric ) {
         Box domain_node_box = phys_domain.domainBox();
         domain_node_box.surroundingNodes();
         for (int dir=0; dir<SpaceDim; dir++) {
            if(node_box.bigEnd(dir)==domain_node_box.bigEnd(dir)) {
               internal_box.growHi(dir,1);
            }
            if(node_box.smallEnd(dir)==domain_node_box.smallEnd(dir)) {
               internal_box.growLo(dir,1);
            }
         }
      }

      FArrayBox mask(node_box,1);
      if(SpaceDim==1) {
         mask.setVal(0.5);
         mask.plus(0.5,internal_box,0,1);
      }
      if(SpaceDim==2) { // corners get 1/4, edges get 1/2
         mask.setVal(0.25);
         for(int dir=0; dir<SpaceDim; dir++) {
            Box node_box0 = node_box;
            node_box0.grow(dir,-1);
            mask.plus(0.25,node_box0,0,1);
         }
         mask.plus(0.25,internal_box,0,1);
      }
      if(SpaceDim==3) { // corners get 1/8, edges get 1/4, faces get 1/2
         mask.setVal(0.125);
         for(int dir=0; dir<SpaceDim; dir++) {
            Box node_box0 = node_box;
            node_box0.grow(dir,-1);
            mask.plus(0.125,node_box0,0,1);
            //
            Box internal_box0 = internal_box;
            internal_box0.grow(dir,1);
            mask.plus(0.125,internal_box0,0,1);
         }
         mask.plus(0.125,internal_box,0,1);
      }

      m_masked_Jnc[dit].getFab().mult(mask);
   }

}

void DomainGrid::defineBoundaryBoxLayout()
{

   const IntVect ghostVect = m_ghosts*IntVect::Unit;
   const Box domain_box = m_domain.domainBox();

   for(int dir=0; dir<SpaceDim; dir++) {
      if(m_domain.isPeriodic(dir)) {
         for(SideIterator si; si.ok(); ++si) {
            Side::LoHiSide side( si() );
            m_periodic_bdry_layout.push_back(
            BoundaryBoxLayoutPtr( new BoundaryBoxLayout( m_grids,
                                                         domain_box,
                                                         dir,
                                                         side,
                                                         ghostVect )));
         }
      }
      else {
         for(SideIterator si; si.ok(); ++si) {
            Side::LoHiSide side( si() );
            m_domain_bdry_layout.push_back(
            BoundaryBoxLayoutPtr( new BoundaryBoxLayout( m_grids,
                                                         domain_box,
                                                         dir,
                                                         side,
                                                         ghostVect )));
         }
      }
   }

}


#include "NamespaceFooter.H"

