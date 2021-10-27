
#include "DomainGrid.H"
#include "DomainGridF_F.H"
#include <array>
#include <cmath>
#include "BoxIterator.H"
#include "SpaceUtils.H"

#include "NamespaceHeader.H"

DomainGrid::DomainGrid( ParmParse&          a_ppgrid,
                  const int                 a_numGhosts,
                  const ProblemDomain&      a_domain,
                  const DisjointBoxLayout&  a_grids )
   : m_axisymmetric(false),
     m_ghosts(a_numGhosts),
     m_mapped_cell_volume(1.0)
    // m_mapped_face_area(SpaceDim,1.0)
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

   // compute mapped cell volume and face areas
   //
   for (int dir=0; dir<SpaceDim; ++dir) {
      m_mapped_cell_volume *= m_dX[dir];
      m_mapped_face_area[dir] = 1.0;
      for (int tdir=0; tdir<SpaceDim; ++tdir) {
         if (tdir != dir) m_mapped_face_area[dir] *= m_dX[tdir];
      }
   }


   if (a_ppgrid.contains("axisymmetric")) {
       a_ppgrid.get("axisymmetric", m_axisymmetric);

      if (m_axisymmetric && SpaceDim==3) {
         MayDay::Error("DomainGrid(): Cannot specify axisymmetry in 3D");
      }
   }
   
   int grid_verbosity;
   a_ppgrid.query( "verbosity", grid_verbosity );
   if(!procID() && grid_verbosity) {
      //cout << "====================== Spatial Grid Parameters =====================" << endl;
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
      if(SpaceDim<3) cout << "  axisymmetric = " << m_axisymmetric << endl;
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
  
   // define masked nodal Jacobians used for area/volume integrals  
 
   m_masked_Jfc.define(m_grids,1,IntVect::Zero);
   m_masked_Jec.define(m_grids,1,IntVect::Zero);
   m_masked_Jnc.define(m_grids,1,IntVect::Zero);

   // first define via copy
   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      for(int dir=0; dir<SpaceDim; dir++) {
         m_masked_Jfc[dit][dir].copy(m_Jfc[dit][dir]);
         m_masked_Jec[dit][dir].copy(m_Jec[dit][dir]);
      }
      m_masked_Jnc[dit].getFab().copy(m_Jnc[dit].getFab()); 
   }
   
   // mask the face centered Jacobian
   for(DataIterator dit(m_grids); dit.ok(); ++dit) {
      for(int dir=0; dir<SpaceDim; dir++) {
         Box face_box = m_grids[dit];
         face_box.surroundingNodes(dir);
         Box internal_box = face_box;
         internal_box.grow(dir,-1);
         
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
      for (int adir=0; adir<SpaceDim; adir++) internal_box.grow(adir,-1);

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

