
#include <array>
#include <cmath>
#include "Constants.H"
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "GridFunctionFactory.H"

#include "ElectroMagneticFields.H"
#include "ElectroMagneticFieldsF_F.H"

#include "MathUtils.H"
#include "SpaceUtils.H"

#include "NamespaceHeader.H"

ElectroMagneticFields::ElectroMagneticFields( ParmParse&   a_ppflds,
                                        const DomainGrid&  a_mesh,
                                        const bool&        a_verbosity )
   : m_verbosity(a_verbosity),
     m_mesh(a_mesh)
{
   // parse input file
   a_ppflds.get( "advance", m_advance );

   if ( procID() == 0 ) {
      cout << " advance fields = " << m_advance << endl << endl;
   }
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
   const int ghosts(m_mesh.ghosts());

   m_stable_dt = meshSpacing[0]/Constants::CVAC;
   for (int dir=1; dir<SpaceDim; dir++) {
      m_stable_dt = Min(m_stable_dt,meshSpacing[dir]/Constants::CVAC);
   } 

   //
   // initialize the member LevelDatas
   //
   //
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   
   // initialize the magnetic field containers
   m_magneticField.define(grids,1,ghostVect);     // FluxBox with 1 comp
   m_magneticField_old.define(grids,1,IntVect::Zero);     // FluxBox with 1 comp
   if(SpaceDim<3) {
      m_magneticField_virtual.define(grids,3-SpaceDim,ghostVect); // FArrayBox
      m_magneticField_virtual_old.define(grids,3-SpaceDim,IntVect::Zero); // FArrayBox
   }

   // initialize the electric field containers
   m_electricField.define(grids,1,ghostVect);     // EdgeDataBox with 1 comp
   m_electricField_old.define(grids,1,IntVect::Zero);     // EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_electricField_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
      m_electricField_virtual_old.define(grids,3-SpaceDim,IntVect::Zero); // NodeFArrayBox
   }
   
   // initialize the current density containers
   m_currentDensity.define(grids,1,ghostVect);    // EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
   }
   
   // initialize the curlE containers
   m_curlE.define(grids,1,IntVect::Zero);     // FluxBox with 1 comp
   if(SpaceDim<3) {
      m_curlE_virtual.define(grids,3-SpaceDim,IntVect::Zero); // FArrayBox
   }
   
   // initialize the curlB containers
   m_curlB.define(grids,1,IntVect::Zero);     // EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_curlB_virtual.define(grids,3-SpaceDim,IntVect::Zero); // NodeFArrayBox
   }
   if(SpaceDim==1) {
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_curlB[dit][0].setVal(0.0);
      }
   }

}

ElectroMagneticFields::~ElectroMagneticFields()
{
}

void ElectroMagneticFields::initialize()
{
   if(!procID()) cout << "Initializing electromagnetic fields..." << endl;
   
   //
   // parse the initial profiles for the fields
   //
   
   GridFunctionFactory  gridFactory;
   const Real this_time = 0.0;
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   LevelData<FArrayBox> temporaryProfile;
   temporaryProfile.define(grids,1,m_magneticField_virtual.ghostVect());
   LevelData<NodeFArrayBox> temporaryNodeProfile;
   temporaryNodeProfile.define(grids,1,m_electricField_virtual.ghostVect());
   
   // set initial magnetic field profile from ICs
   for (int dir=0; dir<SpaceDim; dir++) {
      stringstream sB; 
      sB << "IC.em_fields.magnetic." << dir; 
      ParmParse ppBfieldIC( sB.str().c_str() );
      RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppBfieldIC,m_verbosity);
      gridFunction->assign( m_magneticField, dir, m_mesh, this_time );
   }
   SpaceUtils::exchangeFluxBox(m_magneticField);
   
   // set initial virtual magnetic field profile from ICs
   if(SpaceDim<3) {
      int numVirt = 3-SpaceDim;
      for( int comp=0; comp<numVirt; comp++) {
         int field_dir = comp + SpaceDim;
         stringstream sBv; 
         sBv << "IC.em_fields.magnetic." << field_dir; 
         ParmParse ppBvfieldIC( sBv.str().c_str() );
         RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppBvfieldIC,m_verbosity);
         gridFunction->assign( temporaryProfile, m_mesh, this_time );
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            m_magneticField_virtual[dit].copy(temporaryProfile[dit],0,comp,1);
         }
      }
      m_magneticField_virtual.exchange();
   }
   
   // set initial electric field profile from ICs
   for (int dir=0; dir<SpaceDim; dir++) {
      stringstream sE; 
      sE << "IC.em_fields.electric." << dir; 
      ParmParse ppEfieldIC( sE.str().c_str() );
      RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppEfieldIC,m_verbosity);
      gridFunction->assign( m_electricField, dir, m_mesh, this_time );
   }
   SpaceUtils::exchangeEdgeDataBox(m_electricField);
   
   // set initial virtual electric field profile from ICs
   if(SpaceDim<3) {
      int numVirt = 3-SpaceDim;
      for( int comp=0; comp<numVirt; comp++) {
         int field_dir = comp + SpaceDim;
         stringstream sEv; 
         sEv << "IC.em_fields.electric." << field_dir; 
         ParmParse ppEvfieldIC( sEv.str().c_str() );
         RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppEvfieldIC,m_verbosity);
         gridFunction->assign( temporaryNodeProfile, m_mesh, this_time );
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox& this_Ev  = m_electricField_virtual[dit].getFab();
            this_Ev.copy(temporaryNodeProfile[dit].getFab(),0,comp,1);
            //Box dst_box = grids[dit];
            //dst_box.surroundingNodes();
            //this_Ev.copy(temporaryNodeProfile[dit].getFab(),dst_box,0,dst_box,comp,1);
         }
      }
      SpaceUtils::exchangeNodeFArrayBox(m_electricField_virtual,m_mesh); // experimental
   }
 
   // initialize the old values via copy
   updateOldFieldValues();

   if(!procID()) cout << "Finished initializing electromagnetic fields " << endl << endl;

}

void ElectroMagneticFields::updateOldFieldValues()
{
   CH_TIME("ElectroMagneticFields::updateOldFieldValues()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {

         FArrayBox& Bdir_old = m_magneticField_old[dit][dir];
         const Box& face_box = Bdir_old.box();
         Bdir_old.copy(m_magneticField[dit][dir],face_box);
         
         FArrayBox& Edir_old = m_electricField_old[dit][dir];
         const Box& edge_box = Edir_old.box();
         Edir_old.copy(m_electricField[dit][dir],edge_box);
         
      }
   }

   if(SpaceDim<3) {
      for(DataIterator dit(grids); dit.ok(); ++dit) {

         FArrayBox& Bv_old = m_magneticField_virtual_old[dit];
         const Box& cell_box = Bv_old.box();
         Bv_old.copy(m_magneticField_virtual[dit],cell_box);

         FArrayBox& Ev_old = m_electricField_virtual_old[dit].getFab();
         Box node_box = Ev_old.box();
         node_box.surroundingNodes();
         Ev_old.copy(m_electricField_virtual[dit].getFab(),node_box);

      }
   }

}

void ElectroMagneticFields::advanceMagneticField( const Real& a_dt )
{
   CH_TIME("ElectroMagneticFields::advanceMagneticField()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
    
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir=0; dir<SpaceDim; dir++) {

         const FArrayBox& Bdir_old = m_magneticField_old[dit][dir];
         const FArrayBox& curlEdir = m_curlE[dit][dir];
               FArrayBox& Bdir     = m_magneticField[dit][dir];

         Box face_box = grids[dit];
         face_box.surroundingNodes(dir);

         FORT_ADVANCE_MAGNETIC_FIELD_COMP( CHF_BOX(face_box),
                                           CHF_CONST_REAL(a_dt),
                                           CHF_CONST_FRA1(Bdir_old,0), 
                                           CHF_CONST_FRA1(curlEdir,0), 
                                           CHF_FRA1(Bdir,0) );
        
      }
   
   }
   SpaceUtils::exchangeFluxBox(m_magneticField);
 
   if(SpaceDim<3) {
         
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         
         Box cell_box = grids[dit];
         
         const int Ncomp = m_magneticField_virtual.nComp();
         for (int comp=0; comp<Ncomp; comp++){

            const FArrayBox& Bv_old = m_magneticField_virtual_old[dit];
            const FArrayBox& curlEv = m_curlE_virtual[dit];
                  FArrayBox& Bvdir  = m_magneticField_virtual[dit];

            FORT_ADVANCE_MAGNETIC_FIELD_COMP( CHF_BOX(cell_box),
                                              CHF_CONST_REAL(a_dt),
                                              CHF_CONST_FRA1(Bv_old,comp), 
                                              CHF_CONST_FRA1(curlEv,comp), 
                                              CHF_FRA1(Bvdir,comp) );
         }

      }
      m_magneticField_virtual.exchange();
   
   }
   
}

void ElectroMagneticFields::advanceElectricField( const Real& a_dt )
{
   CH_TIME("ElectroMagneticFields::advanceElectricField()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   const Real dtcvac2 = Constants::CVAC*Constants::CVAC*a_dt;
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir=0; dir<SpaceDim; dir++) {

         const FArrayBox& Edir_old = m_electricField_old[dit][dir];
         const FArrayBox& curlBdir = m_curlB[dit][dir];
               FArrayBox& Edir     = m_electricField[dit][dir];

         Box edge_box = grids[dit];
         edge_box.surroundingNodes();
         edge_box.enclosedCells(dir);

         FORT_ADVANCE_ELECTRIC_FIELD_COMP( CHF_BOX(edge_box),
                                           CHF_CONST_REAL(dtcvac2),
                                           CHF_CONST_FRA1(Edir_old,0), 
                                           CHF_CONST_FRA1(curlBdir,0), 
                                           CHF_FRA1(Edir,0) );
        
      }
   }
   SpaceUtils::exchangeEdgeDataBox(m_electricField);

   if(SpaceDim<3) {
   
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         
         Box node_box = grids[dit];
         node_box.surroundingNodes();
         
         const int Ncomp = m_electricField_virtual.nComp();
         for (int comp=0; comp<Ncomp; comp++){

            const FArrayBox& Ev_old = m_electricField_virtual_old[dit].getFab();
            const FArrayBox& curlBv = m_curlB_virtual[dit].getFab();
                  FArrayBox& Evdir  = m_electricField_virtual[dit].getFab();

            FORT_ADVANCE_ELECTRIC_FIELD_COMP( CHF_BOX(node_box),
                                              CHF_CONST_REAL(dtcvac2),
                                              CHF_CONST_FRA1(Ev_old,comp), 
                                              CHF_CONST_FRA1(curlBv,comp), 
                                              CHF_FRA1(Evdir,comp) );
         }

      }
      SpaceUtils::exchangeNodeFArrayBox(m_electricField_virtual,m_mesh); // experimental
   
   }

}

void ElectroMagneticFields::setCurlB()
{
   CH_TIME("ElectroMagneticFields::setCurlB()");
    
   // get some Grid info
   //
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const RealVect& dX(m_mesh.getdX());
  
   if(SpaceDim==3) { // FluxBox ==> EdgeDataBox

      for(DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
           
            int dirj, dirk;
            dirj = dir + 1;
            dirj = dirj % SpaceDim;
            dirk = dir + 2;
            dirk = dirk % SpaceDim;

            FArrayBox& curlB_dir = m_curlB[dit][dir];
            const FArrayBox& B_dirj = m_magneticField[dit][dirj];           
            const FArrayBox& B_dirk = m_magneticField[dit][dirk];

            int additive = 0;
            SpaceUtils::simpleStagGradComp( curlB_dir, grids[dit], 0,
                                            B_dirk, 0,
                                            dX[dirj], dirj,
                                            additive ); 
            additive = 1;
            SpaceUtils::simpleStagGradComp( curlB_dir, grids[dit], 0,
                                            B_dirj, 0,
                                            -dX[dirk], dirk,
                                            additive ); 
            
         }
      }

   }
   
   if(SpaceDim==2) { // FluxBox ==> NodeFArrayBox and FArrayBox ==> EdgeDataBox

      for(DataIterator dit(grids); dit.ok(); ++dit) {
   
         //        
         // FluxBox ==> NodeFArrayBox
         //
 
         Box node_box = grids[dit];

         int dirj, dirk;
         dirj = 0;
         dirk = 1;

         FArrayBox& curlBv = m_curlB_virtual[dit].getFab();
         const FArrayBox& B_dirj = m_magneticField[dit][dirj];           
         const FArrayBox& B_dirk = m_magneticField[dit][dirk];

         int additive = 0;
         SpaceUtils::simpleStagGradComp( curlBv, grids[dit], 0,
                                         B_dirk, 0,
                                         dX[dirj], dirj,
                                         additive ); 
         additive = 1;
         SpaceUtils::simpleStagGradComp( curlBv, grids[dit], 0,
                                         B_dirj, 0,
                                         -dX[dirk], dirk,
                                         additive ); 
            
         //
         // FArrayBox ==> EdgeDataBox
         //
         
         for (int dir=0; dir<SpaceDim; dir++) {
           
            int dirj;
            dirj = dir + 1;
            dirj = dirj % SpaceDim;

            FArrayBox& curlB_dir = m_curlB[dit][dir];
            const FArrayBox& B_dir2 = m_magneticField_virtual[dit];           

            int additive = 0;
            int sign = 1-2*dir;
            SpaceUtils::simpleStagGradComp( curlB_dir, grids[dit], 0,
                                            B_dir2, 0,
                                            sign*dX[dirj], dirj,
                                            additive ); 
            
         }

      }

   }
   
   if(SpaceDim==1) { // 2 comp FArrayBox ==> 2 comp NodeFArrayBox

      for(DataIterator dit(grids); dit.ok(); ++dit) {
   
         int dirj = 0;
         int dirk = 1;

         FArrayBox& curlB_dir = m_curlB_virtual[dit].getFab();
         const FArrayBox& B_virt = m_magneticField_virtual[dit];           

         int additive = 0;
         SpaceUtils::simpleStagGradComp( curlB_dir, grids[dit], 0,
                                         B_virt, dirk,
                                         -dX[0], 0,
                                         additive ); 

         SpaceUtils::simpleStagGradComp( curlB_dir, grids[dit], 1,
                                         B_virt, dirj,
                                         dX[0], 0,
                                         additive ); 

      }

   }  
   
}

void ElectroMagneticFields::setCurlE()
{
   CH_TIME("ElectroMagneticFields::setCurlE()");
    
   // get some Grid info
   //
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const RealVect& dX(m_mesh.getdX());
  
   if(SpaceDim==3) { // EdgeDataBox ==> FluxBox

      for(DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
           
            int dirj, dirk;
            dirj = dir + 1;
            dirj = dirj % SpaceDim;
            dirk = dir + 2;
            dirk = dirk % SpaceDim;

            FArrayBox& curlE_dir = m_curlE[dit][dir];
            const FArrayBox& E_dirj = m_electricField[dit][dirj];           
            const FArrayBox& E_dirk = m_electricField[dit][dirk];

            int additive = 0;
            SpaceUtils::simpleStagGradComp( curlE_dir, grids[dit], 0,
                                            E_dirk, 0,
                                            dX[dirj], dirj,
                                            additive ); 
            additive = 1;
            SpaceUtils::simpleStagGradComp( curlE_dir, grids[dit], 0,
                                            E_dirj, 0,
                                            -dX[dirk], dirk,
                                            additive ); 
            
         }
      }

   }
   
   if(SpaceDim==2) { // EdgeDataBox ==> FArrayBox and NodeFArrayBox ==> FluxBox

      for(DataIterator dit(grids); dit.ok(); ++dit) {
   
         //        
         // EdgeDataBox ==> FArrayBox
         //

         int dirj, dirk;
         dirj = 0;
         dirk = 1;

         FArrayBox& curlEv = m_curlE_virtual[dit];
         const FArrayBox& E_dirj = m_electricField[dit][dirj];           
         const FArrayBox& E_dirk = m_electricField[dit][dirk];

         int additive = 0;
         SpaceUtils::simpleStagGradComp( curlEv, grids[dit], 0,
                                         E_dirk, 0,
                                         dX[dirj], dirj,
                                         additive ); 
         additive = 1;
         SpaceUtils::simpleStagGradComp( curlEv, grids[dit], 0,
                                         E_dirj, 0,
                                         -dX[dirk], dirk,
                                         additive ); 
         //
         // NodeFArrayBox ==> FluxBox
         //
         
         for (int dir=0; dir<SpaceDim; dir++) {
           
            int dirj;
            dirj = dir + 1;
            dirj = dirj % SpaceDim;

            FArrayBox& curlE_dir = m_curlE[dit][dir];
            const FArrayBox& E_dir2 = m_electricField_virtual[dit].getFab();           

            int additive = 0;
            int sign = 1-2*dir;
            SpaceUtils::simpleStagGradComp( curlE_dir, grids[dit], 0,
                                            E_dir2, 0,
                                            sign*dX[dirj], dirj,
                                            additive ); 
            
         }

      }

   }
   
   if(SpaceDim==1) { // 2 comp NodeFArrayBox ==> 2 comp FArrayBox

      for(DataIterator dit(grids); dit.ok(); ++dit) {
   
         int dirj = 0;
         int dirk = 1;

         FArrayBox& curlE_dir = m_curlE_virtual[dit];
         const FArrayBox& E_virt = m_electricField_virtual[dit].getFab();           

         int additive = 0;
         SpaceUtils::simpleStagGradComp( curlE_dir, grids[dit], 0,
                                         E_virt, dirk,
                                         -dX[0], 0,
                                         additive ); 

         SpaceUtils::simpleStagGradComp( curlE_dir, grids[dit], 1,
                                         E_virt, dirj,
                                         dX[0], 0,
                                         additive ); 

      }

   }  
   
}


#include "NamespaceFooter.H"

