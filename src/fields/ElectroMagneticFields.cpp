
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
                                        const CodeUnits&   a_units,
                                        const bool&        a_verbosity )
   : m_verbosity(a_verbosity),
     m_mesh(a_mesh)
{
   // parse input file
   a_ppflds.get( "advance", m_advance );

   // set the Courant time step associated with the speed of light
   const RealVect& meshSpacing(m_mesh.getdX());
   Real invDXsq = 0.0;
   for (int dir=0; dir<SpaceDim; dir++) {
      invDXsq = invDXsq + 1.0/meshSpacing[dir]/meshSpacing[dir];
   } 
   Real DXeff = 1.0/pow(invDXsq,0.5);
   m_stable_dt = DXeff/a_units.CvacNorm();

   // set the normalization factor for J in Ampere's law
   const Real mu0  = Constants::MU0;  // SI units
   const Real qe   = Constants::QE;   // SI units
   const Real cvac = Constants::CVAC; // SI units
   const Real Xscale = a_units.getScale(a_units.LENGTH);
   const Real Bscale = a_units.getScale(a_units.MAGNETIC_FIELD);
   m_Jnorm_factor = mu0*qe*cvac*Xscale/Bscale;
   
   if(!procID()){
      cout << " m_Jnorm_factor = " << m_Jnorm_factor << endl;
      cout << " advance fields = " << m_advance << endl;
      cout << " stable dt = " << m_stable_dt << endl << endl;
   }
   
   //
   // initialize the member LevelDatas
   //
  
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   
   // initialize the magnetic field containers
   m_magneticField.define(grids,1,ghostVect);       // FluxBox with 1 comp
   m_magneticField_old.define(grids,1,ghostVect);   // FluxBox with 1 comp
   if(SpaceDim<3) {
      m_magneticField_virtual.define(grids,3-SpaceDim,ghostVect);      // FArrayBox
      m_magneticField_virtual_old.define(grids,3-SpaceDim,ghostVect);  // FArrayBox
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
   updateOldElectricField();
   updateOldMagneticField();

   if(!procID()) cout << "Finished initializing electromagnetic fields " << endl << endl;

}

void ElectroMagneticFields::advanceElectricField_2ndHalf()
{
   CH_TIME("ElectroMagneticFields::advanceElectricField_2ndHalf()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& Edir = m_electricField[dit][dir];
         const FArrayBox& Edir_old = m_electricField_old[dit][dir];
         const Box& edge_box = Edir.box();
         Edir.mult(2.0);
         Edir.minus(Edir_old,edge_box,0,0,1);
      }
   
      if(SpaceDim<3) {
         FArrayBox& Ev = m_electricField_virtual[dit].getFab();
         const FArrayBox& Ev_old = m_electricField_virtual_old[dit].getFab();
         const Box& node_box = Ev.box();
         const int Ncomp = Ev.nComp();
         Ev.mult(2.0);
         Ev.minus(Ev_old,node_box,0,0,Ncomp);
      }

   }
   SpaceUtils::exchangeEdgeDataBox(m_electricField);
   if(SpaceDim<3) SpaceUtils::exchangeNodeFArrayBox(m_electricField_virtual,m_mesh); // experimental

}

void ElectroMagneticFields::advanceMagneticField_2ndHalf()
{
   CH_TIME("ElectroMagneticFields::advanceMagneticField_2ndHalf()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& Bdir = m_magneticField[dit][dir];
         const FArrayBox& Bdir_old = m_magneticField_old[dit][dir];
         const Box& face_box = Bdir.box();
         Bdir.mult(2.0);
         Bdir.minus(Bdir_old,face_box,0,0,1);
      }
   
      if(SpaceDim<3) {
         FArrayBox& Bv = m_magneticField_virtual[dit];
         const FArrayBox& Bv_old = m_magneticField_virtual_old[dit];
         const Box& cell_box = Bv.box();
         const int Ncomp = Bv.nComp();
         Bv.mult(2.0);
         Bv.minus(Bv_old,cell_box,0,0,Ncomp);
      }

   }
   SpaceUtils::exchangeFluxBox(m_magneticField);
   if(SpaceDim<3) m_magneticField_virtual.exchange();

}

void ElectroMagneticFields::updateOldElectricField()
{
   CH_TIME("ElectroMagneticFields::updateOldElectricField()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& Edir_old = m_electricField_old[dit][dir];
         const Box& edge_box = Edir_old.box();
         Edir_old.copy(m_electricField[dit][dir],edge_box);
      }
      
      FArrayBox& Ev_old = m_electricField_virtual_old[dit].getFab();
      const Box& node_box = Ev_old.box();
      Ev_old.copy(m_electricField_virtual[dit].getFab(),node_box);

   }

}

void ElectroMagneticFields::updateOldMagneticField()
{
   CH_TIME("ElectroMagneticFields::updateOldMagneticField()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());

   for(DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& Bdir_old = m_magneticField_old[dit][dir];
         const Box& face_box = Bdir_old.box();
         Bdir_old.copy(m_magneticField[dit][dir],face_box);
      }   
      
      if(SpaceDim<3) {
         FArrayBox& Bv_old = m_magneticField_virtual_old[dit];
         const Box& cell_box = Bv_old.box();
         Bv_old.copy(m_magneticField_virtual[dit],cell_box);
      }

   }

}

void ElectroMagneticFields::advanceMagneticField( const Real& a_cnormDt )
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
                                           CHF_CONST_REAL(a_cnormDt),
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
                                              CHF_CONST_REAL(a_cnormDt),
                                              CHF_CONST_FRA1(Bv_old,comp), 
                                              CHF_CONST_FRA1(curlEv,comp), 
                                              CHF_FRA1(Bvdir,comp) );
         }

      }
      m_magneticField_virtual.exchange();
   
   }
   
}

void ElectroMagneticFields::advanceElectricField( const Real& a_cnormDt )
{
   CH_TIME("ElectroMagneticFields::advanceElectricField()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir=0; dir<SpaceDim; dir++) {

         const FArrayBox& Edir_old = m_electricField_old[dit][dir];
         const FArrayBox& curlBdir = m_curlB[dit][dir];
         const FArrayBox& Jdir     = m_currentDensity[dit][dir];
               FArrayBox& Edir     = m_electricField[dit][dir];

         Box edge_box = grids[dit];
         edge_box.surroundingNodes();
         edge_box.enclosedCells(dir);

         FORT_ADVANCE_ELECTRIC_FIELD_COMP( CHF_BOX(edge_box),
                                           CHF_CONST_REAL(a_cnormDt),
                                           CHF_CONST_FRA1(Edir_old,0), 
                                           CHF_CONST_FRA1(curlBdir,0), 
                                           CHF_CONST_FRA1(Jdir,0), 
                                           CHF_CONST_REAL(m_Jnorm_factor),
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
            const FArrayBox& Jv     = m_currentDensity_virtual[dit].getFab();
                  FArrayBox& Evdir  = m_electricField_virtual[dit].getFab();

            FORT_ADVANCE_ELECTRIC_FIELD_COMP( CHF_BOX(node_box),
                                              CHF_CONST_REAL(a_cnormDt),
                                              CHF_CONST_FRA1(Ev_old,comp), 
                                              CHF_CONST_FRA1(curlBv,comp), 
                                              CHF_CONST_FRA1(Jv,comp), 
                                              CHF_CONST_REAL(m_Jnorm_factor),
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

   SpaceUtils::exchangeEdgeDataBox(m_curlB);   
   if(SpaceDim<3) SpaceUtils::exchangeNodeFArrayBox(m_curlB_virtual,m_mesh); // experimental
   
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
   
   SpaceUtils::exchangeFluxBox(m_curlE);
   if(SpaceDim<3) m_curlE_virtual.exchange();   
   
}
  
void ElectroMagneticFields::setCurrentDensity( const PicSpeciesPtrVect&  a_pic_species_ptr_vect )
{
   CH_TIME("ElectroMagneticFields::setCurrentDensity()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
    
   // set the current density member to zero
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         m_currentDensity[dit][dir].setVal(0.0);
      }
      if(SpaceDim<3) {
         FArrayBox& this_Jv_nodes = m_currentDensity_virtual[dit].getFab();
         this_Jv_nodes.setVal(0.0);
      }
   }

   // loop over all pic species and add the current density 
   for (int s=0; s<a_pic_species_ptr_vect.size(); ++s) {
      const PicSpeciesPtr this_picSpecies(a_pic_species_ptr_vect[s]);
      int this_charge = this_picSpecies->charge();
      if(this_charge != 0) {

         //this_picSpecies->setCurrentDensity(); // set J for this_picSpecies prior to call here
         const LevelData<EdgeDataBox>& species_J = this_picSpecies->getCurrentDensity();
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir=0; dir<SpaceDim; dir++) {
               const FArrayBox& this_species_Jdir = species_J[dit][dir];
               m_currentDensity[dit][dir].plus(this_species_Jdir,0,0,1);
            }
         }
         if(SpaceDim<3) {
            const LevelData<NodeFArrayBox>& species_Jv = this_picSpecies->getCurrentDensity_virtual();
            for(DataIterator dit(grids); dit.ok(); ++dit) {
               const FArrayBox& this_species_Jv = species_Jv[dit].getFab();
               FArrayBox& this_Jv = m_currentDensity_virtual[dit].getFab();
               this_Jv.plus(this_species_Jv,0,0,m_currentDensity_virtual.nComp());
            }
         }

      }
   }

   // call exchange
   //SpaceUtils::exchangeEdgeDataBox(m_currentDensity);
   //if(SpaceDim<3) SpaceUtils::exchangeNodeFArrayBox(m_currentDensity_virtual,m_mesh);

}


#include "NamespaceFooter.H"

