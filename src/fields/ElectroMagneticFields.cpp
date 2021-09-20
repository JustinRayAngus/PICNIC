
#include <array>
#include <cmath>
#include "Constants.H"
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "GridFunctionFactory.H"

#include "ElectroMagneticFields.H"
#include "ElectroMagneticFieldsF_F.H"

#include "MathUtils.H"

#include "NamespaceHeader.H"

ElectroMagneticFields::ElectroMagneticFields( ParmParse&   a_ppflds,
                                        const DomainGrid&  a_mesh,
                                        const CodeUnits&   a_units,
                                        const bool&        a_verbosity )
    : m_mesh(a_mesh),
      m_verbosity(a_verbosity),
      m_use_poisson(false),
      m_field_bc(NULL),
      m_writeRho(false),
      m_writeDivs(false),
      m_writeCurls(false)
{
  
   // parse input file
   a_ppflds.get( "advance", m_advance );
   
   if(m_advance) {
      m_advanceE_inPlane = true;
      m_advanceE_virtual = true;
      m_advanceE_comp = {true,true,true};
      if(a_ppflds.contains("advance_electric_field")) {
         vector<int> this_advanceE_comp(3);
         a_ppflds.getarr( "advance_electric_field", this_advanceE_comp, 0, 3 );
         for (int dim=0; dim<3; dim++) m_advanceE_comp[dim] = (this_advanceE_comp[dim] == 1);
         m_advanceE_inPlane = false;
         for (int dim=0; dim<SpaceDim; dim++) m_advanceE_inPlane |= m_advanceE_comp[dim];
         m_advanceE_virtual = false;
         for (int dim=SpaceDim; dim<3; dim++) m_advanceE_virtual |= m_advanceE_comp[dim];
      }

      m_advanceB_inPlane = true;
      m_advanceB_virtual = true;
      m_advanceB_comp = {true,true,true};
      if(a_ppflds.contains("advance_magnetic_field")) {
         vector<int> this_advanceB_comp(3);
         a_ppflds.getarr( "advance_magnetic_field", this_advanceB_comp, 0, 3 );
         for (int dim=0; dim<3; dim++) m_advanceB_comp[dim] = (this_advanceB_comp[dim] == 1);
         m_advanceB_inPlane = false;
         for (int dim=0; dim<SpaceDim; dim++) m_advanceB_inPlane |= m_advanceB_comp[dim];
         m_advanceB_virtual = false;
         for (int dim=SpaceDim; dim<3; dim++) m_advanceB_virtual |= m_advanceB_comp[dim];
      }
   }
   else {
      m_advanceE_inPlane = false;
      m_advanceE_virtual = false;
      m_advanceE_comp = {false,false,false};
      m_advanceB_inPlane = false;
      m_advanceB_virtual = false;
      m_advanceB_comp = {false,false,false};
   }   
   m_advanceE = m_advanceE_inPlane;
   m_advanceE |= m_advanceE_virtual;
   m_advanceB = m_advanceB_inPlane;
   m_advanceB |= m_advanceB_virtual;

   a_ppflds.query( "write_rho", m_writeRho );
   a_ppflds.query( "write_divs", m_writeDivs );
   a_ppflds.query( "write_curls", m_writeCurls );

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
   
   const Real ep0  = Constants::EP0;  // SI units
   const Real Escale = a_units.getScale(a_units.ELECTRIC_FIELD);
   m_rhoCnorm_factor = qe/ep0*Xscale/Escale;
   
   if(!procID()){
      cout << " m_Jnorm_factor = " << m_Jnorm_factor << endl;
      cout << " advance fields = " << m_advance << endl;
      if(m_advance) {
         cout << " advance_E = " << m_advanceE_comp[0] << " " 
                                 << m_advanceE_comp[1] << " " 
                                 << m_advanceE_comp[2] << endl;
         cout << " advance_B = " << m_advanceB_comp[0] << " " 
                                 << m_advanceB_comp[1] << " " 
                                 << m_advanceB_comp[2] << endl;
      }
      cout << " write curls = " << m_writeCurls << endl;
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
   m_magneticField_diff.define(grids,1);            // FluxBox with 1 comp
   if(SpaceDim<3) {
      m_magneticField_virtual.define(grids,3-SpaceDim,ghostVect);      // FArrayBox
      m_magneticField_virtual_old.define(grids,3-SpaceDim,ghostVect);  // FArrayBox
      m_magneticField_virtual_diff.define(grids,3-SpaceDim); // FArrayBox
   }

   // initialize the electric field containers
   m_electricField.define(grids,1,ghostVect);         // EdgeDataBox with 1 comp
   m_electricField_old.define(grids,1,IntVect::Zero); // EdgeDataBox with 1 comp
   m_electricField_diff.define(grids,1,IntVect::Zero);// EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_electricField_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
      m_electricField_virtual_old.define(grids,3-SpaceDim,IntVect::Zero); // NodeFArrayBox
      m_electricField_virtual_diff.define(grids,3-SpaceDim,IntVect::Zero);// NodeFArrayBox
   }
   
   // initialize the current density containers
   m_currentDensity.define(grids,1,ghostVect);    // EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
   }
   m_chargeDensity.define(grids,1,IntVect::Zero);
   
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

   // set in-plane curls to zero for 1D
   if(SpaceDim==1) {
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_curlB[dit][0].setVal(0.0);
         m_curlE[dit][0].setVal(0.0);
      }
   }
   
   // initialize the divE and divB containers
   m_divE.define(grids,1,IntVect::Zero);
   m_divB.define(grids,1,IntVect::Zero);
   
   // define pointer to the BC object
   if(m_advance) m_field_bc = new FieldBC( m_mesh, m_verbosity );
   
   // initialize poisson solver
   a_ppflds.query( "use_poisson", m_use_poisson );
   if(m_use_poisson) {

      // parse the input file for the poisson solver
      getPoissonParameters(m_poisson_params, m_mesh);

      // setup the dbl for the elliptic solver
      int numlevels = m_poisson_params.numLevels;
      //Vector<DisjointBoxLayout> vectGrids(numlevels,grids);
      Vector<DisjointBoxLayout> vectGrids;
      setGrids(vectGrids,  m_poisson_params);

      for (int ilev = 0; ilev < vectGrids.size(); ilev++) { 
         //vectGrids[ilev] = grids;
         if(!procID()) {
            cout << "using poisson corrector" << endl;
            cout << "grid at level " << ilev << " = " << vectGrids[ilev] << endl;
         }
      }

      // define the data containers for the potential solver
      m_phi_vector.resize(numlevels);
      m_rhs_vector.resize(numlevels);
      int ncomp  = 1;
      for (int ilev = 0; ilev < numlevels; ilev++) {
         m_phi_vector[ilev] = new LevelData<NodeFArrayBox>(vectGrids[ilev], ncomp, ghostVect);
         m_rhs_vector[ilev] = new LevelData<NodeFArrayBox>(vectGrids[ilev], ncomp, ghostVect);
         LevelData<NodeFArrayBox>& phifabs= *m_phi_vector[ilev];
         LevelData<NodeFArrayBox>& rhsfabs= *m_rhs_vector[ilev];
      }
      m_electricField_correction.define(grids,1,ghostVect);  // EdgeDataBox with 1 comp
      m_phi0.define(grids,1,ghostVect);
      m_rhs0.define(grids,1,ghostVect);

      // define the elliptic solver
      m_bottomSolver.m_verbosity = 1;
      nodeDefineSolver(m_poisson_solver, vectGrids, m_bottomSolver, m_poisson_params); 
   
   }

}

ElectroMagneticFields::~ElectroMagneticFields()
{

   if(m_field_bc!=NULL) {
      delete m_field_bc;
      m_field_bc = NULL;
   }
   for (int n=0; n<m_phi_vector.size(); n++) {
      delete m_phi_vector[n];
      delete m_rhs_vector[n];
   }

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
   }
   
   // set initial electric field profile from ICs
   for (int dir=0; dir<SpaceDim; dir++) {
      stringstream sE; 
      sE << "IC.em_fields.electric." << dir; 
      ParmParse ppEfieldIC( sE.str().c_str() );
      RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppEfieldIC,m_verbosity);
      gridFunction->assign( m_electricField, dir, m_mesh, this_time );
   }
   
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
   }
   
   // apply BCs to initial profiles
   applyBCs_electricField(0.0);
   applyBCs_magneticField(0.0);
 
   // set the initial divs
   setDivE();
   setDivB();

   // set the initial curls
   setCurlE();
   setCurlB();

   // initialize the old values via copy
   updateOldElectricField();
   updateOldMagneticField();

   if(!procID()) cout << "Finished initializing electromagnetic fields " << endl << endl;

}

void ElectroMagneticFields::advanceElectricField_2ndHalf( const Real&  a_theta )
{
   CH_TIME("ElectroMagneticFields::advanceElectricField_2ndHalf()");
   
   // E^{n+1} = C1*E^{n+theta} + C0*E^n

   Real C1 = 1.0/a_theta; 
   Real C0 = (a_theta - 1.0)/a_theta; 

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   if(m_advanceE_inPlane) {
   
      for(DataIterator dit(grids); dit.ok(); ++dit) {

         for (int dir=0; dir<SpaceDim; dir++) {
            if(m_advanceE_comp[dir]) {
               FArrayBox& Edir = m_electricField[dit][dir];
               const FArrayBox& Edir_old = m_electricField_old[dit][dir];
               Edir.mult(C1);
               Edir.plus(Edir_old,C0,0,0,1);
            }
         }

      }

      SpaceUtils::exchangeEdgeDataBox(m_electricField);
   }
   
   if(SpaceDim<3 && m_advanceE_virtual) {
      
      for(DataIterator dit(grids); dit.ok(); ++dit) {
      
         FArrayBox& Ev = m_electricField_virtual[dit].getFab();
         const FArrayBox& Ev_old = m_electricField_virtual_old[dit].getFab();
         const int Ncomp = Ev.nComp();
         for (int comp=0; comp<Ncomp; comp++) {
            if(m_advanceE_comp[SpaceDim+comp]) {
               Ev.mult(C1,comp,1);
               Ev.plus(Ev_old,C0,comp,comp,1);
            }
         }

      }

      SpaceUtils::exchangeNodeFArrayBox(m_electricField_virtual,m_mesh); // experimental
   }
   
   // apply the BCs and perform exchange   
   //Real dummy_time = 0.0;
   //applyBCs_electricField(dummy_time);

}

void ElectroMagneticFields::advanceMagneticField_2ndHalf( const Real&  a_theta )
{
   CH_TIME("ElectroMagneticFields::advanceMagneticField_2ndHalf()");
   
   // B^{n+1} = C1*B^{n+theta} + C0*B^n
   
   Real C1 = 1.0/a_theta; 
   Real C0 = (a_theta - 1.0)/a_theta; 
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   if(m_advanceB_inPlane) {
 
      for(DataIterator dit(grids); dit.ok(); ++dit) {

         for (int dir=0; dir<SpaceDim; dir++) {
            if(m_advanceB_comp[dir]) {
               FArrayBox& Bdir = m_magneticField[dit][dir];
               const FArrayBox& Bdir_old = m_magneticField_old[dit][dir];
               Bdir.mult(C1);
               Bdir.plus(Bdir_old,C0,0,0,1);
            }
         }

      }

   }  
   
   if(SpaceDim<3 && m_advanceB_virtual) {
      
      for(DataIterator dit(grids); dit.ok(); ++dit) {

         FArrayBox& Bv = m_magneticField_virtual[dit];
         const FArrayBox& Bv_old = m_magneticField_virtual_old[dit];
         const int Ncomp = Bv.nComp();
         for (int comp=0; comp<Ncomp; comp++) {
            if(m_advanceB_comp[SpaceDim+comp]) {
               Bv.mult(C1,comp,1);
               Bv.plus(Bv_old,C0,comp,comp,1);
            }
         }

      }

   }

   // apply the BCs and perform exchange   
   //Real dummy_time = 0.0;
   //applyBCs_magneticField(dummy_time);

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

   if(m_advanceB_inPlane) {

      for(DataIterator dit(grids); dit.ok(); ++dit) {
      
         for (int dir=0; dir<SpaceDim; dir++) {
            if(m_advanceB_comp[dir]) {
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
   
      }

   }
 
   if(SpaceDim<3 && m_advanceB_virtual) {
         
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         
         Box cell_box = grids[dit];
         
         const int Ncomp = m_magneticField_virtual.nComp();
         for (int comp=0; comp<Ncomp; comp++){
            if(m_advanceB_comp[SpaceDim+comp]) {
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

      }

   }
   
   // apply the BCs and perform exchange   
   //Real dummy_time = 0.0;
   //applyBCs_magneticField(dummy_time);

}

void ElectroMagneticFields::advanceElectricField( const Real& a_cnormDt )
{
   CH_TIME("ElectroMagneticFields::advanceElectricField()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());

   if(m_advanceE_inPlane) {

      for(DataIterator dit(grids); dit.ok(); ++dit) {
      
         for (int dir=0; dir<SpaceDim; dir++) {
            if(m_advanceE_comp[dir]) {
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

      }

   }

   if(SpaceDim<3 && m_advanceE_virtual) {
   
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         
         Box node_box = grids[dit];
         node_box.surroundingNodes();
         
         const int Ncomp = m_electricField_virtual.nComp();
         for (int comp=0; comp<Ncomp; comp++){
            if(m_advanceE_comp[SpaceDim+comp]) {
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

      }

   }
   
   // apply the BCs and perform exchange   
   //Real dummy_time = 0.0;
   //applyBCs_electricField(dummy_time);

}

void ElectroMagneticFields::applyBCs_electricField( const Real  a_time )
{
   CH_TIME("ElectroMagneticFields::applyBCs_electricField()");
   
   if(m_advanceE_inPlane) {
     
      // exchange to fill internal ghost cells 
      SpaceUtils::exchangeEdgeDataBox(m_electricField);
     
      // apply the boundary bcs 
      m_field_bc->applyEdgeBC( m_electricField, a_time );

   }
   
   if(SpaceDim<3 && m_advanceE_virtual) {
      
      // exchange to fill internal ghost cells 
      SpaceUtils::exchangeNodeFArrayBox(m_electricField_virtual,m_mesh); // experimental
      
      // apply the boundary bcs 
      m_field_bc->applyNodeBC( m_electricField_virtual, a_time );
   
   }

}

void ElectroMagneticFields::applyBCs_magneticField( const Real  a_time )
{
   CH_TIME("ElectroMagneticFields::applyBCs_magneticField()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());

   if(m_advanceB_inPlane) {
      
      // exchange to fill internal ghost cells 
      SpaceUtils::exchangeFluxBox(m_magneticField);
      
      // apply the boundary bcs 
      m_field_bc->applyFluxBC( m_magneticField, a_time );

   }
 
   if(SpaceDim<3 && m_advanceB_virtual) {
         
      // exchange to fill internal ghost cells 
      m_magneticField_virtual.exchange();
      
      // apply the boundary bcs 
      m_field_bc->applyCellBC( m_magneticField_virtual, a_time );

   }

}

void ElectroMagneticFields::saveElectricField()
{
   CH_TIME("ElectroMagneticFields::saveElectricField()");

   SpaceUtils::copyLevelData( m_electricField, 
                              m_electricField_diff );
   if(SpaceDim<3) {
      SpaceUtils::copyLevelData(  m_electricField_virtual, 
                                  m_electricField_virtual_diff );
   }
}

void ElectroMagneticFields::saveMagneticField()
{
   CH_TIME("ElectroMagneticFields::saveMagneticField()");

   SpaceUtils::copyLevelData( m_magneticField, 
                              m_magneticField_diff );
   if(SpaceDim<3) {
      SpaceUtils::copyLevelData(  m_magneticField_virtual, 
                                  m_magneticField_virtual_diff );
   }
}
   
Real ElectroMagneticFields::diffElectricField()
{
   CH_TIME("ElectroMagneticFields::diffElectricField()");

   SpaceUtils::addToLevelData(  m_electricField, 
                                m_electricField_diff, 
                                -1  );
   if(SpaceDim<3) {
      SpaceUtils::addToLevelData( m_electricField_virtual, 
                                  m_electricField_virtual_diff, 
                                  -1  );
   }

   return computeEFieldDiffNorm();
}

Real ElectroMagneticFields::diffMagneticField()
{
   CH_TIME("ElectroMagneticFields::diffMagneticField()");

   SpaceUtils::addToLevelData(  m_magneticField, 
                                m_magneticField_diff, 
                                -1  );
   if(SpaceDim<3) {
      SpaceUtils::addToLevelData( m_magneticField_virtual, 
                                  m_magneticField_virtual_diff, 
                                  -1  );
   }

   return computeBFieldDiffNorm();
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

void ElectroMagneticFields::setDivE()
{
   CH_TIME("ElectroMagneticFields::setDivE()");
    
   // get some Grid info
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const RealVect& dX(m_mesh.getdX());
  
   for(DataIterator dit(grids); dit.ok(); ++dit) {
         
      FArrayBox& this_divE( m_divE[dit].getFab() );
      Box node_box = grids[dit];
      node_box.surroundingNodes();
         
      this_divE.setVal(0.0, node_box, 0, this_divE.nComp());
      for (int dir=0; dir<SpaceDim; dir++) {
         int additive = 1;
         const FArrayBox& this_E( m_electricField[dit][dir] );
         SpaceUtils::simpleStagGradComp( this_divE, node_box, 0,
                                         this_E, 0,
                                         dX[dir], dir,
                                         additive ); 
      }
   }
   
}

void ElectroMagneticFields::setDivB()
{
   CH_TIME("ElectroMagneticFields::setDivB()");
    
   // get some Grid info
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const RealVect& dX(m_mesh.getdX());
   
   for(DataIterator dit(grids); dit.ok(); ++dit) {
         
      FArrayBox& this_divB( m_divB[dit] );
      const Box& cell_box = grids[dit];
         
      this_divB.setVal(0.0, cell_box, 0, this_divB.nComp());
      for (int dir=0; dir<SpaceDim; dir++) {
         int additive = 1;
         const FArrayBox& this_B( m_magneticField[dit][dir] );
         SpaceUtils::simpleStagGradComp( this_divB, cell_box, 0,
                                         this_B, 0,
                                         dX[dir], dir,
                                         additive ); 
      }
   }
   
}

void ElectroMagneticFields::solvePoisson( const Real&  a_time )
{
   CH_TIME("ElectroMagneticFields::solvePoisson()");

   //
   // 1) solve for correction potential: nabla^2(phi) = div(E) - rhoC/ep0
   // 2) compute the correction field:   E_corr = -nabla(phi)
   //

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit; 
  
   /* 
   GridFunctionFactory  gridFactory;
   LevelData<NodeFArrayBox> temporaryNodeProfile;
   temporaryNodeProfile.define(grids,1,m_electricField_virtual.ghostVect());
   
   string prhs("poisson.rhs"); 
   ParmParse pp( prhs );
   RefCountedPtr<GridFunction> gridFunction = gridFactory.create(pp,m_verbosity);
   gridFunction->assign( temporaryNodeProfile, m_mesh, a_time );
   */

   // set the rhs for container with default dbl
   int numlevels = m_poisson_params.numLevels;
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      Box dst_box = grids[dit];
      dst_box.surroundingNodes();
      FArrayBox& this_rhs = m_rhs0[dit].getFab();
      const FArrayBox& this_divE = m_divE[dit].getFab();
      const FArrayBox& this_rhoC = m_chargeDensity[dit].getFab();
      this_rhs.copy(this_rhoC,dst_box,0,dst_box,0,1);
      this_rhs.mult(-m_rhoCnorm_factor);
      this_rhs.plus(this_divE,dst_box,0,0,1);
      //this_rhs.copy(temporaryNodeProfile[dit].getFab(),dst_box,0,dst_box,0,1);
   }
   SpaceUtils::exchangeNodeFArrayBox(m_rhs0,m_mesh);

   // use copyTo to put rhs on container with dbl used for poisson solver
   LevelData<NodeFArrayBox>& rhs_vect0= *m_rhs_vector[0];
   m_rhs0.copyTo(rhs_vect0);

   // solve for phi
   int lbase = 0; 
   m_poisson_solver.solve(m_phi_vector, m_rhs_vector, numlevels-1, lbase);
   
   // use copyTo to fill phi container with default dbl
   LevelData<NodeFArrayBox>& phi_vect0= *m_phi_vector[0];
   phi_vect0.copyTo(m_phi0);
      
   // compute E_corr = -nabla(phi)
   const RealVect& dX(m_mesh.getdX());
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& phi = m_phi0[dit].getFab();
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& deltaE = m_electricField_correction[dit][dir];
         const Box& edge_box = deltaE.box();
            
         int additive = 0;
         SpaceUtils::simpleStagGradComp( deltaE, edge_box, 0,
                                         phi, 0,
                                         -dX[dir], dir,
                                         additive ); 
      }
   }
   SpaceUtils::exchangeEdgeDataBox(m_electricField_correction);

}

void ElectroMagneticFields::correctElectricField()
{
   CH_TIME("ElectroMagneticFields::correctElectricField()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
  
   // add E_corr = -nabla(phi) to E
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& Edir = m_electricField[dit][dir];
         const FArrayBox& deltaE = m_electricField_correction[dit][dir];
         const Box& edge_box = deltaE.box();
         Edir.plus(deltaE,edge_box,0,0,1);
      }
   }

}   
   
#include "NamespaceFooter.H"

