#include <array>
#include <cmath>
#include "Constants.H"
#include "BoxIterator.H"
#include "GridFunctionFactory.H"
#include "CH_HDF5.H"
#include "MathUtils.H"
#include "ElectroMagneticFields.H"

#include "NamespaceHeader.H"

ElectroMagneticFields::ElectroMagneticFields( ParmParse&    a_ppflds,
                                        const DomainGrid&   a_mesh,
                                        const CodeUnits&    a_units,
                                        const bool&         a_verbosity,
                                        const EMVecType&    a_vec_type )
    : m_verbosity(a_verbosity),
      m_use_poisson(false),
      m_writeDivs(false),
      m_writeCurls(false),
      m_writeRho(false),
      m_mesh(a_mesh),
      m_field_bc(NULL),
      m_vec_size_bfield(0),
      m_vec_size_bfield_virtual(0),
      m_vec_size_efield(0),
      m_vec_size_efield_virtual(0),
      m_vec_offset_bfield(-1),
      m_vec_offset_bfield_virtual(-1),
      m_vec_offset_efield(-1),
      m_vec_offset_efield_virtual(-1)
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
   
   const Real volume_scale = a_units.getScale(a_units.VOLUME);
   const Real dVolume = m_mesh.getMappedCellVolume()*volume_scale; // SI units
   m_energyE_factor = 0.5*Escale*Escale*ep0*dVolume;
   m_energyB_factor = 0.5*Bscale*Bscale/mu0*dVolume;
   
   if(!procID()) {
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
   // define the member LevelDatas
   //
  
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   
   // define the magnetic field containers
   m_magneticField.define(grids,1,ghostVect);       // FluxBox with 1 comp
   m_magneticField_rhs.define(grids,1);            // FluxBox with 1 comp
   if(SpaceDim<3) {
      m_magneticField_virtual.define(grids,3-SpaceDim,ghostVect);      // FArrayBox
      m_magneticField_virtual_rhs.define(grids,3-SpaceDim); // FArrayBox
   }

   // define the electric field containers
   m_electricField.define(grids,1,ghostVect);         // EdgeDataBox with 1 comp
   m_electricField_rhs.define(grids,1,IntVect::Zero);// EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_electricField_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
      m_electricField_virtual_rhs.define(grids,3-SpaceDim,IntVect::Zero);// NodeFArrayBox
   }
   
   // define the current density containers
   m_currentDensity.define(grids,1,ghostVect);    // EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
   }
   m_chargeDensity.define(grids,1,IntVect::Zero);
   
   // define the curlE containers
   m_curlE.define(grids,1,IntVect::Zero);     // FluxBox with 1 comp
   if(SpaceDim<3) {
      m_curlE_virtual.define(grids,3-SpaceDim,IntVect::Zero); // FArrayBox
   }
   
   // define the curlB containers
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
   
   // define the divE and divB containers
   m_divE.define(grids,1,IntVect::Zero);
   m_divB.define(grids,1,IntVect::Zero);
   
   // define pointer to the BC object
   if(m_advance) m_field_bc = new FieldBC( m_mesh, m_verbosity );
   
   // define poisson solver
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
      }
      m_electricField_correction.define(grids,1,ghostVect);  // EdgeDataBox with 1 comp
      m_phi0.define(grids,1,ghostVect);
      m_rhs0.define(grids,1,ghostVect);

      // define the elliptic solver
      m_bottomSolver.m_verbosity = 1;
      nodeDefineSolver(m_poisson_solver, vectGrids, m_bottomSolver, m_poisson_params); 
   
   }

   m_vec_size_bfield = SpaceUtils::nDOF( m_magneticField );
   m_vec_size_efield = SpaceUtils::nDOF( m_electricField );
   if (SpaceDim < 3) {
     m_vec_size_bfield_virtual = SpaceUtils::nDOF( m_magneticField_virtual );
     m_vec_size_efield_virtual = SpaceUtils::nDOF( m_electricField_virtual );
   }

   int vec_size_total = 0;

   if ((a_vec_type == e_and_b) || (a_vec_type == b_only)) {

     m_vec_offset_bfield = vec_size_total;
     vec_size_total += m_vec_size_bfield;

     if (SpaceDim < 3) {
       m_vec_offset_bfield_virtual = vec_size_total;
       vec_size_total += m_vec_size_bfield_virtual;
     }

   }

   if ((a_vec_type == e_and_b) || (a_vec_type == e_only)) {
     m_vec_offset_efield = vec_size_total;
     vec_size_total += m_vec_size_efield;

     if (SpaceDim < 3) {
       m_vec_offset_efield_virtual = vec_size_total;
       vec_size_total += m_vec_size_efield_virtual;
     }

   }

   return;
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

void ElectroMagneticFields::initialize( const Real          a_cur_time,
                                        const std::string&  a_restart_file_name )
{
   const DisjointBoxLayout& grids(m_mesh.getDBL());
 
   if(a_restart_file_name.empty()) {
   
      if(!procID()) cout << "Initializing electromagnetic fields from input file ..." << endl;
     
      //
      // parse the initial profiles for the fields
      //
   
      GridFunctionFactory  gridFactory;
   
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
         gridFunction->assign( m_magneticField, dir, m_mesh, a_cur_time );
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
            gridFunction->assign( temporaryProfile, m_mesh, a_cur_time );
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
         gridFunction->assign( m_electricField, dir, m_mesh, a_cur_time );
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
            gridFunction->assign( temporaryNodeProfile, m_mesh, a_cur_time );
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
      applyBCs_electricField(a_cur_time);
      applyBCs_magneticField(a_cur_time);
   
   }   
   else {
   
      if(!procID()) cout << "Initializing electromagnetic fields from restart file..." << endl;

#ifdef CH_USE_HDF5
      HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );

      // read in the magnetic field data
      handle.setGroup("magnetic_field");
      read( handle, m_magneticField, "data", grids );
   
      // read in the electric field data
      handle.setGroup("electric_field");
      read(handle, m_electricField, "data", grids);
   
      // read in the virtual field data
      if(SpaceDim<3) {
         handle.setGroup("virtual_magnetic_field");
         read(handle, m_magneticField_virtual, "data", grids);
      
         handle.setGroup("virtual_electric_field");
         read(handle, m_electricField_virtual, "data", grids);
      }
  
      handle.close();

#else
      MayDay::Error("restart only defined with hdf5");
#endif
 
   }

   // set the initial divs
   setDivE();
   setDivB();

   // set the initial curls
   setCurlE();
   setCurlB();

   if(!procID()) cout << "Finished initializing electromagnetic fields " << endl << endl;

}

void ElectroMagneticFields::computeRHSMagneticField(  const Real a_time,
                                                      const Real a_cnormDt )
{
  CH_TIME("ElectroMagneticFields::computeRHSMagneticField()");
  
  setCurlE();

  SpaceUtils::zero( m_magneticField_rhs );
  if (SpaceDim<3) SpaceUtils::zero( m_magneticField_virtual_rhs );

  if (!advanceB()) return;

  auto grids(m_mesh.getDBL());

  if(m_advanceB_inPlane) {
    for(auto dit(grids.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir < SpaceDim; dir++) {
        if(m_advanceB_comp[dir]) {
           FArrayBox& Brhsdir = m_magneticField_rhs[dit][dir];
           Brhsdir.copy( m_curlE[dit][dir] );
           Brhsdir.mult(-a_cnormDt);
        }
      }
    }
  }

  if(SpaceDim<3 && m_advanceB_virtual) {
    for(auto dit(grids.dataIterator()); dit.ok(); ++dit) {
      Box cell_box = grids[dit];
      const int Ncomp = m_magneticField_virtual.nComp();
      for (int comp=0; comp<Ncomp; comp++){
        FArrayBox& Bvrhsdir  = m_magneticField_virtual_rhs[dit];
        Bvrhsdir.copy(m_curlE_virtual[dit], comp, comp);
        Bvrhsdir.mult(-a_cnormDt*m_advanceB_comp[SpaceDim+comp], comp);
      }
    }
  }

  return;
}

void ElectroMagneticFields::computeRHSElectricField(  const Real a_time, 
                                                      const Real a_cnormDt )
{
  CH_TIME("ElectroMagneticFields::computeRHSElectricField()");
    
  setCurlB();
  solvePoisson(a_time);
  
  SpaceUtils::zero( m_electricField_rhs );
  if (SpaceDim<3) SpaceUtils::zero( m_electricField_virtual_rhs );

  if (!advanceE()) return;

  auto grids(m_mesh.getDBL());

  if(m_advanceE_inPlane) {
    for(auto dit(grids.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
        if(m_advanceE_comp[dir]) {
          const FArrayBox& curlBdir = m_curlB[dit][dir];
          const FArrayBox& Jdir = m_currentDensity[dit][dir];
          FArrayBox& Erhsdir = m_electricField_rhs[dit][dir];

          Erhsdir.copy( Jdir );
          Erhsdir.mult( -m_Jnorm_factor );
          Erhsdir.plus( curlBdir );
          Erhsdir.mult( a_cnormDt );

        }
      }
    }
  }

  if(SpaceDim<3 && m_advanceE_virtual) {
     for(auto dit(grids.dataIterator()); dit.ok(); ++dit) {
        const int Ncomp = m_electricField_virtual.nComp();
        for (int comp=0; comp<Ncomp; comp++){
           if(m_advanceE_comp[SpaceDim+comp]) {
              const FArrayBox& curlBv = m_curlB_virtual[dit].getFab();
              const FArrayBox& Jv = m_currentDensity_virtual[dit].getFab();
              FArrayBox& Evrhsdir = m_electricField_virtual_rhs[dit].getFab();

              Evrhsdir.copy( Jv, comp, comp );
              Evrhsdir.mult( -m_Jnorm_factor, comp );
              Evrhsdir.plus( curlBv, comp, comp );
              Evrhsdir.mult( a_cnormDt, comp );
           }
        }
     }
  }

  return;
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

Real ElectroMagneticFields::electricFieldEnergy()
{
  CH_TIME("ElectroMagneticFields::electricFieldEnergy()");
  const DisjointBoxLayout& grids(m_mesh.getDBL());
  
  const LevelData<NodeFArrayBox>& masked_Ja_nc(m_mesh.getMaskedJnc());
  const LevelData<EdgeDataBox>& masked_Ja_ec(m_mesh.getMaskedJec());
 
  Real energyE_local = 0.0;
  for(DataIterator dit(grids); dit.ok(); ++dit) {
    for(int dir=0; dir<SpaceDim; dir++) {
      Box edge_box = grids[dit];
      edge_box.surroundingNodes();
      edge_box.enclosedCells(dir);
      
      const FArrayBox& this_E  = m_electricField[dit][dir];
      const FArrayBox& this_Ja = masked_Ja_ec[dit][dir];
      for (BoxIterator bit(edge_box); bit.ok(); ++bit) {
        energyE_local += this_E(bit(),0)*this_E(bit(),0)*this_Ja(bit(),0);
      }
    }
  }
  
  if(SpaceDim<3) {
    for(DataIterator dit(grids); dit.ok(); ++dit) {
      Box node_box = grids[dit];
      node_box.surroundingNodes();

      const FArrayBox& this_E  = m_electricField_virtual[dit].getFab();
      const FArrayBox& this_Ja = masked_Ja_nc[dit].getFab();
      for(int n=0; n<this_E.nComp(); n++) {
        for (BoxIterator bit(node_box); bit.ok(); ++bit) {
          energyE_local += this_E(bit(),n)*this_E(bit(),n)*this_Ja(bit(),0);
        }
      }
    }
  }
  
  energyE_local = energyE_local*m_energyE_factor; // [Joules]
  
  Real energyE_global = 0.0;
#ifdef CH_MPI
  MPI_Allreduce(  &energyE_local,
                  &energyE_global,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
  energyE_global = energyE_local;
#endif

  return energyE_global;
}

Real ElectroMagneticFields::magneticFieldEnergy()
{
  CH_TIME("ElectroMagneticFields::magneticFieldEnergy()");
  const DisjointBoxLayout& grids(m_mesh.getDBL());

  const LevelData<FArrayBox>& Ja_cc(m_mesh.getJcc());
  const LevelData<FluxBox>& masked_Ja_fc(m_mesh.getMaskedJfc());
 
  Real energyB_local = 0.0;
  for(DataIterator dit(grids); dit.ok(); ++dit) {
    for(int dir=0; dir<SpaceDim; dir++) {
      Box face_box = grids[dit];
      face_box.surroundingNodes(dir);
      
      const FArrayBox& this_B  = m_magneticField[dit][dir];
      const FArrayBox& this_Ja = masked_Ja_fc[dit][dir];
      for (BoxIterator bit(face_box); bit.ok(); ++bit) {
        energyB_local += this_B(bit(),0)*this_B(bit(),0)*this_Ja(bit(),0);
      }
    }
  }
  
  if(SpaceDim<3) {
    for(DataIterator dit(grids); dit.ok(); ++dit) {
      Box cell_box = grids[dit];

      const FArrayBox& this_B  = m_magneticField_virtual[dit];
      const FArrayBox& this_Ja = Ja_cc[dit];
      for(int n=0; n<this_B.nComp(); n++) {
        for (BoxIterator bit(cell_box); bit.ok(); ++bit) {
          energyB_local += this_B(bit(),n)*this_B(bit(),n)*this_Ja(bit(),0);
        }
      }
    }
  }
   
  energyB_local = energyB_local*m_energyB_factor; // [Joules]

  Real energyB_global = 0.0;
#ifdef CH_MPI
  MPI_Allreduce(  &energyB_local,
                  &energyB_global,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
  energyB_global = energyB_local;
#endif

  return energyB_global;

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

   if (!advanceE()) return;
   if (!usePoisson()) return;

   setDivE();

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

   correctElectricField();
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

