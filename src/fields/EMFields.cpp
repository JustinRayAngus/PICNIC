#include <array>
#include <cmath>
#include "PicnicConstants.H"
#include "BoxIterator.H"
#include "CH_HDF5.H"
#include "MathUtils.H"
#include "EMFields.H"
#include "FieldsF_F.H"

#include "NamespaceHeader.H"

EMFields::EMFields( ParmParse&   a_ppflds,
              const DomainGrid&  a_mesh,
              const CodeUnits&   a_units,
              const bool&        a_verbosity )
    : m_verbosity(a_verbosity),
      m_use_filtering(false),
      m_use_poisson(false),
      m_enforce_gauss_startup(false),
      m_external_fields(false),
      m_writeDivs(false),
      m_writeCurls(false),
      m_writeRho(false),
      m_writeSigma(false),
      m_writeExB(false),
      m_mesh(a_mesh),
      m_field_bc(NULL),
      m_vec_size_bfield(0),
      m_vec_size_bfield_virtual(0),
      m_vec_size_efield(0),
      m_vec_size_efield_virtual(0),
      m_vec_offset_bfield(-1),
      m_vec_offset_bfield_virtual(-1),
      m_vec_offset_efield(-1),
      m_vec_offset_efield_virtual(-1),
      m_pc_mass_matrix_width(1),
      m_pc_mass_matrix_include_ij(false),
      m_pc_diag_only(false),
      m_vec_type(e_and_b)
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

   a_ppflds.query( "write_divs", m_writeDivs );
   a_ppflds.query( "write_curls", m_writeCurls );
   a_ppflds.query( "write_rho", m_writeRho );
   a_ppflds.query( "write_sigma", m_writeSigma );
   a_ppflds.query( "write_exb", m_writeExB );

   a_ppflds.query( "pc_diagonal_only", m_pc_diag_only);
   a_ppflds.query( "pc_mass_matrix_width", m_pc_mass_matrix_width);
   if(m_pc_mass_matrix_width>0) {
      a_ppflds.query( "pc_mass_matrix_include_ij", m_pc_mass_matrix_include_ij);
   }

   a_ppflds.query( "use_filtering", m_use_filtering );
   if (m_use_filtering) {
     m_filterE_inPlane = m_advanceE_inPlane;
     m_filterE_virtual = m_advanceE_virtual;
     a_ppflds.query( "filterE_inPlane", m_filterE_inPlane );
     a_ppflds.query( "filterE_virtual", m_filterE_virtual );
   }
   else {
     m_filterE_inPlane = false;
     m_filterE_virtual = false;
   }

   // set the Courant time step associated with the speed of light
   const RealVect& meshSpacing(m_mesh.getdX());
   Real invDXsq = 0.0;
   for (int dir=0; dir<SpaceDim; dir++) {
      invDXsq = invDXsq + 1.0/meshSpacing[dir]/meshSpacing[dir];
   }
   Real DXeff = 1.0/pow(invDXsq,0.5);
   m_cvacNorm = a_units.CvacNorm();
   m_stable_dt = DXeff/m_cvacNorm;

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

   const Real volume_scale = m_mesh.getVolumeScale();
   const Real dVolume = m_mesh.getMappedCellVolume()*volume_scale; // SI units
   m_energyE_factor = 0.5*Escale*Escale*ep0*dVolume;
   m_energyB_factor = 0.5*Bscale*Bscale/mu0*dVolume;

   const Real area_scale = m_mesh.getAreaScale();
   m_SdA_factor = Escale*Bscale/mu0*area_scale;
   m_time_scale = a_units.getScale(a_units.TIME);

   // define profiles for external fields
   a_ppflds.query( "external_fields", m_external_fields );
   if(m_external_fields) initExternalFields();

   //
   // define the member LevelDatas
   //

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit;

   // define the magnetic field containers
   m_magneticField.define(grids,1,ghostVect); // FluxBox with 1 comp
   m_magneticField_old.define(grids,1);
   m_magneticField_rhs.define(grids,1);
   SpaceUtils::zero( m_magneticField_rhs );
   if(SpaceDim<3) {
      m_magneticField_virtual.define(grids,3-SpaceDim,ghostVect); // FArrayBox
      m_magneticField_virtual_old.define(grids,3-SpaceDim);
      m_magneticField_virtual_rhs.define(grids,3-SpaceDim);
      SpaceUtils::zero( m_magneticField_virtual_rhs );
   }

   // define the electric field containers
   m_electricField.define(grids,1,ghostVect); // EdgeDataBox with 1 comp
   m_electricField_old.define(grids,1,IntVect::Zero);
   m_electricField_rhs.define(grids,1,IntVect::Zero);
   SpaceUtils::zero( m_electricField_rhs );
   if(SpaceDim<3) {
      m_electricField_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
      m_electricField_virtual_old.define(grids,3-SpaceDim,IntVect::Zero);
      m_electricField_virtual_rhs.define(grids,3-SpaceDim,IntVect::Zero);
      SpaceUtils::zero( m_electricField_virtual_rhs );
   }

   // define the filtered field containers
   m_electricField_filtered.define(grids,1,ghostVect); // EdgeDataBox with 1 comp
   m_electricField_virtual_filtered.define(grids,3-SpaceDim,ghostVect);

   // define the current density containers
   m_currentDensity.define(grids,1,ghostVect);    // EdgeDataBox with 1 comp
   if(SpaceDim<3) {m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect);} // NodeFArrayBox
   m_chargeDensity.define(grids,1,IntVect::Zero);
   m_surfaceCharge.define(grids,1,IntVect::Zero);


   // define the curlE containers
   m_curlE.define(grids,1,IntVect::Zero);     // FluxBox with 1 comp
   if(SpaceDim<3) {m_curlE_virtual.define(grids,3-SpaceDim,IntVect::Zero);} // FArrayBox

   // define the curlB containers
   m_curlB.define(grids,1,IntVect::Zero);     // EdgeDataBox with 1 comp
   if(SpaceDim<3) {m_curlB_virtual.define(grids,3-SpaceDim,IntVect::Zero);} // NodeFArrayBox

   // define ExB containers (poynting flux)
   m_ExB.define(grids,2,IntVect::Zero);
   if(SpaceDim<3) {m_EvxB.define(grids,6-2*SpaceDim,IntVect::Zero);}

   // set in-plane curl to zero
   SpaceUtils::zero( m_curlE );
   SpaceUtils::zero( m_curlB );

   // define the divE and divB containers
   m_divE.define(grids,1,IntVect::Zero);
   m_divB.define(grids,1,IntVect::Zero);

   // define pointer to the BC object
   m_field_bc = new FieldBC( m_mesh, m_verbosity );

   // define poisson solver
   a_ppflds.query( "use_poisson", m_use_poisson );
   a_ppflds.query( "enforce_gauss_startup", m_enforce_gauss_startup );
   if(m_use_poisson || m_enforce_gauss_startup) {

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
      cout << " use binomial filtering = " << m_use_filtering << endl;
      if (m_use_filtering) {
        cout << " filterE_inPlane = " << (m_filterE_inPlane?"true":"false") << endl;
        cout << " filterE_virtual = " << (m_filterE_virtual?"true":"false") << endl;
      }
      cout << " use poisson = " << m_use_poisson << endl;
      cout << " enforce gauss startup = " << m_enforce_gauss_startup << endl;
      cout << " external fields = " << m_external_fields << endl;
      cout << " write curls = " << (m_writeCurls?"true":"false") << endl;
      cout << " write divs  = " << (m_writeDivs?"true":"false") << endl;
      cout << " write rho   = " << (m_writeRho?"true":"false") << endl;
      cout << " write sigma = " << (m_writeSigma?"true":"false") << endl;
      cout << " write ExB   = " << (m_writeExB?"true":"false") << endl;
      cout << " stable dt = " << m_stable_dt << endl << endl;
   }

   return;
}

void EMFields::defineVectorsAndDOFs( const EMVecType& a_vec_type)
{
   m_vec_type = a_vec_type;

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

   if ((a_vec_type == e_and_b) || (a_vec_type == e_only) || (a_vec_type == curl2)) {
     m_vec_offset_efield = vec_size_total;
     vec_size_total += m_vec_size_efield;

     if (SpaceDim < 3) {
       m_vec_offset_efield_virtual = vec_size_total;
       vec_size_total += m_vec_size_efield_virtual;
     }

   }

   if (a_vec_type == e_and_b) {
      m_gdofs.define(  vec_size_total,
                       m_magneticField,
                       m_magneticField_virtual,
                       m_electricField,
                       m_electricField_virtual );
   } else if ((a_vec_type == e_only) || (a_vec_type == curl2)) {
      m_gdofs.define(  vec_size_total,
                       m_electricField,
                       m_electricField_virtual );
   }

   // define and initialize BC masking containers

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit;

   if (a_vec_type == curl2) {

     m_PC_mask_E.define(grids,2+2*SpaceDim,ghostVect);
     for (int d=0; d < m_PC_mask_E.nComp(); d++) {
       SpaceUtils::setVal(m_PC_mask_E, d, 1.0);
     }
     if(SpaceDim<3) {
       m_PC_mask_Ev.define(grids,2+2*SpaceDim,ghostVect);
       for (int d=0; d < m_PC_mask_Ev.nComp(); d++) {
         SpaceUtils::setVal(m_PC_mask_Ev, d, 1.0);
       }
     }

   } else {

     m_PC_mask_B.define(grids,4,ghostVect);
     SpaceUtils::setVal(m_PC_mask_B, 0, 1.0);
     SpaceUtils::setVal(m_PC_mask_B, 1, 1.0);
     SpaceUtils::setVal(m_PC_mask_B, 2, 0.0);
     SpaceUtils::setVal(m_PC_mask_B, 3, 0.0);
     m_PC_mask_E.define(grids,4,ghostVect);
     SpaceUtils::setVal(m_PC_mask_E, 0, 1.0);
     SpaceUtils::setVal(m_PC_mask_E, 1, 1.0);
     SpaceUtils::setVal(m_PC_mask_E, 2, 0.0);
     SpaceUtils::setVal(m_PC_mask_E, 3, 0.0);
     if(SpaceDim<3) {
       m_PC_mask_Bv.define(grids,3*SpaceDim,ghostVect);
       m_PC_mask_Ev.define(grids,3*SpaceDim,ghostVect);
       for (int d=0; d < 2*SpaceDim; d++) {
         SpaceUtils::setVal(m_PC_mask_Bv, d, 1.0);
         SpaceUtils::setVal(m_PC_mask_Ev, d, 1.0);
       }
       for (int d=2*SpaceDim; d < 3*SpaceDim; d++) {
         SpaceUtils::setVal(m_PC_mask_Bv, d, 0.0);
         SpaceUtils::setVal(m_PC_mask_Ev, d, 0.0);
       }
     }


   }

   return;
}

EMFields::~EMFields()
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

void EMFields::initialize( const Real          a_cur_time,
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
         RefCountedPtr<GridFunction> gridFunction(NULL);
         if(ppBfieldIC.contains("type")) {
            gridFunction = gridFactory.create(ppBfieldIC,m_mesh,m_verbosity);
         }
         else {
            gridFunction = gridFactory.createZero();
         }
         gridFunction->assign( m_magneticField, dir, m_mesh, a_cur_time );
      }

      // set initial virtual magnetic field profile from ICs
#if CH_SPACEDIM<3
      int numVirt = 3-SpaceDim;
      for( int comp=0; comp<numVirt; comp++) {
         int field_dir = comp + SpaceDim;
         stringstream sBv;
         sBv << "IC.em_fields.magnetic." << field_dir;
         ParmParse ppBvfieldIC( sBv.str().c_str() );
         RefCountedPtr<GridFunction> gridFunction(NULL);
         if(ppBvfieldIC.contains("type")) {
            gridFunction = gridFactory.create(ppBvfieldIC,m_mesh,m_verbosity);
         }
         else {
            gridFunction = gridFactory.createZero();
         }
         gridFunction->assign( temporaryProfile, m_mesh, a_cur_time );
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            m_magneticField_virtual[dit].copy(temporaryProfile[dit],0,comp,1);
         }
      }
#endif

      // set initial electric field profile from ICs
      for (int dir=0; dir<SpaceDim; dir++) {
         stringstream sE;
         sE << "IC.em_fields.electric." << dir;
         ParmParse ppEfieldIC( sE.str().c_str() );
         RefCountedPtr<GridFunction> gridFunction(NULL);
         if(ppEfieldIC.contains("type")) {
            gridFunction = gridFactory.create(ppEfieldIC,m_mesh,m_verbosity);
         }
         else {
            gridFunction = gridFactory.createZero();
         }
         gridFunction->assign( m_electricField, dir, m_mesh, a_cur_time );
      }

      // set initial virtual electric field profile from ICs
#if CH_SPACEDIM<3
      for( int comp=0; comp<numVirt; comp++) {
         int field_dir = comp + SpaceDim;
         stringstream sEv;
         sEv << "IC.em_fields.electric." << field_dir;
         ParmParse ppEvfieldIC( sEv.str().c_str() );
         RefCountedPtr<GridFunction> gridFunction(NULL);
         if(ppEvfieldIC.contains("type")) {
            gridFunction = gridFactory.create(ppEvfieldIC,m_mesh,m_verbosity);
         }
         else {
            gridFunction = gridFactory.createZero();
         }
         gridFunction->assign( temporaryNodeProfile, m_mesh, a_cur_time );
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox& this_Ev  = m_electricField_virtual[dit].getFab();
            this_Ev.copy(temporaryNodeProfile[dit].getFab(),0,comp,1);
         }
      }
#endif

      SpaceUtils::exchangeFluxBox(m_magneticField);
      SpaceUtils::exchangeEdgeDataBox(m_electricField);
#if CH_SPACEDIM<3
      m_magneticField_virtual.exchange();
      SpaceUtils::exchangeNodeFArrayBox(m_electricField_virtual);
#endif

      // apply BCs to initial profiles
      applyBCs_electricField(a_cur_time);
      applyBCs_magneticField(a_cur_time);

      if(m_enforce_gauss_startup) {
        solvePoisson( a_cur_time );
        correctElectricField();
      }

      // initialize old values by copy
      for(DataIterator dit(grids); dit.ok(); ++dit) {

         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& Edir_old = m_electricField_old[dit][dir];
            const Box& edge_box = Edir_old.box();
            Edir_old.copy(m_electricField[dit][dir],edge_box);
            //
            FArrayBox& Bdir_old = m_magneticField_old[dit][dir];
            const Box& face_box = Bdir_old.box();
            Bdir_old.copy(m_magneticField[dit][dir],face_box);
         }

#if CH_SPACEDIM<3
         FArrayBox& Ev_old = m_electricField_virtual_old[dit].getFab();
         const Box& node_box = Ev_old.box();
         Ev_old.copy(m_electricField_virtual[dit].getFab(),node_box);
         //
         FArrayBox& Bv_old = m_magneticField_virtual_old[dit];
         const Box& cell_box = Bv_old.box();
         Bv_old.copy(m_magneticField_virtual[dit],cell_box);
#endif

      }

      for (int dir=0; dir<SpaceDim; dir++) {
         m_intSdAdt_lo[dir] = 0.0;
         m_intSdAdt_hi[dir] = 0.0;
      }

   }
   else {

      if(!procID()) cout << "Initializing electromagnetic fields from restart file..." << endl;

#ifdef CH_USE_HDF5
      HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );

      //
      //  read cummulative boundary diagnostic data
      //
      HDF5HeaderData header;
      handle.setGroup("field_bdry_data");
      header.readFromFile( handle );

      m_intSdAdt_lo[0] = header.m_real["intSdAdt_lo0"];
      m_intSdAdt_hi[0] = header.m_real["intSdAdt_hi0"];
#if CH_SPACEDIM>1
      m_intSdAdt_lo[1] = header.m_real["intSdAdt_lo1"];
      m_intSdAdt_hi[1] = header.m_real["intSdAdt_hi1"];
#endif
#if CH_SPACEDIM==3
      m_intSdAdt_lo[2] = header.m_real["intSdAdt_lo2"];
      m_intSdAdt_hi[2] = header.m_real["intSdAdt_hi2"];
#endif

      // read in the magnetic field data
      handle.setGroup("magnetic_field");
      LevelData<FluxBox> magField_temp;
      read(handle, magField_temp, "data", grids);
      magField_temp.copyTo(m_magneticField);

      // read in the electric field data
      handle.setGroup("electric_field");
      LevelData<EdgeDataBox> elecField_temp;
      read(handle, elecField_temp, "data", grids);
      elecField_temp.copyTo(m_electricField);

#if CH_SPACEDIM<3
      // read in the virtual magnetic field data
      handle.setGroup("virtual_magnetic_field");
      LevelData<FArrayBox> magField_virt_temp;
      read(handle, magField_virt_temp, "data", grids);
      magField_virt_temp.copyTo(m_magneticField_virtual);

      // read in the virtual electric field data
      handle.setGroup("virtual_electric_field");
      LevelData<NodeFArrayBox> elecField_virt_temp;
      read(handle, elecField_virt_temp, "data", grids);
      elecField_virt_temp.copyTo(m_electricField_virtual);
#endif

      // call exchange
      SpaceUtils::exchangeFluxBox(m_magneticField);
      SpaceUtils::exchangeEdgeDataBox(m_electricField);
#if CH_SPACEDIM<3
      m_magneticField_virtual.exchange();
      SpaceUtils::exchangeNodeFArrayBox(m_electricField_virtual);
#endif

      // apply BCs to initial profiles
      applyBCs_electricField(a_cur_time);
      applyBCs_magneticField(a_cur_time);

      //if(procID()==0) {SpaceUtils::inspectFluxBox(m_magneticField,0);}
      //if(procID()==0) {SpaceUtils::inspectEdgeDataBox(m_electricField,0);}
      //if(procID()==0) {SpaceUtils::inspectNodeFArrayBox(m_electricField_virtual,0);}
      //if(procID()==0) {SpaceUtils::inspectFArrayBox(m_magneticField_virtual,0,1);}

      //
      //  read in the old field data if it exist
      //  This data is needed for proper restart of time advance
      //  schemes that use staggering
      //

      // read in the old magnetic field data
      int group_id = handle.setGroup("magnetic_field_old");
      if(group_id<0) { // group not found, initialize by copy
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir=0; dir<SpaceDim; dir++) {
               FArrayBox& Bdir_old = m_magneticField_old[dit][dir];
               const Box& face_box = Bdir_old.box();
               Bdir_old.copy(m_magneticField[dit][dir],face_box);
            }
         }
      }
      else read(handle, m_magneticField_old, "data", grids);

      // read in the old electric field data
      group_id = handle.setGroup("electric_field_old");
      if(group_id<0) { // group not found, initialize by copy
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir=0; dir<SpaceDim; dir++) {
               FArrayBox& Edir_old = m_electricField_old[dit][dir];
               const Box& edge_box = Edir_old.box();
               Edir_old.copy(m_electricField[dit][dir],edge_box);
            }
         }
      }
      else read(handle, m_electricField_old, "data", grids);

#if CH_SPACEDIM<3
      // read in the old virtual field data
      group_id = handle.setGroup("virtual_magnetic_field_old");
      if(group_id<0) { // group not found, initialize by copy
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox& Bv_old = m_magneticField_virtual_old[dit];
            const Box& cell_box = Bv_old.box();
            Bv_old.copy(m_magneticField_virtual[dit],cell_box);
         }
      }
      else read(handle, m_magneticField_virtual_old, "data", grids);

      group_id = handle.setGroup("virtual_electric_field_old");
      if(group_id<0) { // group not found, initialize by copy
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox& Ev_old = m_electricField_virtual_old[dit].getFab();
            const Box& node_box = Ev_old.box();
            Ev_old.copy(m_electricField_virtual[dit].getFab(),node_box);
         }
      }
      else read(handle, m_electricField_virtual_old, "data", grids);
#endif

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

void EMFields::initializeMassMatricesForPC( const IntVect&  a_ncomp_xx,
		                            const IntVect&  a_ncomp_xy,
		                            const IntVect&  a_ncomp_xz,
		                            const IntVect&  a_ncomp_yx,
		                            const IntVect&  a_ncomp_yy,
		                            const IntVect&  a_ncomp_yz,
		                            const IntVect&  a_ncomp_zx,
		                            const IntVect&  a_ncomp_zy,
		                            const IntVect&  a_ncomp_zz )
{
   CH_TIME("EMFields::initializeMassMatricesForPC()");

   //////////////////////////////////////////////////////////////////
   //
   // This is called from PicSpeciesInterface::initializeMassMatrices()
   // The passed number of components are those corresponding to the full
   // mass matrices used to compute J. They are only used here to ensure
   // that the mass matrices used for the PC do not have more componets
   // than the mass matrices used to compute J
   //
   //////////////////////////////////////////////////////////////////

   if (!procID()) { cout << "Initializing mass matrices for PC..." << endl; }

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit;

   int ncomp_xx;
   if(m_use_filtering) ncomp_xx = 1 + 2*(m_pc_mass_matrix_width+2);
   else ncomp_xx = 1 + 2*m_pc_mass_matrix_width;
   for(int dir=0; dir<SpaceDim; dir++) {
      m_ncomp_xx_pc[dir] = std::min(ncomp_xx,a_ncomp_xx[dir]);
   }
   m_sigma_xx_pc.define(grids,m_ncomp_xx_pc.product(),ghostVect);
   SpaceUtils::zero( m_sigma_xx_pc );

   if(m_pc_mass_matrix_include_ij) {
     int ncomp_xy;
     if(m_use_filtering) ncomp_xy = 2*(m_pc_mass_matrix_width+2);
     else ncomp_xy = 2*m_pc_mass_matrix_width;
     m_ncomp_xy_pc[0] = ncomp_xy;
#if CH_SPACEDIM==2
     m_ncomp_xy_pc[1] = ncomp_xy;
#elif CH_SPACEDIM==3
     m_ncomp_xy_pc[2] = ncomp_xy + 1;
#endif
     for (int dir=0; dir<SpaceDim; dir++) {
       if(m_ncomp_xy_pc[dir]>=a_ncomp_xy[dir]) {
         m_ncomp_xy_pc = a_ncomp_xy;
         break;
       }
     }

     int ncomp_xz;
     if(m_use_filtering) ncomp_xz = 2*(m_pc_mass_matrix_width+2);
     else ncomp_xz = 2*m_pc_mass_matrix_width;
     m_ncomp_xz_pc[0] = ncomp_xz;
#if CH_SPACEDIM==2
     m_ncomp_xz_pc[1] = ncomp_xz + 1;
#elif CH_SPACEDIM==3
     m_ncomp_xz_pc[2] = ncomp_xz ;
#endif
     for (int dir=0; dir<SpaceDim; dir++) {
       if(m_ncomp_xz_pc[dir]>=a_ncomp_xz[dir]) {
         m_ncomp_xz_pc = a_ncomp_xz;
         break;
       }
     }
     m_sigma_xy_pc.define(grids,m_ncomp_xy_pc.product(),ghostVect);
     m_sigma_xz_pc.define(grids,m_ncomp_xz_pc.product(),ghostVect);
     SpaceUtils::zero( m_sigma_xy_pc );
     SpaceUtils::zero( m_sigma_xz_pc );
   }
   else{
     m_ncomp_xy_pc = IntVect::Zero;
     m_ncomp_xz_pc = IntVect::Zero;
   }

#if CH_SPACEDIM==1
   // define the PC mass matrices corresponding to Jy
   m_ncomp_yy_pc[0] = 1 + 2*m_pc_mass_matrix_width;
   if(m_ncomp_yy_pc[0]>a_ncomp_yy[0]) {m_ncomp_yy_pc[0] = a_ncomp_yy[0];}
   m_sigma_yy_pc.define(grids,m_ncomp_yy_pc[0],ghostVect);
   SpaceUtils::zero( m_sigma_yy_pc );

   if(m_pc_mass_matrix_include_ij) {
     m_ncomp_yx_pc[0] = 2*m_pc_mass_matrix_width;
     m_ncomp_yz_pc[0] = 1 + 2*m_pc_mass_matrix_width;
     if(m_ncomp_yx_pc[0]>a_ncomp_yx[0]) {m_ncomp_yx_pc = a_ncomp_yx;}
     if(m_ncomp_yz_pc[0]>a_ncomp_yz[0]) {m_ncomp_yz_pc = a_ncomp_yz;}
     m_sigma_yx_pc.define(grids,m_ncomp_yx_pc[0],ghostVect);
     m_sigma_yz_pc.define(grids,m_ncomp_yz_pc[0],ghostVect);
     SpaceUtils::zero( m_sigma_yx_pc );
     SpaceUtils::zero( m_sigma_yz_pc );
   }
   else{
     m_ncomp_yx_pc = IntVect::Zero;
     m_ncomp_yz_pc = IntVect::Zero;
   }
#elif CH_SPACEDIM==2
   int ncomp_yy = 1 + 2*m_pc_mass_matrix_width;
   for(int dir=0; dir<SpaceDim; dir++) {
      m_ncomp_yy_pc[dir] = std::min(ncomp_yy,a_ncomp_yy[dir]);
   }

   if(m_pc_mass_matrix_include_ij) {
     m_ncomp_yx_pc[0] = 2*m_pc_mass_matrix_width;
     m_ncomp_yx_pc[1] = 2*m_pc_mass_matrix_width;
     m_ncomp_yz_pc[0] = 1 + 2*m_pc_mass_matrix_width;
     m_ncomp_yz_pc[1] = 2*m_pc_mass_matrix_width;
     for (int dir=0; dir<SpaceDim; dir++) {
       if(m_ncomp_yx_pc[dir]>=a_ncomp_yx[dir]) {
         m_ncomp_yx_pc = a_ncomp_yx;
         break;
       }
     }
     for (int dir=0; dir<SpaceDim; dir++) {
       if(m_ncomp_yz_pc[dir]>=a_ncomp_yz[dir]) {
         m_ncomp_yz_pc = a_ncomp_yz;
         break;
       }
     }
   }
   else{
     m_ncomp_yx_pc = IntVect::Zero;
     m_ncomp_yz_pc = IntVect::Zero;
   }
#endif

   // define the PC mass matrices corresponding to Jz
   int ncomp_zz = 1 + 2*m_pc_mass_matrix_width;
   for(int dir=0; dir<SpaceDim; dir++) { m_ncomp_zz_pc[dir] = ncomp_zz; }
   if(m_ncomp_zz_pc[0]>a_ncomp_zz[0]) {m_ncomp_zz_pc = a_ncomp_zz;}
   m_sigma_zz_pc.define(grids,m_ncomp_zz_pc.product(),ghostVect);
   SpaceUtils::zero( m_sigma_zz_pc );

   if(m_pc_mass_matrix_include_ij) {
     m_ncomp_zx_pc[0] = 2*m_pc_mass_matrix_width;
     m_ncomp_zy_pc[0] = 1+2*m_pc_mass_matrix_width;
#if CH_SPACEDIM==2
      m_ncomp_zx_pc[1] = 1+2*m_pc_mass_matrix_width;
     m_ncomp_zy_pc[1] = 2*m_pc_mass_matrix_width;
#elif CH_SPACEDIM==3
     m_ncomp_zx_pc[2] = 2*m_pc_mass_matrix_width;
     m_ncomp_zy_pc[2] = 2*m_pc_mass_matrix_width;
#endif
     for (int dir=0; dir<SpaceDim; dir++) {
       if(m_ncomp_zx_pc[dir]>=a_ncomp_zx[dir]) {
         m_ncomp_zx_pc = a_ncomp_zx;
         break;
       }
     }
     for (int dir=0; dir<SpaceDim; dir++) {
       if(m_ncomp_zy_pc[dir]>=a_ncomp_zy[dir]) {
         m_ncomp_zy_pc = a_ncomp_zy;
         break;
       }
     }
     m_sigma_zx_pc.define(grids,m_ncomp_zx_pc.product(),ghostVect);
     m_sigma_zy_pc.define(grids,m_ncomp_zy_pc.product(),ghostVect);
     SpaceUtils::zero( m_sigma_zx_pc );
     SpaceUtils::zero( m_sigma_zy_pc );
   }
   else{
     m_ncomp_zx_pc = IntVect::Zero;
     m_ncomp_zy_pc = IntVect::Zero;
   }

   if(!procID()) {
      if(m_use_filtering) cout << " filtering is being used ... " << endl;
      cout << " pc_mass_matrix_width = " << m_pc_mass_matrix_width << endl;
      cout << " pc_mass_matrix_include_ij = " << m_pc_mass_matrix_include_ij << endl;
      cout << " ncomp_xx_pc = " << m_ncomp_xx_pc << endl;
      cout << " ncomp_xy_pc = " << m_ncomp_xy_pc << endl;
      cout << " ncomp_xz_pc = " << m_ncomp_xz_pc << endl;
      cout << " ncomp_yx_pc = " << m_ncomp_yx_pc << endl;
      cout << " ncomp_yy_pc = " << m_ncomp_yy_pc << endl;
      cout << " ncomp_yz_pc = " << m_ncomp_yz_pc << endl;
      cout << " ncomp_zx_pc = " << m_ncomp_zx_pc << endl;
      cout << " ncomp_zy_pc = " << m_ncomp_zy_pc << endl;
      cout << " ncomp_zz_pc = " << m_ncomp_zz_pc << endl;
      cout << "Finished initializing mass matrices for PC " << endl << endl;
#if CH_SPACEDIM==1
      if (m_use_filtering && m_pc_mass_matrix_width<3) {
         cout << "Warning: filtering being used and pc_mass_matrix_width = " << m_pc_mass_matrix_width << endl;
         cout << "         Consider setting pc_mass_matrix_width = 3 and " << endl;
         cout << "         num_ghosts = 3 if GMRES is converging slow or not at all. " << endl << endl;
      }
#endif
   }

}


void EMFields::initExternalFields()
{

   // define grid function describing external fields

   GridFunctionFactory  gridFactory;

   //
   //  define external electric field profiles
   //

   stringstream sE0;
   sE0 << "em_fields.external.electric.0";
   ParmParse ppExtE0( sE0.str().c_str() );
   m_gridFunction_extE0 = gridFactory.create(ppExtE0,m_mesh,m_verbosity);

   stringstream sE1;
   sE1 << "em_fields.external.electric.1";
   ParmParse ppExtE1( sE1.str().c_str() );
   m_gridFunction_extE1 = gridFactory.create(ppExtE1,m_mesh,m_verbosity);

   stringstream sE2;
   sE2 << "em_fields.external.electric.2";
   ParmParse ppExtE2( sE2.str().c_str() );
   m_gridFunction_extE2 = gridFactory.create(ppExtE2,m_mesh,m_verbosity);

   //
   //  define external magnetic field profiles
   //

   stringstream sB0;
   sB0 << "em_fields.external.magnetic.0";
   ParmParse ppExtB0( sB0.str().c_str() );
   m_gridFunction_extB0 = gridFactory.create(ppExtB0,m_mesh,m_verbosity);

   stringstream sB1;
   sB1 << "em_fields.external.magnetic.1";
   ParmParse ppExtB1( sB1.str().c_str() );
   m_gridFunction_extB1 = gridFactory.create(ppExtB1,m_mesh,m_verbosity);

   stringstream sB2;
   sB2 << "em_fields.external.magnetic.2";
   ParmParse ppExtB2( sB2.str().c_str() );
   m_gridFunction_extB2 = gridFactory.create(ppExtB2,m_mesh,m_verbosity);

}

void EMFields::computeRHSMagneticField( const Real  a_dt )
{
  CH_TIME("EMFields::computeRHSMagneticField()");

  setCurlE();

  if (!advanceB()) return;

  const Real cnormDt = a_dt*m_cvacNorm;
  auto grids(m_mesh.getDBL());

  if(m_advanceB_inPlane) {
    SpaceUtils::zero( m_magneticField_rhs );
    for(auto dit(grids.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir < SpaceDim; dir++) {
        if(m_advanceB_comp[dir]) {
           FArrayBox& Brhsdir = m_magneticField_rhs[dit][dir];
           Brhsdir.copy( m_curlE[dit][dir] );
           Brhsdir.mult(-cnormDt);
        }
      }
    }
  }

#if CH_SPACEDIM<3
  if(m_advanceB_virtual) {
    SpaceUtils::zero( m_magneticField_virtual_rhs );
    for(auto dit(grids.dataIterator()); dit.ok(); ++dit) {
      Box cell_box = grids[dit];
      const int Ncomp = m_magneticField_virtual.nComp();
      for (int comp=0; comp<Ncomp; comp++){
        FArrayBox& Bvrhsdir  = m_magneticField_virtual_rhs[dit];
        Bvrhsdir.copy(m_curlE_virtual[dit], comp, comp);
        Bvrhsdir.mult(-cnormDt*m_advanceB_comp[SpaceDim+comp], comp);
      }
    }
  }
#endif

  return;
}

void EMFields::computeRHSElectricField( const Real  a_dt )
{
  CH_TIME("EMFields::computeRHSElectricField()");

  setCurlB();

  if(!advanceE()) return;

  const Real cnormDt = a_dt*m_cvacNorm;
  auto grids(m_mesh.getDBL());

  if(m_advanceE_inPlane) {
    SpaceUtils::zero( m_electricField_rhs );
    for(auto dit(grids.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
        if(m_advanceE_comp[dir]) {
          const FArrayBox& curlBdir = m_curlB[dit][dir];
          const FArrayBox& Jdir = m_currentDensity[dit][dir];
          FArrayBox& Erhsdir = m_electricField_rhs[dit][dir];

          Erhsdir.copy( Jdir );
          Erhsdir.mult( -m_Jnorm_factor );
          Erhsdir.plus( curlBdir );
          Erhsdir.mult( cnormDt );
        }
      }
    }
  }

#if CH_SPACEDIM<3
  if(m_advanceE_virtual) {
    SpaceUtils::zero( m_electricField_virtual_rhs );
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
          Evrhsdir.mult( cnormDt, comp );
        }
      }
    }
  }
#endif

  return;
}

void EMFields::applyBCs_electricField( const Real  a_time )
{
   CH_TIME("EMFields::applyBCs_electricField()");

   m_field_bc->applyEdgeBC( m_electricField, a_time );
   if (m_PC_mask_E.isDefined()) {
      m_field_bc->applyEdgePCMask( m_PC_mask_E, a_time, m_vec_type );
      SpaceUtils::exchangeEdgeDataBox( m_PC_mask_E );
   }
#if CH_SPACEDIM<3
   m_field_bc->applyNodeBC( m_electricField_virtual, a_time );
   if (m_PC_mask_Ev.isDefined()) {
     m_field_bc->applyNodePCMask( m_PC_mask_Ev, a_time, m_vec_type );
     SpaceUtils::exchangeNodeFArrayBox( m_PC_mask_Ev );
   }
#endif

}

void EMFields::applyBCs_magneticField( const Real  a_time )
{
   CH_TIME("EMFields::applyBCs_magneticField()");

   m_field_bc->applyFluxBC( m_magneticField, a_time );
   if (m_PC_mask_B.isDefined()) {
     m_field_bc->applyFluxPCMask( m_PC_mask_B, a_time );
     SpaceUtils::exchangeFluxBox( m_PC_mask_B );
   }
#if CH_SPACEDIM<3
   m_field_bc->applyCellBC( m_magneticField_virtual, a_time );
   if (m_PC_mask_B.isDefined()) {
     m_field_bc->applyCellPCMask( m_PC_mask_Bv, a_time );
     m_PC_mask_Bv.exchange();
   }
#endif

}

void EMFields::applyAbsorbingBCs( const Real  a_time,
                                  const Real  a_dt )
{
    CH_TIME("EMFields::applyAbsorbingBCs()");

#if CH_SPACEDIM==1
    const Real cnormDt = a_dt*m_cvacNorm;
    m_field_bc->applyAbsorbingBCs( m_magneticField_virtual,
                                   m_electricField_virtual_old,
                                   m_currentDensity_virtual, m_Jnorm_factor,
                                   a_time, cnormDt );
#endif

}

void EMFields::updateBoundaryProbes( const Real  a_dt )
{
   CH_TIME("EMFields::updateBoundaryProbes()");

   if(m_field_bc==NULL) return;

   // update surface integrals of poynting flux for diagnostics
   // by default, value for periodic BCs will be zero

   RealVect intSdA_lo_local(D_DECL6(0,0,0,0,0,0));
   RealVect intSdA_hi_local(D_DECL6(0,0,0,0,0,0));
   m_field_bc->computeIntSdA( intSdA_lo_local, intSdA_hi_local,
                              m_electricField,
                              m_magneticField,
                              m_electricField_virtual,
                              m_magneticField_virtual );

   intSdA_lo_local = intSdA_lo_local*m_SdA_factor; // [Joules/s]
   intSdA_hi_local = intSdA_hi_local*m_SdA_factor; // [Joules/s]
   m_intSdA_lo = intSdA_lo_local;
   m_intSdA_hi = intSdA_hi_local;
#ifdef CH_MPI
   MPI_Allreduce( &intSdA_lo_local,
                  &m_intSdA_lo,
                  SpaceDim,
                  MPI_CH_REAL,
                  MPI_SUM,
                  MPI_COMM_WORLD );
   MPI_Allreduce( &intSdA_hi_local,
                  &m_intSdA_hi,
                  SpaceDim,
                  MPI_CH_REAL,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#endif

   // update running integral over time
   for (int dir=0; dir<SpaceDim; ++dir) {
      m_intSdAdt_lo[dir] += m_intSdA_lo[dir]*a_dt*m_time_scale; // [Joules]
      m_intSdAdt_hi[dir] += m_intSdA_hi[dir]*a_dt*m_time_scale; // [Joules]
   }

}

void EMFields::setPoyntingFlux()
{
   CH_TIME("EMFields::setPoyntingFlux()");

   const DisjointBoxLayout& grids(m_mesh.getDBL());

   if(SpaceDim==3) {

      // worry about it later

   }

   if(SpaceDim==2) {

      for(DataIterator dit(grids); dit.ok(); ++dit) {

               FArrayBox& ExB_0 = m_ExB[dit][0];
               FArrayBox& ExB_1 = m_ExB[dit][1];
               FArrayBox& ExB_v = m_EvxB[dit].getFab();
         const FArrayBox& E0 = m_electricField[dit][0];
         const FArrayBox& E1 = m_electricField[dit][1];
         const FArrayBox& Ev = m_electricField_virtual[dit].getFab();
         const FArrayBox& B0 = m_magneticField[dit][0];
         const FArrayBox& B1 = m_magneticField[dit][1];
         const FArrayBox& Bv = m_magneticField_virtual[dit];

         Box node_box = surroundingNodes(grids[dit]);
         Box edge_box_0 = enclosedCells(node_box,0);
         Box edge_box_1 = enclosedCells(node_box,1);

         // ExB_0(0) = Ex*By (or -Ex*Bz)
         ExB_0.copy(E0,edge_box_0,0,edge_box_0,0,1);
         ExB_0.mult(B1,edge_box_0,0,0,1);
	 if(m_mesh.anticyclic()) ExB_0.negate(0,1);

         // ExB_0(1) = -Ex*Bz (or Ex*By)
         SpaceUtils::interpolateStag(ExB_0,edge_box_0,1,Bv,0,1);
         ExB_0.mult(E0,edge_box_0,0,1,1);
	 if(!m_mesh.anticyclic()) ExB_0.negate(1,1);

         // ExB_1(0) = Ey*Bz (or -Ez*By)
         SpaceUtils::interpolateStag(ExB_1,edge_box_1,0,Bv,0,0);
         ExB_1.mult(E1,edge_box_1,0,0,1);
	 if(m_mesh.anticyclic()) ExB_1.negate(0,1);

         // ExB_1(1) = -Ey*Bx (or Ez*Bx)
         ExB_1.copy(B0,edge_box_1,0,edge_box_1,1,1);
         ExB_1.mult(E1,edge_box_1,0,1,1);
	 if(!m_mesh.anticyclic()) ExB_1.negate(1,1);

         // ExB_v(0) = Ez*Bx (or -Ey*Bx)
         SpaceUtils::interpolateStag(ExB_v,node_box,0,B0,0,1);
         ExB_v.mult(Ev,node_box,0,0,1);
         if(m_mesh.anticyclic()) ExB_v.negate(0,1);

         // ExB_v(1) = -Ez*By (or Ey*Bz)
         SpaceUtils::interpolateStag(ExB_v,node_box,1,B1,0,0);
         ExB_v.mult(Ev,node_box,0,1,1);
         if(!m_mesh.anticyclic()) ExB_v.negate(1,1);

      }

   }

   if(SpaceDim==1) {

      for(DataIterator dit(grids); dit.ok(); ++dit) {

               FArrayBox& ExB = m_ExB[dit][0];
               FArrayBox& ExB_v = m_EvxB[dit].getFab();
         const FArrayBox& E0 = m_electricField[dit][0];
         const FArrayBox& Ev = m_electricField_virtual[dit].getFab();
         const FArrayBox& B0 = m_magneticField[dit][0];
         const FArrayBox& Bv = m_magneticField_virtual[dit];

         Box cell_box = grids[dit];
         Box node_box = surroundingNodes(grids[dit]);

         // ExB(0) = Ex*By
         ExB.copy(E0,cell_box,0,cell_box,0,1);
         ExB.mult(Bv,cell_box,0,0,1);

         // ExB(1) = -Ex*Bz
         ExB.copy(E0,cell_box,0,cell_box,1,1);
         ExB.mult(Bv,cell_box,1,1,1);
         ExB.negate(1,1);

         // ExB_v(0) = Ey*Bz
         SpaceUtils::interpolateStag(ExB_v,node_box,0,Bv,1,0);
         ExB_v.mult(Ev,node_box,0,0,1);

         // ExB_v(1) = -Ey*Bx
         ExB_v.copy(Ev,node_box,0,node_box,1,1);
         ExB_v.mult(B0,node_box,0,1,1);
         ExB_v.negate(1,1);

         // ExB_v(2) = Ez*Bx
         ExB_v.copy(Ev,node_box,1,node_box,2,1);
         ExB_v.mult(B0,node_box,0,2,1);

         // ExB_v(3) = -Ez*By
         SpaceUtils::interpolateStag(ExB_v,node_box,3,Bv,0,0);
         ExB_v.mult(Ev,node_box,1,3,1);
         ExB_v.negate(3,1);

      }

   }

}

Real EMFields::fieldEnergyMod( const Real  a_dt ) const
{
  CH_TIME("EMFields::fieldEnergyMod()");
  const DisjointBoxLayout& grids(m_mesh.getDBL());

  Real energyE_local = 0.0;
  const Real cnormDt = a_dt*m_cvacNorm;
  const Real factor = cnormDt*m_Jnorm_factor;

  const LevelData<EdgeDataBox>& masked_Ja_ec(m_mesh.getMaskedJec());
  for(DataIterator dit(grids); dit.ok(); ++dit) {
    for(int dir=0; dir<SpaceDim; dir++) {
      Box edge_box = grids[dit];
      edge_box.surroundingNodes();
      edge_box.enclosedCells(dir);

      const FArrayBox& this_E  = m_electricField[dit][dir];
      const FArrayBox& this_J  = m_currentDensity[dit][dir];
      const FArrayBox& this_Ja = masked_Ja_ec[dit][dir];
      for (BoxIterator bit(edge_box); bit.ok(); ++bit) {
        energyE_local += this_E(bit(),0)*this_J(bit(),0)*this_Ja(bit(),0);
      }
    }
  }

#if CH_SPACEDIM<3
  const LevelData<NodeFArrayBox>& masked_Ja_nc(m_mesh.getMaskedJnc());
  for(DataIterator dit(grids); dit.ok(); ++dit) {
    Box node_box = grids[dit];
    node_box.surroundingNodes();

    const FArrayBox& this_E  = m_electricField_virtual[dit].getFab();
    const FArrayBox& this_J  = m_currentDensity_virtual[dit].getFab();
    const FArrayBox& this_Ja = masked_Ja_nc[dit].getFab();
    for(int n=0; n<this_E.nComp(); n++) {
      for (BoxIterator bit(node_box); bit.ok(); ++bit) {
        energyE_local += this_E(bit(),n)*this_J(bit(),n)*this_Ja(bit(),0);
      }
    }
  }
#endif

  energyE_local = energyE_local*factor*m_energyE_factor; // [Joules]

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

Real EMFields::electricFieldEnergy( const bool  a_is_stag ) const
{
  CH_TIME("EMFields::electricFieldEnergy()");
  const DisjointBoxLayout& grids(m_mesh.getDBL());

  Real energyE_local = 0.0;

  const LevelData<EdgeDataBox>& masked_Ja_ec(m_mesh.getMaskedJec());
  for(DataIterator dit(grids); dit.ok(); ++dit) {
    for(int dir=0; dir<SpaceDim; dir++) {
      Box edge_box = grids[dit];
      edge_box.surroundingNodes();
      edge_box.enclosedCells(dir);

      const FArrayBox& this_Ja = masked_Ja_ec[dit][dir];
      const FArrayBox& this_E  = m_electricField[dit][dir];
      const FArrayBox& this_Eold  = m_electricField_old[dit][dir];
      for (BoxIterator bit(edge_box); bit.ok(); ++bit) {
        if (a_is_stag) { // E lives at half steps
          energyE_local += this_E(bit(),0)*this_Eold(bit(),0)*this_Ja(bit(),0);
        }
        else {
          energyE_local += this_E(bit(),0)*this_E(bit(),0)*this_Ja(bit(),0);
        }
      }
    }
  }

#if CH_SPACEDIM<3
  const LevelData<NodeFArrayBox>& masked_Ja_nc(m_mesh.getMaskedJnc());
  for(DataIterator dit(grids); dit.ok(); ++dit) {
    Box node_box = grids[dit];
    node_box.surroundingNodes();

    const FArrayBox& this_Ja = masked_Ja_nc[dit].getFab();
    const FArrayBox& this_E  = m_electricField_virtual[dit].getFab();
    const FArrayBox& this_Eold  = m_electricField_virtual_old[dit].getFab();
    for(int n=0; n<this_E.nComp(); n++) {
      for (BoxIterator bit(node_box); bit.ok(); ++bit) {
        if (a_is_stag) { // E lives at half steps
          energyE_local += this_E(bit(),n)*this_Eold(bit(),n)*this_Ja(bit(),0);
        }
        else {
          energyE_local += this_E(bit(),n)*this_E(bit(),n)*this_Ja(bit(),0);
        }
      }
    }
  }
#endif

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

Real EMFields::magneticFieldEnergy( const bool  a_is_stag ) const
{
  CH_TIME("EMFields::magneticFieldEnergy()");
  const DisjointBoxLayout& grids(m_mesh.getDBL());

  Real energyB_local = 0.0;

  const LevelData<FluxBox>& masked_Ja_fc(m_mesh.getMaskedJfc());
  for(DataIterator dit(grids); dit.ok(); ++dit) {
    for(int dir=0; dir<SpaceDim; dir++) {
      Box face_box = grids[dit];
      face_box.surroundingNodes(dir);
      const FArrayBox& this_Ja = masked_Ja_fc[dit][dir];
      const FArrayBox& this_B  = m_magneticField[dit][dir];
      const FArrayBox& this_Bold  = m_magneticField_old[dit][dir];
      for (BoxIterator bit(face_box); bit.ok(); ++bit) {
        if (a_is_stag) { // B lives at half steps
          energyB_local += this_B(bit(),0)*this_Bold(bit(),0)*this_Ja(bit(),0);
        }
        else {
          energyB_local += this_B(bit(),0)*this_B(bit(),0)*this_Ja(bit(),0);
        }
      }
    }
  }

#if CH_SPACEDIM<3
  const LevelData<FArrayBox>& Ja_cc(m_mesh.getJcc());
  for(DataIterator dit(grids); dit.ok(); ++dit) {
    Box cell_box = grids[dit];

    const FArrayBox& this_Ja = Ja_cc[dit];
    const FArrayBox& this_B  = m_magneticField_virtual[dit];
    const FArrayBox& this_Bold  = m_magneticField_virtual_old[dit];
    for(int n=0; n<this_B.nComp(); n++) {
      for (BoxIterator bit(cell_box); bit.ok(); ++bit) {
        if (a_is_stag) { // B lives at half steps
          energyB_local += this_B(bit(),n)*this_Bold(bit(),n)*this_Ja(bit(),0);
        }
        else {
          energyB_local += this_B(bit(),n)*this_B(bit(),n)*this_Ja(bit(),0);
        }
      }
    }
  }
#endif

  energyB_local = energyB_local*m_energyB_factor; // [Joules]

  Real energyB_global = 0.0;
#ifdef CH_MPI
  MPI_Allreduce(  &energyB_local,
                  &energyB_global,
                  1,
                  MPI_CH_REAL,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
  energyB_global = energyB_local;
#endif

  return energyB_global;

}

Real EMFields::max_wc0dt( const CodeUnits&  a_units,
                          const Real&       a_dt ) const
{
   Real max_absB_local = 0.0;

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for(int dir=0; dir<SpaceDim; dir++) {
         Box face_box = grids[dit];
         face_box.surroundingNodes(dir);
         Real box_max = m_magneticField[dit][dir].max(face_box,0);
         Real box_min = m_magneticField[dit][dir].min(face_box,0);
         box_max = Max(box_max,abs(box_min));
         max_absB_local = Max(max_absB_local,box_max);
      }
      if(SpaceDim<3) {
         Box cell_box = grids[dit];
         for(int n=0; n<m_magneticField_virtual.nComp(); n++) {
            Real box_max = m_magneticField_virtual[dit].max(cell_box,n);
            Real box_min = m_magneticField_virtual[dit].min(cell_box,n);
            box_max = Max(box_max,abs(box_min));
            max_absB_local = Max(max_absB_local,box_max);
         }
      }
   }

   Real max_absB_global = max_absB_local;
#ifdef CH_MPI
   MPI_Allreduce( &max_absB_local,
                  &max_absB_global,
                  1,
                  MPI_CH_REAL,
                  MPI_MAX,
                  MPI_COMM_WORLD );
#endif

  const Real dt_sec = a_dt*a_units.getScale(a_units.TIME);
  const Real B_T = max_absB_global*a_units.getScale(a_units.MAGNETIC_FIELD);

  const Real max_wc0dt = B_T*Constants::QE/Constants::ME*dt_sec;
  return max_wc0dt;

}

void EMFields::setCurlB()
{
   CH_TIME("EMFields::setCurlB()");

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const RealVect& dX(m_mesh.getdX());

#if CH_SPACEDIM==3

   for(DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; dir++) {

         int dirj = dir + 1;
         dirj = dirj % SpaceDim;
         int dirk = dir + 2;
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

   SpaceUtils::exchangeEdgeDataBox(m_curlB);

#elif CH_SPACEDIM==2

   const LevelData<FArrayBox>& Jcc = m_mesh.getJcc();
   const LevelData<EdgeDataBox>& Jec = m_mesh.getJec();
   const bool anticyclic = m_mesh.anticyclic();

   for(DataIterator dit(grids); dit.ok(); ++dit) {

      // FluxBox ==> NodeFArrayBox

      Box node_box = grids[dit];

      int dirj = 0;
      int dirk = 1;

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
      if(anticyclic) curlBv.negate();

      // FArrayBox ==> EdgeDataBox

      for (int dir=0; dir<SpaceDim; dir++) {

         int dirj = dir + 1;
         dirj = dirj % SpaceDim;

         FArrayBox& curlB_dir = m_curlB[dit][dir];
         const FArrayBox& B_virt = m_magneticField_virtual[dit];

         int additive = 0;
         int sign = 1-2*dir;
         if(anticyclic) sign = -sign;
         if(m_mesh.axisymmetric() && dir==1) {
            FArrayBox JaBv;
            JaBv.define(B_virt.box(),1);
            JaBv.copy(B_virt,0,0,1);
            JaBv.mult(Jcc[dit],0,0,1);
            SpaceUtils::simpleStagGradComp( curlB_dir, grids[dit], 0,
                                            JaBv, 0,
                                            sign*dX[dirj], dirj,
                                            additive );
            curlB_dir.divide(Jec[dit][dir],0,0,1); // NEED TO APPLY ON-AXIS CORRECTION
         }
         else {
            SpaceUtils::simpleStagGradComp( curlB_dir, grids[dit], 0,
                                            B_virt, 0,
                                            sign*dX[dirj], dirj,
                                            additive );
         }

      }

   }

   m_field_bc->applyOnAxisCurlBC( m_curlB, m_magneticField_virtual );
   SpaceUtils::exchangeEdgeDataBox(m_curlB);
   SpaceUtils::exchangeNodeFArrayBox(m_curlB_virtual);

#elif CH_SPACEDIM==1

   const string& geom_type = m_mesh.geomType();
   const LevelData<FArrayBox>& Xcc = m_mesh.getXcc();
   const LevelData<NodeFArrayBox>& Xnc = m_mesh.getXnc();

   for(DataIterator dit(grids); dit.ok(); ++dit) {

      int dirj = 0;
      int dirk = 1;

      FArrayBox& curlB_virt = m_curlB_virtual[dit].getFab();
      const FArrayBox& B_virt = m_magneticField_virtual[dit];

      int additive = 0;
      if(m_mesh.axisymmetric() && geom_type=="sph_R") {
         FArrayBox rBz(B_virt.box(),1);
         rBz.copy(B_virt,dirk,0,1);
         rBz.mult(Xcc[dit],0,0,1);
         SpaceUtils::simpleStagGradComp( curlB_virt, grids[dit], dirj,
                                         rBz, 0,
                                         -dX[0], 0,
                                         additive );
         curlB_virt.divide(Xnc[dit].getFab(),0,dirj,1); // NEED TO APPLY ON AXIS CORRECTION
      }
      else {
         SpaceUtils::simpleStagGradComp( curlB_virt, grids[dit], dirj,
                                         B_virt, dirk,
                                         -dX[0], 0,
                                         additive );
      }

      if(m_mesh.axisymmetric()) {
         FArrayBox rBy(B_virt.box(),1);
         rBy.copy(B_virt,dirj,0,1);
         rBy.mult(Xcc[dit],0,0,1);
         SpaceUtils::simpleStagGradComp( curlB_virt, grids[dit], dirk,
                                         rBy, 0,
                                         dX[0], 0,
                                         additive );
         curlB_virt.divide(Xnc[dit].getFab(),0,dirk,1); // NEED TO APPLY ON AXIS CORRECTION
      }
      else {
         SpaceUtils::simpleStagGradComp( curlB_virt, grids[dit], dirk,
                                         B_virt, dirj,
                                         dX[0], 0,
                                         additive );
      }

   }

   m_field_bc->applyOnAxisCurlBC( m_curlB_virtual, m_magneticField_virtual );
   SpaceUtils::exchangeEdgeDataBox(m_curlB);
   SpaceUtils::exchangeNodeFArrayBox(m_curlB_virtual);

#endif

}

void EMFields::setCurlE()
{
   CH_TIME("EMFields::setCurlE()");

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const RealVect& dX(m_mesh.getdX());

#if CH_SPACEDIM==3 // EdgeDataBox ==> FluxBox

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

#elif CH_SPACEDIM==2 // EdgeDataBox ==> FArrayBox and NodeFArrayBox ==> FluxBox

    const bool anticyclic = m_mesh.anticyclic();
    const LevelData<NodeFArrayBox>& Jnc = m_mesh.getJnc();
    const LevelData<FluxBox>& Jfc = m_mesh.getJfc();

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
      if(anticyclic) curlEv.negate();

      //
      // NodeFArrayBox ==> FluxBox
      //

      for (int dir=0; dir<SpaceDim; dir++) {

         int dirj;
         dirj = dir + 1;
         dirj = dirj % SpaceDim;

         FArrayBox& curlE_dir = m_curlE[dit][dir];
         const FArrayBox& E_virt = m_electricField_virtual[dit].getFab();

         int additive = 0;
         int sign = 1-2*dir;
         if(anticyclic) sign = -sign;
         if(m_mesh.axisymmetric() && dir==1) {
            FArrayBox JaEy;
            JaEy.define(E_virt.box(),1);
            JaEy.copy(E_virt,0,0,1);
            JaEy.mult(Jnc[dit].getFab(),0,0,1);
            SpaceUtils::simpleStagGradComp( curlE_dir, grids[dit], 0,
                                            JaEy, 0,
                                            sign*dX[dirj], dirj,
                                            additive );
            curlE_dir.divide(Jfc[dit][dir],0,0,1);
         }
         else {
            SpaceUtils::simpleStagGradComp( curlE_dir, grids[dit], 0,
                                            E_virt, 0,
                                            sign*dX[dirj], dirj,
                                            additive );
	 }

      }

   }

#elif CH_SPACEDIM==1 // 2 comp NodeFArrayBox ==> 2 comp FArrayBox

   const string& geom_type = m_mesh.geomType();
   const LevelData<FArrayBox>& Xcc = m_mesh.getXcc();
   const LevelData<NodeFArrayBox>& Xnc = m_mesh.getXnc();

   for(DataIterator dit(grids); dit.ok(); ++dit) {

      int dirj = 0;
      int dirk = 1;

      FArrayBox& curlE_virt = m_curlE_virtual[dit];
      const FArrayBox& E_virt = m_electricField_virtual[dit].getFab();

      int additive = 0;
      if(m_mesh.axisymmetric() && geom_type=="sph_R") {
         FArrayBox rEz(E_virt.box(),1);
         rEz.copy(E_virt,dirk,0,1);
         rEz.mult(Xnc[dit].getFab(),0,0,1);
         SpaceUtils::simpleStagGradComp( curlE_virt, grids[dit], dirj,
                                         rEz, 0,
                                         -dX[0], 0,
                                         additive );
         curlE_virt.divide(Xcc[dit],0,dirj,1);
      }
      else {
         SpaceUtils::simpleStagGradComp( curlE_virt, grids[dit], dirj,
                                         E_virt, dirk,
                                         -dX[0], 0,
                                         additive );
      }

      if(m_mesh.axisymmetric()) {
         FArrayBox rEy(E_virt.box(),1);
         rEy.copy(E_virt,dirj,0,1);
         rEy.mult(Xnc[dit].getFab(),0,0,1);
         SpaceUtils::simpleStagGradComp( curlE_virt, grids[dit], dirk,
                                         rEy, 0,
                                         dX[0], 0,
                                         additive );
         curlE_virt.divide(Xcc[dit],0,dirk,1);
      }
      else {
         SpaceUtils::simpleStagGradComp( curlE_virt, grids[dit], dirk,
                                         E_virt, dirj,
                                         dX[0], 0,
                                         additive );
      }

   }

#endif

   SpaceUtils::exchangeFluxBox(m_curlE);
#if CH_SPACEDIM<3
   m_curlE_virtual.exchange();
#endif

}

void EMFields::setDivE()
{
   CH_TIME("EMFields::setDivE()");

   // get some Grid info
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const LevelData<EdgeDataBox>& Jec = m_mesh.getJec();
   const LevelData<NodeFArrayBox>& Jnc = m_mesh.getJnc();
   const RealVect& dX(m_mesh.getdX());

   LevelData<EdgeDataBox> EforDiv;
   EforDiv.define(grids,1,m_electricField.ghostVect());
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         EforDiv[dit][dir].copy(m_electricField[dit][dir]);
      }
   }
   m_field_bc->applyToEforDiv( EforDiv );

   // first do dir=0 to address singularity on axis for cyl/sph geoms
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& this_divE( m_divE[dit].getFab() );
      Box node_box = grids[dit];
      node_box.surroundingNodes();

      this_divE.setVal(0.0, node_box, 0, this_divE.nComp());
      int additive = 0;
      int dir = 0;
      const FArrayBox& this_E( EforDiv[dit][dir] );
      FArrayBox JaE0;
      JaE0.define(this_E.box(),1);
      JaE0.copy(this_E,0,0,1);
      JaE0.mult(Jec[dit][dir],0,0,1);
      SpaceUtils::simpleStagGradComp( this_divE, node_box, 0,
                                      JaE0, 0,
                                      dX[dir], dir,
                                      additive );
      this_divE.divide(Jnc[dit].getFab(),0,0,1);

   }
   m_field_bc->applyOnAxisDivBC( m_divE, m_electricField );

   // now do dir=1:SpaceDim
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& this_divE( m_divE[dit].getFab() );
      Box node_box = grids[dit];
      node_box.surroundingNodes();

      int additive = 1;
      for (int dir=1; dir<SpaceDim; dir++) {
         const FArrayBox& this_E( EforDiv[dit][dir] );
         SpaceUtils::simpleStagGradComp( this_divE, node_box, 0,
                                         this_E, 0,
                                         dX[dir], dir,
                                         additive );
      }

   }

}

void EMFields::setDivB()
{
   CH_TIME("EMFields::setDivB()");

   // get some Grid info
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const LevelData<FluxBox>& Jfc = m_mesh.getJfc();
   const LevelData<FArrayBox>& Jcc = m_mesh.getJcc();
   const RealVect& dX(m_mesh.getdX());

   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& this_divB( m_divB[dit] );
      const Box& cell_box = grids[dit];

      this_divB.setVal(0.0, cell_box, 0, this_divB.nComp());
      int additive = 0;
      for (int dir=0; dir<SpaceDim; dir++) {
         const FArrayBox& this_B( m_magneticField[dit][dir] );
         if(dir==0) {
            FArrayBox JaB0;
            JaB0.define(this_B.box(),1);
            JaB0.copy(this_B,0,0,1);
            JaB0.mult(Jfc[dit][dir],0,0,1);
            SpaceUtils::simpleStagGradComp( this_divB, cell_box, 0,
                                            JaB0, 0,
                                            dX[dir], dir,
                                            additive );
            this_divB.divide(Jcc[dit],0,0,1);
         }
         else {
            additive = 1;
            SpaceUtils::simpleStagGradComp( this_divB, cell_box, 0,
                                            this_B, 0,
                                            dX[dir], dir,
                                            additive );
         }
      }
   }

}

void EMFields::enforceGaussLaw( const LevelData<NodeFArrayBox>&  a_rho,
		                const Real                       a_time )
{
   CH_TIME("EMFields::enforceGaussLaw()");

   setChargeDensity( a_rho );
   solvePoisson( a_time );
   correctElectricField();

}

void EMFields::solvePoisson( const Real  a_time )
{
   CH_TIME("EMFields::solvePoisson()");

   setDivE();

   // 1) solve for correction potential: nabla^2(phi) = div(E) - rhoC/ep0
   // 2) compute the correction field:   E_corr = -nabla(phi)

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
   SpaceUtils::exchangeNodeFArrayBox(m_rhs0);

   // use copyTo to put rhs on container with dbl used for poisson solver
   LevelData<NodeFArrayBox>& rhs_vect0= *m_rhs_vector[0];
   m_rhs0.copyTo(rhs_vect0);

   // solve for phi
   int lbase = 0;
   m_poisson_solver.solve(m_phi_vector, m_rhs_vector, numlevels-1, lbase);

   // use copyTo to fill phi container with default dbl
   LevelData<NodeFArrayBox>& phi_vect0= *m_phi_vector[0];
   phi_vect0.copyTo(m_phi0);
#if CH_SPACEDIM==1
   // kluge to get fields correct at boundaries. Needed because Chombo's
   // interpretation of where boundaries are is different than used
   // here. The first ghost layer is not properly filled to maintain
   // nabla^2(phi) = -rho on boundaries. Also, the neumann BC is applied
   // at location of half cell inside the boundary.
   m_field_bc->applyPhiBC( m_phi0, m_chargeDensity, m_rhoCnorm_factor );
#endif

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

void EMFields::correctElectricField()
{
   CH_TIME("EMFields::correctElectricField()");
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

