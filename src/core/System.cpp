#include "System.H"
#include "CH_Timer.H"

#include "BoxIterator.H"
#include "DomainGrid.H"
#include "SpaceUtils.H"

#include "dataFileIO.H"
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "SpecialOperatorFactory.H"

#include "NamespaceHeader.H"


System::System( ParmParse&  a_pp )
   :
     m_iter_max(0),
     m_iter_max_particles(0),
     m_part_order_swap(false),
     m_freeze_particles_jacobian(true),
     m_quasi_freeze_particles_jacobian(false),
     m_use_mass_matrices(false),
     m_theta(0.5),
     m_advance_method(PIC_DSMC),
     m_mesh(nullptr),
     m_units(nullptr),
     m_dataFile(nullptr),
     m_meshInterp(nullptr),
     m_scattering(nullptr),
     m_verbosity(0),
     m_use_specialOps(false),
     m_rtol(1e-6),
     m_atol(1e-12),
     m_time_integrator(nullptr)
{
   ParmParse ppsys("system");
   parseParameters(ppsys);

   //m_units = new CodeUnits( ppsys );
   //if(!procID()) m_units->printParameters( procID() );
 
   createProblemDomain();          // create the problem domain
  
   DisjointBoxLayout grids;
   getDisjointBoxLayout( grids );  // define the disjointBoxLayout
   
   // initialize the coordinates and grid
   ParmParse ppgrid( "grid" );
   m_mesh = new DomainGrid( ppgrid, m_num_ghosts, m_domain, grids ); 
   
   m_units = new CodeUnits( ppsys, m_mesh->axisymmetric() );
   if(!procID()) m_units->printParameters( procID() );

   m_dataFile = new dataFileIO( ppsys, *m_mesh, *m_units );
     
   createMeshInterp();

   createState( a_pp );
   
   m_scattering = new ScatteringInterface( *m_units, m_pic_species_ptr_vect ); 

   createSpecialOperators();

   if(m_advance_method == PIC_DSMC) {
     m_time_integrator = new PICTimeIntegrator_DSMC;
   } else if(m_advance_method == PIC_EM_EXPLICIT) {
     m_time_integrator = new PICTimeIntegrator_EM_Explicit;
   } else if (m_advance_method == PIC_EM_SEMI_IMPLICIT) {
     m_time_integrator = new PICTimeIntegrator_EM_SemiImplicit;
   } else if(m_advance_method == PIC_EM_THETA_IMPLICIT) {
     m_time_integrator = new PICTimeIntegrator_EM_ThetaImplicit;
     dynamic_cast<PICTimeIntegrator_EM_ThetaImplicit*>(m_time_integrator)->theta(m_theta);
   }
   m_time_integrator->define( this,
                              m_pic_species_ptr_vect, 
                              m_emfields ); // m_emfields can't be NULL for implicit solvers
   m_time_integrator->setSolverParams( m_rtol, m_atol, m_iter_max );
}

System::~System()
{
   delete m_time_integrator;
   delete m_mesh;
   delete m_units;
   delete m_dataFile;
   if(m_meshInterp!=NULL) {
      delete m_meshInterp;
      m_meshInterp = NULL;
   }
   delete m_scattering;
}

void System::initialize( const int           a_cur_step,
                         const Real          a_cur_time,
                         const std::string&  a_restart_file_name )
{
   CH_TIME("System::initialize()");
   
   // read restart file
   if(!a_restart_file_name.empty()) {
      readCheckpointFile( a_restart_file_name );
   }

   // initialize the pic species
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      this_picSpecies->initialize(*m_units,a_cur_time,a_restart_file_name);
   }
   
   // initialize the scattering operators
    
   // set moments for initial mean free path calculation
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->scatter()) {
         this_picSpecies->binTheParticles();
         this_picSpecies->setNumberDensity();
         this_picSpecies->setEnergyDensity();
      }
   }
   m_scattering->initialize( m_pic_species_ptr_vect, *m_mesh, a_restart_file_name );
      
   // initialize the electromagnetic fields
   if (!m_emfields.isNull()) {
      m_emfields->initialize(a_cur_time,a_restart_file_name);
      setChargeDensity();
      setCurrentDensity();
   }
     
   // assert that m_emfields is not NULL if using forces 
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->forces()) {
         CH_assert(!m_emfields.isNull());
         break;
      }
   }

   m_time_integrator->initialize();
   
   // need to get time step from restart file to initialize MassMatrix 
   // for PC on restart (maybe better to do this in time integrator?) 
   if(m_use_mass_matrices && !a_restart_file_name.empty()) {
     HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );
     HDF5HeaderData header;
     header.readFromFile( handle );
     Real a_cur_dt   = header.m_real["cur_dt"];
     handle.close();
     setMassMatrices( a_cur_dt );
   }

}

void System::createProblemDomain()
{
   CH_TIME("System::createProblemDomain()");

   ParmParse ppgrid( "grid" );
   ppgrid.query( "num_ghosts", m_num_ghosts );
   int this_DIM = SpaceDim;

   // Set the grid size
   m_num_cells.resize( this_DIM );
   for (int i=0; i<this_DIM; ++i) m_num_cells[i] = 0;
   ppgrid.getarr( "num_cells", m_num_cells, 0, this_DIM );
   for (int i=0; i<this_DIM; ++i) CH_assert( m_num_cells[i]>0 );

   // Determine which spatial directions are periodic
   m_is_periodic.resize(this_DIM);
   vector<int> isPeriodic( this_DIM ); // why should I have to do this?
   ppgrid.getarr( "is_periodic", isPeriodic, 0, this_DIM );
   for (int dim=0; dim<SpaceDim; dim++)  {
      m_is_periodic[dim] = (isPeriodic[dim] == 1);
   }

   // Get the domain box decomposition parameters
   if (ppgrid.contains("config_decomp")) {
      m_config_decomp.resize( this_DIM );
      for (int i=0; i<this_DIM; ++i) m_config_decomp[i] = 0;
      ppgrid.getarr( "config_decomp", m_config_decomp, 0, this_DIM );
      for (int i=0; i<this_DIM; ++i) CH_assert( m_config_decomp[i]>0 );
   }

   int grid_verbosity;
   ppgrid.query( "verbosity", grid_verbosity );
   if (procID() == 0 && grid_verbosity) {
      cout << "====================== Spatial Grid Parameters =====================" << endl;
      cout << "  ghost layers = " << m_num_ghosts << endl;
      cout << "  number of cells = ";
      for (int i=0; i<SpaceDim; i++) cout << m_num_cells[i] << " ";
         cout << endl;
         cout << "  is periodic = ";
      for (int i=0; i<SpaceDim; i++) cout << m_is_periodic[i] << " ";
         cout << endl;
         if (m_config_decomp.size() > 0) {
            cout << "  configuration decomposition = ";
         for (int i=0; i<m_config_decomp.size(); i++)
            cout << m_config_decomp[i] << " ";
            cout << endl << endl;
         }
   }

   if(!procID()) cout << "Constructing ProblemDomain" << endl;

   IntVect hiEnd; 
   for (int dir=0; dir<SpaceDim; ++dir) hiEnd[dir] = m_num_cells[dir]-1;
   Box level0Domain(IntVect::Zero, hiEnd);

   bool isPeriodic_array[SpaceDim];
   for (int dir=0; dir<SpaceDim; ++dir) isPeriodic_array[dir] = (m_is_periodic[dir] == 1);

   m_domain.define( level0Domain.smallEnd(),
                    level0Domain.bigEnd(),
                    isPeriodic_array );

   if(!procID()) cout << "Done constructing ProblemDomain" << endl << endl;
  
}

void System::getDisjointBoxLayout( DisjointBoxLayout&  a_grids )
{
   CH_TIME("System::getDisjointBoxLayout()");

   if(!procID()) cout << "Constructing DisjointBoxLayout" << endl;

   Vector<Box> boxes;
   const Box& domain_box = m_domain.domainBox();
   
   // some AMR stuff and the mpi stuff for one of the particle handling methods
   // requires using boxes of a fixed length in each direction. Ensure that is the case
   IntVect boxSize;
   boxSize[0] = domain_box.size(0)/m_config_decomp[0];
   for (int dir=1; dir<SpaceDim; ++dir) {
      boxSize[dir] = domain_box.size(dir)/m_config_decomp[dir];
      CH_assert(boxSize[dir]==boxSize[0]);
   }
   
   // Chop up the configuration space domain box over the number of processors specified
   // for this block (single block here).  At this point, we insist that the box 
   // decomposes uniformly, or an error is thrown.
   IntVect n_loc = IntVect::Zero;
   for (int dir=0; dir<SpaceDim; ++dir) {
      int decomp_dir = m_config_decomp[dir];
      if (domain_box.size(dir)%decomp_dir != 0) {
	//stringstream msg("Decomposition in configuration direction ", ios_base::out|ios_base::ate);
        //msg << dir << " does not evenly divide domain dimension";
        //MayDay::Error( msg.str().c_str() );
      }
      else {
	n_loc[dir] = domain_box.size(dir) / decomp_dir;
      }
   }

   int box_cell_num(1);
   for (int dir=0; dir<SpaceDim; ++dir) {
      box_cell_num *= n_loc[dir];
   }
    
   if (box_cell_num > 0) {
      IntVect box_size(n_loc);
      Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size-1);
      IntVect lo = IntVect::Zero;
      IntVect hi;
      for (int dir=0; dir<SpaceDim; ++dir) {
         hi[dir] = domain_box.size(dir)/n_loc[dir]-1;
      }
      Box skeleton(lo, hi);
      BoxIterator bit(skeleton);
      for (bit.begin();bit.ok();++bit) {
         Box thisBox = patch + bit()*box_size;
         boxes.push_back(thisBox);
      }

   }
   else {
      MayDay::Error( "Configuration domain box cannot be load balanced" );
   }
  
   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
   LoadBalance( procMap, boxes );

   a_grids.define( boxes, procMap, m_domain );
   a_grids.close();
   
   if(!procID()) cout << "Done constructing DisjointBoxLayout" << endl << endl;

   if(!procID() && m_verbosity) {
      for (int n=0; n<boxes.size(); n++) {
         const Box& local_box = boxes[n];
         cout << " box " << local_box << " is assigned to process " << procMap[n] << endl;
      }
      cout << endl;
   }

}


void System::createMeshInterp()
{
   CH_TIME("System::createMeshInterp()");

   // get some mesh information
   const ProblemDomain& domain(m_mesh->getDomain()); 
   const RealVect& meshSpacing(m_mesh->getdX());
   const RealVect& meshOrigin(m_mesh->getXmin());
  
   if(m_meshInterp!=NULL) {
      delete m_meshInterp;
   }
   
   m_meshInterp = static_cast<MeshInterp*> (new MeshInterp( domain.domainBox(),
                                                            meshSpacing,
                                                            meshOrigin ));

}


void System::createState( ParmParse&  a_pp )
{

   createPICspecies();

   createEMfields();

}

void System::createEMfields()
{
   
   if(!procID()) {
      cout << "Creating Electromagnetic fields object..." << endl << endl;
   }
   
   ParmParse ppflds( "em_fields" );

   EMVecType em_vec_type;
   if (m_advance_method == PIC_EM_SEMI_IMPLICIT) {
     em_vec_type = e_only;
   } else {
     em_vec_type = e_and_b;
   }
   
   bool use_fields = false;
   ppflds.query("use",use_fields);
   if(use_fields) {
      bool verbose = true;

      m_emfields = RefCountedPtr<ElectroMagneticFields>
                                  (new ElectroMagneticFields( ppflds,
                                                              *m_mesh,
                                                              *m_units,
                                                              verbose,
                                                              em_vec_type)  ); 
      
      // number of components for mass matrices depends on particle weighting scheme.
      // I don't like having this here.... Also need to check that all charged 
      // particles are using the same weighting scheme ...
      bool init_mass_matrices = false;
      if (m_use_mass_matrices) init_mass_matrices = true;
#ifdef MASS_MATRIX_TEST
      CH_assert(!m_use_mass_matrices);
      init_mass_matrices = true;
#endif
      if (init_mass_matrices) {
         InterpType sp_interp = UNKNOWN;
         for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
            PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
            if(this_picSpecies->charge() == 0) continue;
            InterpType this_interp = this_picSpecies->getInterpType();
            if(sp_interp==UNKNOWN) sp_interp = this_interp;
            else CH_assert(sp_interp==this_interp);
         }
         m_emfields->initializeMassMatrices(sp_interp, m_use_mass_matrices);
      }

      if(m_advance_method == PIC_DSMC ) {
         if(!procID()) cout << "advance_method PIC_DSMC cannot be used with electromagnetic fields on" << endl;
         exit(EXIT_FAILURE);
      }
      if(!procID()) cout << "Finished creating Electromagnetic fields object" << endl << endl;
   }
   else {
      m_emfields = RefCountedPtr<ElectroMagneticFields>(NULL); // don't have to do this, but for clarity 
      if(!procID()) cout << "Electromagnetic fields are not being used" << endl << endl;
   }

}


void System::createPICspecies()
{
   // Create all PIC species
   
   if(!procID()) {
      cout << "Creating PIC species..." << endl << endl;
   }

   bool more_vars(true);
   int species;
   while(more_vars) { // look for pic species...
 
      species = m_pic_species_ptr_vect.size();

      stringstream s;
      s << "pic_species." << species; 

      ParmParse ppspc( s.str().c_str() );
     
      string name;
      if(ppspc.contains("name")) {
         ppspc.get("name",name);
      } 
      else {
         more_vars = false;
      }
   
      if(more_vars) {
      
         // Create the pic species object
         PicSpecies* picSpecies = new PicSpecies( ppspc, species, name, *m_meshInterp, *m_mesh );

         // Add the new pic species object to the vector of pointers to PicSpecies
         m_pic_species_ptr_vect.push_back(PicSpeciesPtr(picSpecies));

      }

   }

   if(!procID()) {
      cout << "Finished creating " << m_pic_species_ptr_vect.size() << " PIC species" << endl << endl;
   }

}

void System::createSpecialOperators()
{
   // Create the vector of special operator (pointers)
   
   SpecialOperatorFactory  specialOpFactory;
   
   if(!procID()) {
      cout << "Adding special operators..." << endl;
   }

   bool more_ops(true);
   string name0;
   int special_op_num = -1;
   while(more_ops) { // look for special operator...
 
      special_op_num = special_op_num + 1;

      stringstream s;
      s << "special_operator." << special_op_num; 
      ParmParse ppspop( s.str().c_str() );
     
      if(ppspop.contains("model")) {
         m_use_specialOps = true;
         m_specialOps = specialOpFactory.create( ppspop, *m_mesh, *m_units, 1 );
         m_specialOps->updateOp(0.0); // should make an initializeOp() function
      } 
      else {
         more_ops = false;
      }
   
   }

   if(!procID()) {
      cout << "Done adding special operators" << endl << endl;
   }


}


void System::writePlotFile( const int     a_cur_step,
                            const double  a_cur_time,
                            const double  a_cur_dt )
{
   CH_TIME("System::writePlotFile()");
   
   if (!m_emfields.isNull()) {
      if(m_emfields->writeRho()) setChargeDensity();
      if(m_emfields->writeExB()) {
         m_emfields->setPoyntingFlux();
      }
      if(m_emfields->writeDivs()) {
         m_emfields->setDivE();
         m_emfields->setDivB();
      }
      m_dataFile->writeElectroMagneticFieldsDataFile( *m_emfields,
                                                      a_cur_step, a_cur_time ); 
   }

   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {

      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
  
      this_picSpecies->setNumberDensity();
      this_picSpecies->setMomentumDensity();
      this_picSpecies->setEnergyDensity();

      if(this_picSpecies->charge() == 0) {
         m_dataFile->writeNeutralSpeciesDataFile( *this_picSpecies, 
                                                  s, a_cur_step, a_cur_time );
      }
      else {
         if(m_writeSpeciesChargeDensity) {
            this_picSpecies->setChargeDensity();      
            this_picSpecies->setChargeDensityOnFaces(); 
            this_picSpecies->setChargeDensityOnNodes(); 
         }     
         if(m_writeSpeciesCurrentDensity) this_picSpecies->setCurrentDensity();      
         m_dataFile->writeChargedSpeciesDataFile( *this_picSpecies, 
                                                  s, a_cur_step, a_cur_time,
                                                  m_writeSpeciesChargeDensity,
                                                  m_writeSpeciesCurrentDensity );
      }

   }
   
}

void System::writeHistFile( const int   a_cur_step,
                            const Real  a_cur_time,
                            const Real  a_cur_dt, 
                            const bool  a_startup )
{
   CH_TIME("System::writeHistFile()");
   
   if(a_startup) setupHistFile(a_cur_step);
   
   // compute solver probes
   int l_exit_status=0, nl_exit_status=0;
   int l_total_iter=0, nl_total_iter=0;
   int l_last_iter=0, nl_iter=0;
   Real nl_abs_res=0.0, nl_rel_res=0.0;
   if(m_solver_probes && !a_startup) {
      m_time_integrator->getConvergenceParams( l_exit_status, l_last_iter, l_total_iter,
                                              nl_exit_status, nl_iter, nl_total_iter,
                                              nl_abs_res, nl_rel_res );
   }
   
   // compute the field probes
   Real energyE_joules, energyB_joules;
   Real max_wc0dt = 0.0;
   if(m_field_probes) {
      energyE_joules = m_emfields->electricFieldEnergy();
      energyB_joules = m_emfields->magneticFieldEnergy();
      max_wc0dt = m_emfields->max_wc0dt(*m_units,a_cur_dt);
   }
   
   // compute the scattering probes
   Real energyIzn_joules, energyExc_joules;
   if(m_scattering_probes) {
      m_scattering->setScatteringProbes();
      energyIzn_joules = m_scattering->ionizationEnergy();
      energyExc_joules = m_scattering->excitationEnergy();
   }
   
   // compute the species probes
   const int numSpecies = m_pic_species_ptr_vect.size();
   std::vector<int> num_particles(numSpecies);
   std::vector<vector<Real>> global_moments(numSpecies);
   std::vector<Real> max_wpdt(numSpecies);
   std::vector<Real> max_wcdt(numSpecies);
   for( int species=0; species<numSpecies; species++) {
      if(!m_species_probes[species]) continue;
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[species]);
      num_particles[species] = this_picSpecies->numParticles();
      this_picSpecies->globalMoments(global_moments.at(species));
      max_wpdt.at(species) = this_picSpecies->max_wpdt(*m_units,a_cur_dt);
      max_wcdt.at(species) = max_wc0dt*abs(this_picSpecies->charge())/this_picSpecies->mass();
   }
   
   //
   // write probe data to the history file 
   //

   //FILE *histFile;
   std::ofstream histFile;
   if(!procID()) {
      //if((histFile = fopen("history.txt", "a")) != NULL) {
      histFile.open("history.txt", std::ios_base::app);
      if(histFile.is_open()) {
         histFile << a_cur_step << " ";
         histFile << std::setprecision(m_history_precision) << std::scientific << a_cur_time << " ";
         if(m_solver_probes) {
            histFile << nl_exit_status << " ";
            histFile << nl_iter << " ";
            histFile << nl_total_iter << " ";
            histFile << nl_abs_res << " ";
            histFile << nl_rel_res << " ";
            histFile << l_exit_status << " ";
            histFile << l_last_iter << " ";
            histFile << l_total_iter << " ";
         }
         if(m_field_probes) {
            histFile << energyE_joules << " ";
            histFile << energyB_joules << " ";
         }
         if(m_field_bdry_probes) {
            const RealVect& intSdA_lo = m_emfields->getIntSdA_lo();  
            const RealVect& intSdA_hi = m_emfields->getIntSdA_hi();  
            const RealVect& intSdAdt_lo = m_emfields->getIntSdAdt_lo();  
            const RealVect& intSdAdt_hi = m_emfields->getIntSdAdt_hi(); 
            for(int dir=0; dir<SpaceDim; ++dir) { 
               if(m_is_periodic[dir]) continue;
               histFile << intSdA_lo[dir] << " ";
               histFile << intSdA_hi[dir] << " ";
               histFile << intSdAdt_lo[dir] << " ";
               histFile << intSdAdt_hi[dir] << " ";
            }
         }
         if(m_scattering_probes) {
            histFile << energyIzn_joules << " ";
            histFile << energyExc_joules << " ";
         }
         for( int species=0; species<numSpecies; species++) {
            if(!m_species_probes[species]) continue;
            histFile << num_particles[species] << " ";
            std::vector<Real>& species_moments = global_moments.at(species);
            for(int mom=0; mom<species_moments.size(); ++mom) {
               histFile << species_moments.at(mom) << " ";
            }
            histFile << max_wpdt.at(species) << " ";
            histFile << max_wcdt.at(species) << " ";
         }
         histFile << "\n";
         histFile.close();
      }
   }  

}

void System::setupHistFile(const int a_cur_step) 
{
   CH_TIME("System::setupHistFile()");
   if(!procID()) cout << "setting up the history file at step " << a_cur_step << endl;

   std::vector<std::string> probe_names;
   probe_names.push_back("#0 step number");
   probe_names.push_back("#1 time [code units]");
   //probe_names.push_back("#2 time step [code units]");

   // parse input for different probes to include in the history file
   ParmParse pphist("history");
   
   m_history_precision = 5;
   pphist.query("precision", m_history_precision);
   
   m_solver_probes = false;
   if(m_time_integrator!=NULL) {
      pphist.query("solver_probes", m_solver_probes);
      if(!procID()) cout << "solver_probes = " << m_solver_probes << endl;
   }
   
   m_field_probes = false;
   if(!m_emfields.isNull()) {
      pphist.query("field_probes", m_field_probes);
      if(!procID()) cout << "field_probes = " << m_field_probes << endl;
   }
   
   m_field_bdry_probes = false;
   if(!m_emfields.isNull()) {
      pphist.query("field_bdry_probes", m_field_bdry_probes);
      if(!procID()) cout << "field_bdry_probes = " << m_field_bdry_probes << endl;
   }
   
   m_scattering_probes = false;
   pphist.query("scattering_probes", m_scattering_probes);
   if(!procID()) cout << "scattering_probes = " << m_scattering_probes << endl;
 
   const int numSpecies = m_pic_species_ptr_vect.size();
   if(numSpecies>0) { 
      m_species_probes.resize(numSpecies,false);
      for (int species=0; species<numSpecies; species++) {
         stringstream ss;
         ss << "species" << species << "_probes";
         if(pphist.contains(ss.str().c_str())) {
            bool this_boolean; 
            pphist.get(ss.str().c_str(), this_boolean);
            m_species_probes[species] = this_boolean;
         }
         if(!procID()) cout << ss.str() << " = " << m_species_probes[species] << endl;
      }
   }
   
   if(m_solver_probes) {
      stringstream ss;
      ss << "#" << probe_names.size() << " nonlinear solver exit status";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " nonlinear solver iterations";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " nonlinear solver total iterations";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " nonlinear solver abs residual";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " nonlinear solver rel residual";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " linear solver exit status";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " linear solver last iteration";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " linear solver total iterations";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
   }
   
   if(m_field_probes) {
      stringstream ss;
      ss << "#" << probe_names.size() << " electric field energy [Joules]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " magnetic field energy [Joules]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
   }
   
   if(m_field_bdry_probes) {
      for (int dir=0; dir<SpaceDim; dir++) {
         if(m_is_periodic[dir]) continue;
         stringstream ss;
         ss << "#" << probe_names.size() << " int S=ExB/mu0 dA from lo-side, dir = " << dir << " [Joules/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " int S=ExB/mu0 dA from hi-side, dir = " << dir << " [Joules/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " int int S dA dt from lo-side, dir = " << dir << " [Joules]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " int int S dA dt from hi-side, dir = " << dir << " [Joules]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
      }
   }
   
   if(m_scattering_probes) {
      stringstream ss;
      ss << "#" << probe_names.size() << " ionization energy [Joules]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " excitation energy [Joules]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
   }
   
   for( int species=0; species<numSpecies; species++) {
      if(!m_species_probes[species]) continue;
      stringstream ss;
      ss << "#" << probe_names.size() << " species " << species << ": macro particles";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " species " << species << ": mass [kg]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " species " << species << ": X-momentum [kg-m/s]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " species " << species << ": Y-momentum [kg-m/s]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " species " << species << ": Z-momentum [kg-m/s]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " species " << species << ": energy [Joules]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " species " << species << ": max wp*dt";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " species " << species << ": max wc*dt";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
   }
   
   // setup or restart the history file
   if(!procID()) {

      const std::string hf_name = "history.txt"; 
      if(a_cur_step==0) {

         std::ofstream histFile(hf_name,std::ofstream::out);
         if(histFile.is_open()) {
            for (auto n=0; n<probe_names.size(); n++) {
               histFile << probe_names.at(n) << " \n";
            }
            histFile.close();
         }
         else {
            cout << "Problem creating/opening history.txt file" << endl;
         }

      } 
      else {
         std::ifstream histFile;
         std::ofstream histFile_temp;
         std::string this_line;
         int this_step;
      
         histFile.open(hf_name,std::ifstream::in);
         histFile_temp.open("history_temp.txt",std::ofstream::out);
         if(histFile.is_open()) {
            while(getline(histFile, this_line)) {
               stringstream ss;
               ss << this_line;
               string temp;
               ss >> temp; // first element is step number
               if( stringstream(temp) >> this_step && this_step >= a_cur_step ) {
                  break;
               }
               else {
                  histFile_temp << this_line; 
                  histFile_temp << "\n"; 
               }
            }
            histFile.close();
            histFile_temp.close();
            remove(hf_name.c_str());
            rename("history_temp.txt",hf_name.c_str());
         }
         else {
            cout << "Problem opening history.txt file on restart" << endl;
         }
      }

      cout << "finished setting up the history file " << endl << endl;
   }

}

void System::writeCheckpointFile( HDF5Handle&  a_handle,
                            const int          a_cur_step,
                            const double       a_cur_time,
                            const double       a_cur_dt )
{
   CH_TIME("System::writeCheckpointFile()");

   MPI_Barrier(MPI_COMM_WORLD);
   if(!procID()) cout << "Writing checkpoint file at step " << a_cur_step << endl;

   HDF5HeaderData header;
   header.m_int["cur_step"]  = a_cur_step;
   header.m_real["cur_time"] = a_cur_time;
   header.m_real["cur_dt"]   = a_cur_dt;
   
   header.writeToFile( a_handle );
   
   if(m_scattering->numScatter()>0) m_scattering->writeCheckpoint( a_handle );
   
   const int write_old_data = m_time_integrator->prepForCheckpoint();
   
   // write the field data
   if (!m_emfields.isNull()) {
      m_dataFile->writeCheckpointEMFields( a_handle, write_old_data, *m_emfields );
   }

   // write the particle data
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      m_dataFile->writeCheckpointParticles( a_handle, *this_picSpecies, s );
   }

}

void System::readCheckpointFile( const std::string&  a_chkpt_fname )
{
   CH_TIME("System::readCheckpointFile()");

   if(!procID()) cout << "Reading restart file " << a_chkpt_fname << endl;

#ifdef CH_USE_HDF5
   HDF5Handle handle( a_chkpt_fname, HDF5Handle::OPEN_RDONLY );

   HDF5HeaderData header;
   header.readFromFile( handle );

   handle.close();

#else
   MayDay::Error("restart only defined with hdf5");
#endif
   
}

void System::parseParameters( ParmParse&  a_ppsys )
{
   a_ppsys.query("advance_method", m_advance_method);
   if(m_advance_method == PICMC_EXPLICIT) m_advance_method = PIC_EM_EXPLICIT;
   if(m_advance_method == PICMC_SEMI_IMPLICIT) m_advance_method = PIC_EM_SEMI_IMPLICIT;
   if(m_advance_method == PICMC_FULLY_IMPLICIT) m_advance_method = PIC_EM_THETA_IMPLICIT;
   if(m_advance_method == PIC_EM_SEMI_IMPLICIT) {
      a_ppsys.get("iter_max",m_iter_max);
      a_ppsys.query("iter_max_particles",m_iter_max_particles);
      a_ppsys.query("part_order_swap",m_part_order_swap);
      a_ppsys.query("abs_tol", m_atol);
      a_ppsys.query("rel_tol", m_rtol);
   } else if(m_advance_method == PIC_EM_THETA_IMPLICIT) {
      a_ppsys.get("iter_max",m_iter_max);
      a_ppsys.query("iter_max_particles",m_iter_max_particles);
      a_ppsys.query("part_order_swap",m_part_order_swap);
      a_ppsys.query("freeze_particles_jacobian",m_freeze_particles_jacobian);
      a_ppsys.query("quasi_freeze_particles_jacobian",m_quasi_freeze_particles_jacobian);
      if(m_quasi_freeze_particles_jacobian) m_freeze_particles_jacobian=false;
      a_ppsys.query("use_mass_matrices",m_use_mass_matrices);
      if(m_use_mass_matrices) {
         m_freeze_particles_jacobian = false;
         m_part_order_swap = true; // for now, so smoke tests don't fail...
      }
      a_ppsys.query("abs_tol", m_atol);
      a_ppsys.query("rel_tol", m_rtol);
      a_ppsys.query("theta_parameter",m_theta);
      CH_assert(m_theta>=0.5 && m_theta<=1.0);
   }

   a_ppsys.query("write_species_charge_density", m_writeSpeciesChargeDensity);
   a_ppsys.query("write_species_current_density", m_writeSpeciesCurrentDensity);
 
   if(!procID()) {
      cout << "advance method  = " << m_advance_method << endl;
      if (m_advance_method == PIC_EM_SEMI_IMPLICIT) {
         cout << "  iter_max = " << m_iter_max << endl;
         cout << "  iter_max_particles = " << m_iter_max_particles << endl;
         cout << "  part_order_swap = " << m_part_order_swap << endl;
         cout << "  abs_tol  = " << m_atol << endl;
         cout << "  rel_tol  = " << m_rtol << endl;
      } else if(m_advance_method == PIC_EM_THETA_IMPLICIT) {
         cout << "  iter_max = " << m_iter_max << endl;
         cout << "  iter_max_particles = " << m_iter_max_particles << endl;
         cout << "  part_order_swap = " << m_part_order_swap << endl;
         cout << "  freeze_particles_jacobian = " << m_freeze_particles_jacobian << endl;
         cout << "  quasi_freeze_particles_jacobian = " << m_quasi_freeze_particles_jacobian << endl;
         cout << "  use_mass_matrices = " << m_use_mass_matrices<< endl;
         cout << "  abs_tol  = " << m_atol << endl;
         cout << "  rel_tol  = " << m_rtol << endl;
         cout << "  theta_parameter = " << m_theta << endl; 
      }
      cout << "write species charge density  = " << m_writeSpeciesChargeDensity << endl;
      cout << "write species current density = " << m_writeSpeciesCurrentDensity << endl << endl;
   }
 
}

void System::preTimeStep( const Real&  a_cur_time,
                          const Real&  a_dt,
                          const int&   a_step_number )
{  
   CH_TIME("System::preTimeStep()");
      
   m_time_integrator->preTimeStep( a_cur_time, a_dt, a_step_number );

   // update old particle values and create inflow particles  
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      this_picSpecies->createInflowParticles( a_cur_time, a_dt );
      this_picSpecies->updateOldParticlePositions();
      this_picSpecies->updateOldParticleVelocities();
   }

}

void System::advance( const Real&  a_cur_time,
                      const Real&  a_dt,
                      const int&   a_step_number )
{  
   CH_TIME("System::advance()");
   m_time_integrator->timeStep( a_cur_time, a_dt, a_step_number );

}

void System::setChargeDensity()
{
   CH_TIME("System::setChargeDensity()");
   const DisjointBoxLayout& grids(m_mesh->getDBL());
 
   LevelData<NodeFArrayBox>& chargeDensity = m_emfields->getChargeDensity();
   
   // set the charge density to zero
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& this_rho_nodes = chargeDensity[dit].getFab();
      this_rho_nodes.setVal(0.0);
   }

   // loop over all pic species and add the charge density 
   for (int s=0; s<m_pic_species_ptr_vect.size(); ++s) {
      const PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      int this_charge = this_picSpecies->charge();
      if(this_charge != 0) {

         this_picSpecies->setChargeDensityOnNodes();
         const LevelData<NodeFArrayBox>& species_rho = this_picSpecies->getChargeDensityOnNodes();
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            const FArrayBox& this_species_rho = species_rho[dit].getFab();
            FArrayBox& this_rho = chargeDensity[dit].getFab();
            this_rho.plus(this_species_rho,0,0,chargeDensity.nComp());
         }

      }
   }

}

void System::setCurrentDensity()
{
   CH_TIME("System::setCurrentDensity()");
   const DisjointBoxLayout& grids(m_mesh->getDBL());
 
   LevelData<EdgeDataBox>&  currentDensity = m_emfields->getCurrentDensity();
   LevelData<NodeFArrayBox>&  currentDensity_virtual = m_emfields->getVirtualCurrentDensity();
   
   // set the current density member to zero
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         currentDensity[dit][dir].setVal(0.0);
      }
      if(SpaceDim<3) {
         FArrayBox& this_Jv_nodes = currentDensity_virtual[dit].getFab();
         this_Jv_nodes.setVal(0.0);
      }
   }

   // loop over all pic species and add the current density 
   for (int s=0; s<m_pic_species_ptr_vect.size(); ++s) {
      const PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      int this_charge = this_picSpecies->charge();
      if(this_charge != 0) {

         this_picSpecies->setCurrentDensity();
         const LevelData<EdgeDataBox>& species_J = this_picSpecies->getCurrentDensity();
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir=0; dir<SpaceDim; dir++) {
               const FArrayBox& this_species_Jdir = species_J[dit][dir];
               currentDensity[dit][dir].plus(this_species_Jdir,0,0,1);
            }
         }
         if(SpaceDim<3) {
            const LevelData<NodeFArrayBox>& species_Jv = this_picSpecies->getCurrentDensity_virtual();
            for(DataIterator dit(grids); dit.ok(); ++dit) {
               const FArrayBox& this_species_Jv = species_Jv[dit].getFab();
               FArrayBox& this_Jv = currentDensity_virtual[dit].getFab();
               this_Jv.plus(this_species_Jv,0,0,currentDensity_virtual.nComp());
               //if(!procID()) SpaceUtils::inspectFArrayBox( this_Jv, this_Jv.box(), 1);
            }
         }

      }
   }

   // call exchange 
   // JRA 1-13-22
   // When I call exchange for virtual J here, the numerical energy regression test 
   // that uses jfnk fails. However, when using jfnk it is possible that a machine-level
   // difference can grow into a tolerance-level difference
   //SpaceUtils::exchangeEdgeDataBox(currentDensity);
   //if(SpaceDim<3) SpaceUtils::exchangeNodeFArrayBox(currentDensity_virtual);

}

void System::setMassMatrices( const Real  a_dt )
{
   CH_TIME("System::setMassMatrices()");
   CH_assert(SpaceDim==1); // only works for 1D so for
   
   LevelData<EdgeDataBox>&  Jtilde = m_emfields->getJtilde();
   LevelData<NodeFArrayBox>&  Jtildev = m_emfields->getVirtualJtilde();
   //
   LevelData<EdgeDataBox>&  sigma_xx = m_emfields->getSigmaxx();
   LevelData<EdgeDataBox>&  sigma_xy = m_emfields->getSigmaxy();
   LevelData<EdgeDataBox>&  sigma_xz = m_emfields->getSigmaxz();
   LevelData<NodeFArrayBox>&  sigma_yx = m_emfields->getSigmayx();
   LevelData<NodeFArrayBox>&  sigma_yy = m_emfields->getSigmayy();
   LevelData<NodeFArrayBox>&  sigma_yz = m_emfields->getSigmayz();
   LevelData<NodeFArrayBox>&  sigma_zx = m_emfields->getSigmazx();
   LevelData<NodeFArrayBox>&  sigma_zy = m_emfields->getSigmazy();
   LevelData<NodeFArrayBox>&  sigma_zz = m_emfields->getSigmazz();
   
   // set all to zero
   SpaceUtils::zero( Jtilde );
   SpaceUtils::zero( Jtildev );
   //
   SpaceUtils::zero( sigma_xx );
   SpaceUtils::zero( sigma_xy );
   SpaceUtils::zero( sigma_xz );
   SpaceUtils::zero( sigma_yx );
   SpaceUtils::zero( sigma_yy );
   SpaceUtils::zero( sigma_yz );
   SpaceUtils::zero( sigma_zx );
   SpaceUtils::zero( sigma_zy );
   SpaceUtils::zero( sigma_zz );

   // loop over all pic species and add their contribution 
   for (int s=0; s<m_pic_species_ptr_vect.size(); ++s) {
      const PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      this_picSpecies->accumulateMassMatrices( sigma_xx, sigma_xy, sigma_xz,
                                               sigma_yx, sigma_yy, sigma_yz,
                                               sigma_zx, sigma_zy, sigma_zz, 
                                               Jtilde, Jtildev, *m_emfields, a_dt );
   }
   
   // add ghost cells to valid cells
   LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   Jtilde.exchange( Jtilde.interval(), m_mesh->reverseCopier(), addEdgeOp );
   sigma_xx.exchange( sigma_xx.interval(), m_mesh->reverseCopier(), addEdgeOp );
   sigma_xy.exchange( sigma_xy.interval(), m_mesh->reverseCopier(), addEdgeOp );
   sigma_xz.exchange( sigma_xz.interval(), m_mesh->reverseCopier(), addEdgeOp );

   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   Jtildev.exchange( Jtildev.interval(), m_mesh->reverseCopier(), addNodeOp );
   sigma_yx.exchange( sigma_yx.interval(), m_mesh->reverseCopier(), addNodeOp );
   sigma_yy.exchange( sigma_yy.interval(), m_mesh->reverseCopier(), addNodeOp );
   sigma_yz.exchange( sigma_yz.interval(), m_mesh->reverseCopier(), addNodeOp );
   sigma_zx.exchange( sigma_zx.interval(), m_mesh->reverseCopier(), addNodeOp );
   sigma_zy.exchange( sigma_zy.interval(), m_mesh->reverseCopier(), addNodeOp );
   sigma_zz.exchange( sigma_zz.interval(), m_mesh->reverseCopier(), addNodeOp );
   
}

void System::scatterParticles( const Real&  a_dt )
{  
   if (m_scattering->numScatter()==0) return;

   CH_TIME("System::scatterParticles()");

   const Real tscale = m_units->getScale(m_units->TIME);
   const Real dt_sec = a_dt*tscale;      

   // prepare all species to be scattered for scattering
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->scatter()) {
         this_picSpecies->binTheParticles();
         this_picSpecies->setNumberDensity();
      }
   }

   m_scattering->applyScattering( m_pic_species_ptr_vect, *m_mesh, dt_sec );
   
}   

void System::postTimeStep( Real&        a_cur_time,
                           const Real&  a_dt,
                           int&         a_step_number )
{  
   CH_TIME("System::postTimeStep()");
   
   if(!m_emfields.isNull()) {
      m_emfields->postTimeStep( a_dt );
   }

   m_time_integrator->postTimeStep( a_cur_time, a_dt );

   // apply special operators
   if(m_specialOps) {
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) { // loop over pic species
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         m_specialOps->applyOp(*this_picSpecies,a_dt);
      }
      m_specialOps->updateOp(a_dt);
   }
   
   // remove particles that have left the physical domain and create new ones
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      this_picSpecies->removeOutflowParticles();
   }
   
   a_cur_time = a_cur_time + a_dt;
   a_step_number = a_step_number + 1;
}

Real System::fieldsDt( const int a_step_number )
{
   Real fldsDt = DBL_MAX;
   if (!m_emfields.isNull() && m_emfields->advance()) {
      fldsDt = m_emfields->stableDt();
      if(!procID()) cout << "electromagnetic fields time step = " << fldsDt << endl;
   }
   return fldsDt;
}


Real System::partsDt( const int a_step_number )
{

   Real stableDt = DBL_MAX;
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      this_picSpecies->setStableDt(*m_units);
      stableDt = min(stableDt,this_picSpecies->stableDt());
   }
  
   return stableDt;

}


Real System::scatterDt( const int a_step_number )
{
   CH_TIME("System::scatterDt()");

   Real scatterDt = DBL_MAX;
   
   if(m_scattering->numScatter()>0) {

      // compute the density and energy moments for mean free time calculation
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->scatter()) {
            this_picSpecies->setNumberDensity();
            this_picSpecies->setEnergyDensity();
         }
      }
      
      scatterDt = m_scattering->scatterDt( m_pic_species_ptr_vect );
 
      // convert to code units
      const Real tscale = m_units->getScale(m_units->TIME);
      scatterDt = scatterDt/tscale;
   
   }
   
   return scatterDt;

}


Real System::specialOpsDt( const int a_step_number )
{
   Real specialOpsDt = DBL_MAX;
   if(m_use_specialOps) specialOpsDt = m_specialOps->specialOpsDt();
   return specialOpsDt;
}

void System::adaptDt( bool&  a_adapt_dt )
{
   if(m_advance_method==PIC_EM_EXPLICIT) {
      if(m_emfields->advance()) a_adapt_dt = false;
   }
   if(m_advance_method==PIC_EM_SEMI_IMPLICIT) {
      if(m_emfields->advance()) a_adapt_dt = false;
   }
}

void System::preRHSOp( const ODEVector<System>&  a_U,
                       const Real                a_dt,
                       const bool                a_from_emjacobian )
{  
  CH_TIME("System::preRHSOp()");
  
  CH_assert(!m_emfields.isNull());
    
  // half dt advance of particle positions and velocities
  // and, if advancing E, compute current density
  //
  // xbar = xn + dt/2*vbar
  // vbar = vn + dt/2*q/m*(E(xbar) + vbar x B(xbar))
  // Jbar = Sum_sSum_p(qp*S(xbar-xg)*vbar)/dV
 
  if(a_from_emjacobian) { // called from linear stage of jfnk

     if (m_use_mass_matrices) {

        m_emfields->computeJfromMassMatrices();

     }
     else {

        CH_assert(m_iter_max_particles>0);
        for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
           auto this_picSpecies(m_pic_species_ptr_vect[s]);
           if (m_quasi_freeze_particles_jacobian) { // freeze particle positions during linear stage
              this_picSpecies->interpolateFieldsToParticles( *m_emfields );
              this_picSpecies->addExternalFieldsToParticles( *m_emfields ); 
              this_picSpecies->advanceVelocities( a_dt, true );
           }
           else {
              this_picSpecies->repositionOutcastsAndApplyForces( *m_emfields, 
                                                                 a_dt, true );
              this_picSpecies->advanceParticlesIteratively( *m_emfields, a_dt, 
                                                            m_part_order_swap, 
                                                            m_iter_max_particles );
           }
        }
        if (m_emfields->advanceE()) setCurrentDensity();
#ifdef MASS_MATRIX_TEST
        setMassMatrices( a_dt );
        m_emfields->computeJfromMassMatrices();
#endif
     }

  }
  else { // called from nonlinear stage of jfnk or from picard solver
  
     for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
        auto this_picSpecies(m_pic_species_ptr_vect[s]);
        this_picSpecies->repositionOutcastsAndApplyForces( *m_emfields, 
                                                           a_dt, true );
        if(m_iter_max_particles>0) {
           this_picSpecies->advanceParticlesIteratively( *m_emfields, a_dt, 
                                                         m_part_order_swap, 
                                                         m_iter_max_particles );
        }
        else {
           this_picSpecies->advanceParticles( *m_emfields, a_dt,
                                              m_part_order_swap ); 
        }
     }
     
     if (m_emfields->advanceE()) {
        if(m_use_mass_matrices) {
           setMassMatrices( a_dt );
           m_emfields->computeJfromMassMatrices();
        }
        else {
           setCurrentDensity();
#ifdef MASS_MATRIX_TEST
           setMassMatrices( a_dt );
           m_emfields->computeJfromMassMatrices();
#endif
        }
     }

  }

  return;    
}

void System::computeRHS( ODEVector<System>&  a_F,
                   const ODEVector<System>&,
                   const Real                a_time,
                   const Real                a_dt,
                   const EMVecType&          a_vec_type )
{
  CH_TIME("System::computeRHS()");
  
  const Real cnormDt = a_dt*m_units->CvacNorm();
  m_emfields->zeroRHS();

  if (a_vec_type == b_only) {

    m_emfields->computeRHSMagneticField( a_time, cnormDt );
    copyBRHSToVec( a_F );

  } else if ( a_vec_type == e_only) {

    m_emfields->computeRHSElectricField( a_time, cnormDt );
    copyERHSToVec( a_F );

  } else if ( a_vec_type == e_and_b ) {

    m_emfields->computeRHSMagneticField( a_time, cnormDt );
    m_emfields->computeRHSElectricField( a_time, cnormDt );
    copyRHSToVec( a_F );

  }

}

void System::computeRHS( ODEVector<System>&  a_F,
                   const ODEVector<System>&,
                   const int                 a_block,
                   const Real                a_time,
                   const Real                a_dt,
                   const int )
{
  CH_TIME("System::computeRHS() from Picard solver");
  
  //m_emfields->zeroRHS();
  const Real cnormDt = a_dt*m_units->CvacNorm();

  if (m_advance_method == PIC_EM_THETA_IMPLICIT) {

    if (a_block == _EM_PICARD_B_FIELD_BLOCK_) {
      m_emfields->computeRHSMagneticField( a_time, cnormDt );
    } else if (a_block == _EM_PICARD_E_FIELD_BLOCK_) {
      m_emfields->computeRHSElectricField( a_time, cnormDt );
    }

  } else if (m_advance_method == PIC_EM_SEMI_IMPLICIT) {

    m_emfields->computeRHSElectricField( a_time, cnormDt );

  }

  copyRHSToVec( a_F );
}

void System::updatePhysicalState(  ODEVector<System>&  a_U,
                          const Real          a_time,
                          const EMVecType&    a_vec_type )
{
  CH_TIME("System::updatePhysicalState()");
  
  if (a_vec_type == b_only) {
    copyBFromVec( a_U );
    m_emfields->applyBCs_magneticField( a_time );
    copyBToVec( a_U );
  } else if (a_vec_type == e_only ) {
    copyEFromVec( a_U );
    m_emfields->applyBCs_electricField( a_time );
    copyEToVec( a_U );
  } else {
    copySolutionFromVec( a_U );
    m_emfields->applyBCs_magneticField( a_time );
    m_emfields->applyBCs_electricField( a_time );
    copySolutionToVec( a_U );
  }
}

void System::updatePhysicalState( ODEVector<System>&  a_U,
                            const int                 a_block,
                            const Real                a_time )
{
  CH_TIME("System::updatePhysicalState() from Picard solver");

  if (m_advance_method == PIC_EM_THETA_IMPLICIT) {

    if (a_block == _EM_PICARD_B_FIELD_BLOCK_) {
      int offset = m_emfields->vecOffsetBField();
      copyBFromVec( a_U, offset );
      m_emfields->applyBCs_magneticField( a_time );
      offset = m_emfields->vecOffsetBField();
      copyBToVec( a_U, offset );
    } else if (a_block == _EM_PICARD_E_FIELD_BLOCK_) {
      int offset = m_emfields->vecOffsetEField();
      copyEFromVec( a_U, offset );
      m_emfields->applyBCs_electricField( a_time );
      offset = m_emfields->vecOffsetEField();
      copyEToVec( a_U, offset );
    }

  } else if (m_advance_method == PIC_EM_SEMI_IMPLICIT) {
    int offset = m_emfields->vecOffsetEField();
    copyEFromVec( a_U, offset );
    m_emfields->applyBCs_electricField( a_time );
    offset = m_emfields->vecOffsetEField();
    copyEToVec( a_U, offset );
  }

}

void System::copyBFromVec(  const ODEVector<System>&  a_vec,
                            int&                      a_offset )
{
  LevelData<FluxBox>& Bfield( m_emfields->getMagneticField() );
  a_offset += SpaceUtils::copyToLevelData( Bfield, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<FArrayBox>& Bfield_virtual( m_emfields->getVirtualMagneticField() );
    a_offset += SpaceUtils::copyToLevelData( Bfield_virtual, a_vec.dataAt(a_offset) );
  }
  return;
}

void System::copyEFromVec( const ODEVector<System>&   a_vec,
                           int&                       a_offset)
{
  //int offset = m_emfields->vecOffsetEField();
  LevelData<EdgeDataBox>& Efield( m_emfields->getElectricField() );
  a_offset += SpaceUtils::copyToLevelData( Efield, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<NodeFArrayBox>& Efield_virtual(m_emfields->getVirtualElectricField());
    a_offset += SpaceUtils::copyToLevelData( Efield_virtual, a_vec.dataAt(a_offset) );
  }
  return;
}

void System::copyBoldFromVec( const ODEVector<System>&  a_vec,
                              int&                      a_offset )
{
  LevelData<FluxBox>& Bfield_old( m_emfields->getMagneticField_old() );
  a_offset += SpaceUtils::copyToLevelData( Bfield_old, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<FArrayBox>& Bfield_virt_old( m_emfields->getVirtualMagneticField_old() );
    a_offset += SpaceUtils::copyToLevelData( Bfield_virt_old, a_vec.dataAt(a_offset) );
  }
  return;
}

void System::copyEoldFromVec( const ODEVector<System>&   a_vec,
                              int&                       a_offset)
{
  LevelData<EdgeDataBox>& Efield_old( m_emfields->getElectricField_old() );
  a_offset += SpaceUtils::copyToLevelData( Efield_old, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<NodeFArrayBox>& Efield_virt_old(m_emfields->getVirtualElectricField_old());
    a_offset += SpaceUtils::copyToLevelData( Efield_virt_old, a_vec.dataAt(a_offset) );
  }
  return;
}

void System::copySolutionFromVec( const ODEVector<System>& a_vec )
{
  CH_TIME("System::copySolutionFromVec()");
  
  int offset(0);
  if (m_advance_method == PIC_EM_THETA_IMPLICIT) {
    copyBFromVec( a_vec, offset );
  }
  copyEFromVec( a_vec, offset );
  return;
}

void System::copyBToVec(  ODEVector<System>&  a_vec,
                          int&                a_offset ) const
{
  const LevelData<FluxBox>& Bfield( m_emfields->getMagneticField() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Bfield );
  if (SpaceDim < 3) {
    const LevelData<FArrayBox>& Bfield_virtual( m_emfields->getVirtualMagneticField() );
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Bfield_virtual );
  }
  return;
}

void System::copyEToVec(  ODEVector<System>&  a_vec,
                          int&                a_offset ) const
{
  //int offset = m_emfields->vecOffsetEField();
  const LevelData<EdgeDataBox>& Efield( m_emfields->getElectricField() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Efield );
  if (SpaceDim < 3) {
    const LevelData<NodeFArrayBox>& Efield_virtual(m_emfields->getVirtualElectricField());
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Efield_virtual );
  }
  return;
}

void System::copyBoldToVec( ODEVector<System>&  a_vec,
                            int&                a_offset ) const
{
  const LevelData<FluxBox>& Bfield_old( m_emfields->getMagneticField_old() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Bfield_old );
  if (SpaceDim < 3) {
    const LevelData<FArrayBox>& Bfield_virt_old( m_emfields->getVirtualMagneticField_old() );
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Bfield_virt_old );
  }
  return;
}

void System::copyEoldToVec( ODEVector<System>&  a_vec,
                            int&                a_offset ) const
{
  const LevelData<EdgeDataBox>& Efield_old( m_emfields->getElectricField_old() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Efield_old );
  if (SpaceDim < 3) {
    const LevelData<NodeFArrayBox>& Efield_virt_old(m_emfields->getVirtualElectricField_old());
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Efield_virt_old );
  }
  return;
}

void System::copySolutionToVec( ODEVector<System>& a_vec ) const
{
  CH_TIME("System::copySolutionToVec()");
  
  int offset(0);
  if (m_advance_method == PIC_EM_THETA_IMPLICIT) {
    copyBToVec( a_vec, offset );
  }
  copyEToVec( a_vec, offset );
  return;
}

void System::copyBRHSFromVec( const ODEVector<System>&  a_vec,
                              int&                      a_offset )
{
  LevelData<FluxBox>& BfieldRHS( m_emfields->getMagneticFieldRHS() );
  a_offset += SpaceUtils::copyToLevelData( BfieldRHS, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<FArrayBox>& BfieldRHS_virtual( m_emfields->getVirtualMagneticFieldRHS() );
    a_offset += SpaceUtils::copyToLevelData( BfieldRHS_virtual, a_vec.dataAt(a_offset) );
  }
  return;
}

void System::copyERHSFromVec( const ODEVector<System>&  a_vec,
                              int&                      a_offset )
{
  LevelData<EdgeDataBox>& EfieldRHS( m_emfields->getElectricFieldRHS() );
  a_offset += SpaceUtils::copyToLevelData( EfieldRHS, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<NodeFArrayBox>& EfieldRHS_virtual(m_emfields->getVirtualElectricFieldRHS());
    a_offset += SpaceUtils::copyToLevelData( EfieldRHS_virtual, a_vec.dataAt(a_offset) );
  }

  return;
}

void System::copyRHSFromVec( const ODEVector<System>& a_vec )
{
  CH_TIME("System::copyRHSFromVec()");
  
  int offset(0);

  if (m_advance_method == PIC_EM_THETA_IMPLICIT) {
    copyBRHSFromVec( a_vec, offset );
  }
  copyERHSFromVec( a_vec, offset );
  return;
}

void System::copyBRHSToVec( ODEVector<System>&  a_vec,
                            int&                a_offset ) const
{
  const LevelData<FluxBox>& BfieldRHS( m_emfields->getMagneticFieldRHS() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), BfieldRHS );
  if (SpaceDim < 3) {
    const LevelData<FArrayBox>& BfieldRHS_virtual( m_emfields->getVirtualMagneticFieldRHS() );
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), BfieldRHS_virtual );
  }
  return;
}

void System::copyERHSToVec( ODEVector<System>&  a_vec,
                            int&                a_offset ) const
{
  const LevelData<EdgeDataBox>& EfieldRHS( m_emfields->getElectricFieldRHS() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), EfieldRHS );
  if (SpaceDim < 3) {
    const LevelData<NodeFArrayBox>& EfieldRHS_virtual(m_emfields->getVirtualElectricFieldRHS());
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), EfieldRHS_virtual );
  }
  return;
}

void System::copyRHSToVec( ODEVector<System>& a_vec ) const
{
  CH_TIME("System::copyRHSToVec()");
  
  int offset(0);
  if (m_advance_method == PIC_EM_THETA_IMPLICIT) {
    copyBRHSToVec( a_vec, offset );
  }
  copyERHSToVec( a_vec, offset );
  return;
}

#include "NamespaceFooter.H"
