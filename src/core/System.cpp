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
     m_advance_method(PIC_DSMC),
     m_mesh(nullptr),
     m_units(nullptr),
     m_dataFile(nullptr),
     m_pic_species(nullptr),
     m_scattering(nullptr),
     m_verbosity(0),
     m_use_specialOps(false),
     m_time_integrator(nullptr)
{
   ParmParse ppsys("system");
   parseParameters(ppsys);

   createProblemDomain();
  
   DisjointBoxLayout grids;
   getDisjointBoxLayout( grids );
   
   // initialize the coordinates and grid
   ParmParse ppgrid( "grid" );
   m_mesh = new DomainGrid( ppgrid, m_num_ghosts, m_domain, grids ); 
   
   m_units = new CodeUnits( ppsys, m_mesh->axisymmetric() );
   m_units->printParameters();

   m_dataFile = new dataFileIO( ppsys, *m_mesh, *m_units );
     
   createState();
   
   m_scattering = new ScatteringInterface( *m_units, m_pic_species->getPtrVect() ); 

   createSpecialOperators();

   m_implicit_advance = false;
   bool em_advance = false;
   if(m_advance_method == PIC_DSMC) {
     m_time_integrator = new PICTimeIntegrator_DSMC;
   } else if(m_advance_method == PIC_EM_EXPLICIT) {
     m_time_integrator = new PICTimeIntegrator_EM_Explicit;
     em_advance = true;
   } else if (m_advance_method == PIC_EM_SEMI_IMPLICIT) {
     m_time_integrator = new PICTimeIntegrator_EM_SemiImplicit;
     m_implicit_advance = true;
     em_advance = true;
   } else if(m_advance_method == PIC_EM_THETA_IMPLICIT) {
     m_time_integrator = new PICTimeIntegrator_EM_ThetaImplicit;
     m_implicit_advance = true;
     em_advance = true;
   }
   else {
      if(!procID()) { 
         cout << "EXIT FAILURE!!!" << endl;
         cout << "m_advance_method = " << m_advance_method;
         cout << " is not a valid option" << endl;
         cout << " valid options: " << PIC_DSMC << ", ";
         cout <<  PIC_EM_EXPLICIT << ", ";
         cout <<  PIC_EM_SEMI_IMPLICIT << ", ";
         cout <<  PIC_EM_THETA_IMPLICIT << " " << endl;
      }
      exit(EXIT_FAILURE);
   }
   if(em_advance && m_emfields.isNull()) {
      if(!procID()) { 
         cout << "EXIT FAILURE!!!" << endl;
         cout << "m_advance_method = " << m_advance_method;
         cout << " requires em_fields.use = true" << endl << endl;
      }
      exit(EXIT_FAILURE);
   }
   m_time_integrator->define( this,
                              m_pic_species->getPtrVect(), 
                              m_emfields ); // m_emfields can't be NULL for implicit solvers
}

System::~System()
{
   delete m_time_integrator;
   delete m_mesh;
   delete m_units;
   delete m_dataFile;
   delete m_pic_species;
   delete m_scattering;
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

void System::createState()
{

   m_pic_species = new PicSpeciesInterface( *m_units, *m_mesh );

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
                                                              em_vec_type )); 
      
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

void System::createSpecialOperators()
{
   // Create the vector of special operator (pointers)
   
   SpecialOperatorFactory  specialOpFactory;
   
   if(!procID()) cout << "Adding special operators..." << endl;

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

   if(!procID()) cout << "Done adding special operators" << endl << endl;


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
   m_pic_species->initialize( *m_units, m_implicit_advance, a_cur_time, a_restart_file_name );
   
   // initialize the scattering operators
   m_scattering->initialize( m_pic_species->getPtrVect(), *m_mesh, a_restart_file_name );
      
   // initialize the electromagnetic fields
   if(m_emfields.isNull()) {CH_assert(!m_pic_species->forces());}
   else {
      setChargeDensity();
      m_emfields->initialize(a_cur_time,a_restart_file_name);
      m_pic_species->setCurrentDensity();
      setCurrentDensity();
   }
     
   m_time_integrator->initialize(a_restart_file_name);
   
}

void System::writePlotFile( const int     a_cur_step,
                            const double  a_cur_time,
                            const double  a_cur_dt )
{
   CH_TIME("System::writePlotFile()");
   
   if (!m_emfields.isNull()) {
      if(m_emfields->writeRho()) setChargeDensity();
      if(m_emfields->writeExB()) m_emfields->setPoyntingFlux();
      if(m_emfields->writeDivs()) {
         m_emfields->setDivE();
         m_emfields->setDivB();
      }
      m_dataFile->writeElectroMagneticFieldsDataFile( *m_emfields,
#ifdef MASS_MATRIX_TEST
                                                      *m_pic_species,
#endif
                                                      a_cur_step, a_cur_time ); 
   }

   m_dataFile->writePicSpecies( *m_pic_species, a_cur_step, a_cur_time);

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
   if(m_solver_probes) {
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
   const int numSpecies = m_pic_species->numSpecies();
   std::vector<int> num_particles(numSpecies);
   std::vector<vector<Real>> global_moments(numSpecies);
   std::vector<vector<Real>> bdry_moments(numSpecies);
   std::vector<vector<Real>> species_solver(numSpecies);
   std::vector<Real> max_wpdt(numSpecies);
   std::vector<Real> max_wcdt(numSpecies);
   PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
   for( int sp=0; sp<numSpecies; sp++) {
      if(!m_species_probes[sp]) continue;
      PicSpeciesPtr species(pic_species_ptr_vect[sp]);
      num_particles[sp] = species->numParticles();
      species->globalMoments(global_moments.at(sp));
      if(m_species_bdry_probes) species->bdryMoments(bdry_moments.at(sp));
      max_wpdt.at(sp) = species->max_wpdt(*m_units,a_cur_dt);
      max_wcdt.at(sp) = max_wc0dt*abs(species->charge())/species->mass();
      if(m_species_solver_probes) species->picardParams(species_solver.at(sp));
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
         for( int sp=0; sp<numSpecies; sp++) {
            if(!m_species_probes[sp]) continue;
            histFile << num_particles[sp] << " ";
            std::vector<Real>& species_moments = global_moments.at(sp);
            for(int mom=0; mom<species_moments.size(); ++mom) {
               histFile << species_moments.at(mom) << " ";
            }
            histFile << max_wpdt.at(sp) << " ";
            histFile << max_wcdt.at(sp) << " ";
            if(m_species_bdry_probes) {
               std::vector<Real>& this_bdry_moment = bdry_moments.at(sp);
               const int num_per_dir = this_bdry_moment.size()/SpaceDim;
               for(int dir=0; dir<SpaceDim; ++dir) { 
                  if(m_is_periodic[dir]) continue;
                  for(int mom=0; mom<num_per_dir; ++mom) {
                     histFile << this_bdry_moment.at(mom+dir*num_per_dir) << " ";
                  }
               }
            }
            if(m_species_solver_probes) {
               std::vector<Real>& this_species_solver = species_solver.at(sp);
               for(int its=0; its<this_species_solver.size(); ++its) {
                  histFile << this_species_solver.at(its) << " ";
               }
            }
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
 
   const int numSpecies = m_pic_species->numSpecies();
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
      m_species_bdry_probes = false;
      pphist.query("species_bdry_probes", m_species_bdry_probes);
      if(!procID()) cout << "species_bdry_probes = " << m_species_bdry_probes << endl;
      m_species_solver_probes = false;
      pphist.query("species_solver_probes", m_species_solver_probes);
      if(!procID()) cout << "species_solver_probes = " << m_species_solver_probes << endl;
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
      if(m_species_bdry_probes) {
      for (int dir=0; dir<SpaceDim; dir++) {
         if(m_is_periodic[dir]) continue;
         stringstream ss;
         ss << "#" << probe_names.size() << " species " << species << ": mass out from lo-side, dir = " << dir << " [kg]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": mass out from hi-side, dir = " << dir << " [kg]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": X-momentum out from lo-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": X-momentum out from hi-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": Y-momentum out from lo-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": Y-momentum out from hi-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": Z-momentum out from lo-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": Z-momentum out from hi-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": energy out from lo-side, dir = " << dir << " [Joules]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": energy out from hi-side, dir = " << dir << " [Joules]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": mass in from lo-side, dir = " << dir << " [kg]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": mass in from hi-side, dir = " << dir << " [kg]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": X-momentum in from lo-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": X-momentum in from hi-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": Y-momentum in from lo-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": Y-momentum in from hi-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": Z-momentum in from lo-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": Z-momentum in from hi-side, dir = " << dir << " [kg-m/s]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": energy in from lo-side, dir = " << dir << " [Joules]";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": energy in from hi-side, dir = " << dir << " [Joules]";
         probe_names.push_back(ss.str()); 
      }
      }
      if(m_species_solver_probes) {
         ss << "#" << probe_names.size() << " species " << species << ": avg picard its";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
         ss << "#" << probe_names.size() << " species " << species << ": max avg picard its";
         probe_names.push_back(ss.str()); 
         ss.str(std::string());
      }
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
   m_dataFile->writeCheckpointPicSpecies( a_handle,  *m_pic_species);
   
   // write solver probes
   int l_exit_status=0, nl_exit_status=0;
   int l_total_iter=0, nl_total_iter=0;
   int l_last_iter=0, nl_iter=0;
   Real nl_abs_res=0.0, nl_rel_res=0.0;
   if(m_solver_probes) {
      m_time_integrator->getConvergenceParams( l_exit_status, l_last_iter, l_total_iter,
                                              nl_exit_status, nl_iter, nl_total_iter,
                                              nl_abs_res, nl_rel_res );
   }
   header.clear();
   const std::string solverGroup = std::string("solver");
   a_handle.setGroup(solverGroup);

   header.m_int["l_exit_status"] = l_exit_status;
   header.m_int["l_last_iter"] = l_last_iter;
   header.m_int["l_total_iter"] = l_total_iter;
   header.m_int["nl_exit_status"] = nl_exit_status;
   header.m_int["nl_iter"] = nl_iter;
   header.m_int["nl_total_iter"] = nl_total_iter;
   header.m_real["nl_abs_res"] = nl_abs_res;
   header.m_real["nl_rel_res"] = nl_rel_res;

   header.writeToFile( a_handle );

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
   if(m_advance_method == "DSMC") m_advance_method = PIC_DSMC;
   if(m_advance_method == PICMC_EXPLICIT) m_advance_method = PIC_EM_EXPLICIT;
   if(m_advance_method == PICMC_SEMI_IMPLICIT) m_advance_method = PIC_EM_SEMI_IMPLICIT;
   if(m_advance_method == PICMC_FULLY_IMPLICIT) m_advance_method = PIC_EM_THETA_IMPLICIT;

   if(!procID()) cout << "advance method  = " << m_advance_method << endl << endl;
 
}

void System::preTimeStep( const Real  a_time,
                          const Real  a_dt,
                          const int   a_step_number )
{  
   CH_TIME("System::preTimeStep()");
   m_time_integrator->preTimeStep( a_time, a_dt, a_step_number );
}

void System::timeStep( const Real  a_cur_time,
                       const Real  a_dt,
                       const int   a_step_number )
{  
   CH_TIME("System::timeStep()");
   m_time_integrator->timeStep( a_cur_time, a_dt, a_step_number );
}

void System::postTimeStep( Real&        a_cur_time,
                           const Real&  a_dt,
                           int&         a_step_number )
{  
   CH_TIME("System::postTimeStep()");
   
   if(!m_emfields.isNull()) m_emfields->postTimeStep( a_dt );

   m_time_integrator->postTimeStep( a_cur_time, a_dt );

   // apply special operators
   if(m_specialOps) {
      m_specialOps->applyOp(m_pic_species->getPtrVect(),a_dt);
      m_specialOps->updateOp(a_dt);
   }
   
   m_pic_species->postTimeStep();
   
   a_cur_time = a_cur_time + a_dt;
   a_step_number = a_step_number + 1;
   m_nonlinear_iter = 0;
}

void System::setChargeDensity()
{
   CH_TIME("System::setChargeDensity()");
 
   // set the pic charge density
   m_pic_species->setChargeDensityOnNodes();
   const LevelData<NodeFArrayBox>& pic_rho = m_pic_species->getChargeDensityOnNodes();
   
   // add it to the total charge density
   LevelData<NodeFArrayBox>& rho = m_emfields->getChargeDensity();
   const DisjointBoxLayout& grids(m_mesh->getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& this_rho = rho[dit].getFab();
      const FArrayBox& this_pic_rho = pic_rho[dit].getFab();
      this_rho.copy(this_pic_rho,0,0,rho.nComp());
   }

}

void System::setCurrentDensity( const bool  a_from_explicit_solver )
{
   CH_TIME("System::setCurrentDensity()");
   
   // set the pic species current density
   if(m_advance_method==PIC_EM_EXPLICIT) m_pic_species->setCurrentDensity( true );
   const LevelData<EdgeDataBox>& pic_J = m_pic_species->getCurrentDensity();
   const LevelData<NodeFArrayBox>& pic_Jv = m_pic_species->getVirtualCurrentDensity();
   
   // add (copy) it to the total current density used to advance E
   LevelData<EdgeDataBox>& J = m_emfields->getCurrentDensity();
   LevelData<NodeFArrayBox>& Jv = m_emfields->getVirtualCurrentDensity();
   const DisjointBoxLayout& grids(m_mesh->getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& this_J = J[dit][dir];
         const FArrayBox& this_pic_J = pic_J[dit][dir];
         this_J.copy(this_pic_J);
      }
#if CH_SPACEDIM<3
      FArrayBox& this_Jv = Jv[dit].getFab();
      const FArrayBox& this_pic_Jv = pic_Jv[dit].getFab();
      this_Jv.copy(this_pic_Jv,0,0,Jv.nComp());
#endif
   }

}

void System::scatterParticles( const Real&  a_dt )
{  
   if (m_scattering->numScatter()==0) return;

   CH_TIME("System::scatterParticles()");

   const Real tscale = m_units->getScale(m_units->TIME);
   const Real dt_sec = a_dt*tscale;      

   m_pic_species->prepForScatter();
   m_scattering->applyScattering( m_pic_species->getPtrVect(), *m_mesh, dt_sec );
   
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
   Real partsDt = DBL_MAX;
   partsDt = m_pic_species->courantDt();
   return partsDt;
}

Real System::scatterDt( const int a_step_number )
{

   Real scatterDt = DBL_MAX;
   
   if(m_scattering->numScatter()>0) {

      scatterDt = m_scattering->scatterDt( m_pic_species->getPtrVect() );
 
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
                       const int                 a_nonlinear_iter,
                       const bool                a_from_emjacobian )
{  
  CH_TIME("System::preRHSOp()");

  m_pic_species->preRHSOp( a_from_emjacobian, *m_emfields, a_dt, m_nonlinear_iter );

  setCurrentDensity();

  // a_nonlinear_iter when using Petsc comes here with zero twice, and passing it
  // to m_pic_species->preRHSOp() doesn't work with new usage of this parameter

  //if(!a_from_emjacobian && !procID()) {
     //cout << "JRA: m_nonlinear_iter = " << m_nonlinear_iter << endl;
     //cout << "JRA: a_nonlinear_iter = " << a_nonlinear_iter << endl << endl;
  //}
  
  if(!a_from_emjacobian) m_nonlinear_iter += 1;

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
