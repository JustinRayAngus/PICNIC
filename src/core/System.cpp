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
#include "ScatteringFactory.H"

#include "NamespaceHeader.H"


System::System( ParmParse&  a_pp )
   :
     m_iter_max(0),
     m_iter_max_particles(0),
     m_freeze_particles_jacobian(true),
     m_theta(0.5),
     m_advance_method(PIC_DSMC),
     m_mesh(nullptr),
     m_units(nullptr),
     m_dataFile(nullptr),
     m_meshInterp(nullptr),
     m_verbosity(0),
     m_rtol(1e-6),
     m_atol(1e-12),
     m_time_integrator(nullptr)
{
   ParmParse ppsys("system");
   parseParameters(ppsys);

   m_units = new CodeUnits( ppsys );
   if(!procID()) m_units->printParameters( procID() );
 
   createProblemDomain();          // create the problem domain
  
   DisjointBoxLayout grids;
   getDisjointBoxLayout( grids );  // define the disjointBoxLayout
   
   // initialize the coordinates and grid
   ParmParse ppgrid( "grid" );
   m_mesh = new DomainGrid( ppgrid, m_num_ghosts, m_domain, grids ); 

   m_dataFile = new dataFileIO( ppsys, *m_mesh, *m_units );
     
   createMeshInterp();

   createState( a_pp );
   
   createScattering();

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
                              m_electromagneticFields );
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
}

void System::initialize( const int           a_cur_step,
                         const double        a_cur_time,
                         const std::string&  a_restart_file_name )
{
   CH_TIME("System::initialize()");
   
   // initialize the pic species
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      this_picSpecies->initialize(*m_units,a_cur_time,a_restart_file_name);
   }
   
   // initialize the scattering operators
   if(m_use_scattering) {
      
      // set moments for initial mean free path calculation
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->scatter()) { // prepare this species for scatter
            this_picSpecies->binTheParticles();
            this_picSpecies->setNumberDensity();
            this_picSpecies->setEnergyDensity();
         }
      }
      
      if(!procID()) cout << "Initializing scattering objects..." << endl << endl;
      for (int sct=0; sct<m_scattering_ptr_vect.size(); sct++) {
         ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
         this_scattering->initialize(*m_mesh, m_pic_species_ptr_vect);
      }
      if(!procID()) {
         cout << "Finished initializing " << m_scattering_ptr_vect.size() << " scattering objects" << endl << endl;
      }

   }

   // initialize the electromagnetic fields
   if (!m_electromagneticFields.isNull()) {
      m_electromagneticFields->initialize(a_cur_time,a_restart_file_name);
      setChargeDensity();
      setCurrentDensity();
   }
     
   // assert that m_electromagneticFields is not NULL if using forces 
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->forces()) {
         CH_assert(!m_electromagneticFields.isNull());
         break;
      }
   }

   m_time_integrator->initialize();

}

void System::createProblemDomain()
{
   CH_TIME("System::createProblemDomain()");

   ParmParse ppgrid( "grid" );
   ppgrid.get( "geometry", m_geom_type );
   ppgrid.query( "num_ghosts", m_num_ghosts );
   int this_DIM = SpaceDim;

   if ( m_geom_type == "cylindrical" || m_geom_type == "cartesian") {
     
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

   }
   else {
      // stringstream msg("m_geom_type ",m_geom_type, " not supported");
      // MayDay::Error( msg.str().c_str() );
      if(!procID()) cout << "m_geom_type " << m_geom_type << " not supported " << endl;
      exit(EXIT_FAILURE);
   }
   
   int grid_verbosity;
   ppgrid.query( "verbosity", grid_verbosity );
   if (procID() == 0 && grid_verbosity) {
      cout << "====================== Spatial Grid Parameters =====================" << endl;
      cout << "  geometry = " << m_geom_type << endl;
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

   bool isPeriodic[SpaceDim];
   for (int dir=0; dir<SpaceDim; ++dir) isPeriodic[dir] = (m_is_periodic[dir] == 1);

   m_domain.define( level0Domain.smallEnd(),
                    level0Domain.bigEnd(),
                    isPeriodic );

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
      m_electromagneticFields = RefCountedPtr<ElectroMagneticFields>
                                  (new ElectroMagneticFields( ppflds,
                                                              *m_mesh,
                                                              *m_units,
                                                              verbose,
                                                              em_vec_type)  ); 
      if(m_advance_method == PIC_DSMC ) {
         if(!procID()) cout << "advance_method PIC_DSMC cannot be used with electromagnetic fields on" << endl;
         exit(EXIT_FAILURE);
      }
      if(!procID()) cout << "Finished creating Electromagnetic fields object" << endl << endl;
   }
   else {
      m_electromagneticFields = RefCountedPtr<ElectroMagneticFields>(NULL); // don't have to do this, but for clarity 
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


void System::createScattering()
{
   // Create all scattering objects
   
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) { // check all species for scattering flag
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      m_use_scattering = this_picSpecies->scatter();
      if(m_use_scattering) break;
   }

   if(m_use_scattering) {

      ScatteringFactory scatteringFactory;
   
      if(!procID()) {
         cout << "Creating scattering objects..." << endl;
      }

      bool more_scattering_models(true);
      while(more_scattering_models) { // look for scattering ops
      
         stringstream s;
         s << "scattering." << m_scattering_ptr_vect.size(); 
         ParmParse pp_scatter( s.str().c_str() );
     
         string model;
         if(pp_scatter.contains("model")) {
         
            // Create this scattering object
            ScatteringPtr this_scattering = scatteringFactory.create( pp_scatter, 1 );
            
            // Add the new scattering object to the vector of pointers to scattering objects
            m_scattering_ptr_vect.push_back(this_scattering);

         } 
         else {
            more_scattering_models = false;
         }
      
      } 
      
      if(!procID()) {
         cout << "Finished creating " << m_scattering_ptr_vect.size() << " scattering objects" << endl << endl;
      }
      
   }
   else {
      if(!procID()) {
         cout << "Scattering flag = false for all pic species" << endl << endl;
      }
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
                            const double  a_cur_time )
{
   CH_TIME("System::writePlotFile()");
   
   if (!m_electromagneticFields.isNull()) {
      if(m_electromagneticFields->writeRho()) setChargeDensity();
      if(m_electromagneticFields->writeDivs()) {
         m_electromagneticFields->setDivE();
         m_electromagneticFields->setDivB();
      }
      m_dataFile->writeElectroMagneticFieldsDataFile( *m_electromagneticFields,
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

void System::writeHistFile( const int     a_cur_step,
                            const double  a_cur_time,
                            const bool    a_startup )
{
   CH_TIME("System::writeHistFile()");
   
   if(a_startup) setupHistFile(a_cur_step);
   
   // compute the field probes
   Real energyE_joules, energyB_joules;
   if(m_field_probes) {
      energyE_joules = m_electromagneticFields->electricFieldEnergy();
      energyB_joules = m_electromagneticFields->magneticFieldEnergy();
   }
   
   // compute the species probes
   const int numSpecies = m_pic_species_ptr_vect.size();
   std::vector<int> num_particles(numSpecies);
   std::vector<vector<Real>> global_moments(numSpecies);
   for( int species=0; species<numSpecies; species++) {
      if(!m_species_probes[species]) continue;
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[species]);
      num_particles[species] = this_picSpecies->numParticles();
      this_picSpecies->globalMoments(global_moments.at(species));
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
         if(m_field_probes) {
            histFile << energyE_joules << " ";
            histFile << energyB_joules << " ";
         }
         for( int species=0; species<numSpecies; species++) {
            if(!m_species_probes[species]) continue;
            histFile << num_particles[species] << " ";
            std::vector<Real>& species_moments = global_moments.at(species);
            for(int mom=0; mom<species_moments.size(); ++mom) {
               histFile << species_moments.at(mom) << " ";
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

   // parse input for different probes to include in the history file
   ParmParse pphist("history");
   
   m_history_precision = 5;
   pphist.query("precision", m_history_precision);
   
   m_field_probes = false;
   if(!m_electromagneticFields.isNull()) {
      pphist.query("field_probes", m_field_probes);
      if(!procID()) cout << "field_probes = " << m_field_probes << endl;
   }
 
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

   if(m_field_probes) {
      stringstream ss;
      ss << "#" << probe_names.size() << " electric field energy [Joules]";
      probe_names.push_back(ss.str()); 
      ss.str(std::string());
      ss << "#" << probe_names.size() << " magnetic field energy [Joules]";
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
   
   // write the field data
   if (!m_electromagneticFields.isNull()) {
      m_dataFile->writeCheckpointEMFields( a_handle, *m_electromagneticFields );
   }

   // write the particle data
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      m_dataFile->writeCheckpointParticles( a_handle, *this_picSpecies, s );
   }

}

void System::readCheckpointFile( const std::string&  a_chkpt_fname,
                                 int&                a_cur_step,
                                 Real&               a_cur_time,
                                 Real&               a_cur_dt )
{
   CH_TIME("System::readCheckpointFile()");

   if(!procID()) cout << "Reading restart file " << a_chkpt_fname << endl;

#ifdef CH_USE_HDF5
   HDF5Handle handle( a_chkpt_fname, HDF5Handle::OPEN_RDONLY );

   HDF5HeaderData header;
   header.readFromFile( handle );

   a_cur_step = header.m_int["cur_step"];
   a_cur_time = header.m_real["cur_time"];
   a_cur_dt   = header.m_real["cur_dt"];

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
      a_ppsys.query("abs_tol", m_atol);
      a_ppsys.query("rel_tol", m_rtol);
   } else if(m_advance_method == PIC_EM_THETA_IMPLICIT) {
      a_ppsys.get("iter_max",m_iter_max);
      a_ppsys.query("iter_max_particles",m_iter_max_particles);
      a_ppsys.query("freeze_particles_jacobian",m_freeze_particles_jacobian);
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
         cout << "  abs_tol  = " << m_atol << endl;
         cout << "  rel_tol  = " << m_rtol << endl;
      } else if(m_advance_method == PIC_EM_THETA_IMPLICIT) {
         cout << "  iter_max = " << m_iter_max << endl;
         cout << "  iter_max_particles = " << m_iter_max_particles << endl;
         cout << "  freeze_particles_jacobian = " << m_freeze_particles_jacobian << endl;
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
 
   LevelData<NodeFArrayBox>& chargeDensity = m_electromagneticFields->getChargeDensity();
   
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
 
   LevelData<EdgeDataBox>&  currentDensity = m_electromagneticFields->getCurrentDensity();
   LevelData<NodeFArrayBox>&  currentDensity_virtual = m_electromagneticFields->getVirtualCurrentDensity();
   
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
            }
         }

      }
   }

   // call exchange
   //SpaceUtils::exchangeEdgeDataBox(m_currentDensity);
   //if(SpaceDim<3) SpaceUtils::exchangeNodeFArrayBox(m_currentDensity_virtual,m_mesh);

}

void System::scatterParticles( const Real&  a_dt )
{  
   if (!m_use_scattering) return;

   CH_TIME("System::scatterParticles()");

   const Real tscale = m_units->getScale(m_units->TIME);
   const Real dt_sec = a_dt*tscale;      

   // prepare all species to be scattered for scattering
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->scatter()) {
         this_picSpecies->binTheParticles();
         this_picSpecies->setNumberDensity(); // only need to set this once as collisions don't change it
      }
   }

   // Should we shuffle the m_scattering_ptr_vect each time ?
   for (int sct=0; sct<m_scattering_ptr_vect.size(); sct++) { // loop over scattering objects
         
      ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
      const int this_species1 = this_scattering->species1();
      const int this_species2 = this_scattering->species2();
         
      PicSpeciesPtr this_picSpecies1(m_pic_species_ptr_vect[this_species1]);
      //this_picSpecies1->setEnergyDensity();
      if(this_species1==this_species2) { // self-species scattering
         this_scattering->applySelfScattering( *this_picSpecies1, *m_mesh, dt_sec );
      }
      else { // inter-species scattering
         PicSpeciesPtr this_picSpecies2(m_pic_species_ptr_vect[this_species2]);
         //this_picSpecies2->setEnergyDensity();
         this_scattering->applyInterScattering( *this_picSpecies1, *this_picSpecies2, *m_mesh, dt_sec );
      }

   }

}   

void System::postTimeStep( Real&        a_cur_time,
                           const Real&  a_dt,
                           int&         a_step_number )
{  
   CH_TIME("System::postTimeStep()");

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
   if (!m_electromagneticFields.isNull() && m_electromagneticFields->advance()) {
      fldsDt = m_electromagneticFields->stableDt();
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
  
   if(m_pic_species_ptr_vect.size()>0) { 
      if(!procID()) cout << "stable particle time step = " << stableDt << endl;
   }
   return stableDt;

}


Real System::scatterDt( const int a_step_number )
{
   Real scatterDt = DBL_MAX;
   const Real tscale = m_units->getScale(m_units->TIME);
   
   if(m_use_scattering) {
  
      // compute the density and energy moments for mean free time calculation
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->scatter()) {
            this_picSpecies->setNumberDensity();
            this_picSpecies->setEnergyDensity();
         }
      }

      for (int sct=0; sct<m_scattering_ptr_vect.size(); sct++) {

         ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
         int sp1 = this_scattering->species1();
         PicSpeciesPtr this_picSpecies1 = m_pic_species_ptr_vect[sp1];
         const LevelData<FArrayBox>& numberDensity1 = this_picSpecies1->getNumberDensity();
         const LevelData<FArrayBox>& energyDensity1 = this_picSpecies1->getEnergyDensity();
         
         int sp2 = this_scattering->species2();
         if(sp1==sp2) {
            this_scattering->setMeanFreeTime(numberDensity1,energyDensity1);
         }
         else {
            PicSpeciesPtr this_picSpecies2 = m_pic_species_ptr_vect[sp2];
            const LevelData<FArrayBox>& numberDensity2 = this_picSpecies2->getNumberDensity();
            const LevelData<FArrayBox>& energyDensity2 = this_picSpecies2->getEnergyDensity();
            this_scattering->setMeanFreeTime( numberDensity1, energyDensity1,
                                              numberDensity2, energyDensity2 );
         }

         scatterDt = Min(scatterDt,this_scattering->scatterDt()); // [s]
         scatterDt = scatterDt/tscale; // convert to code units

      }
      if(!procID()) cout << "mean free scattering time = " << scatterDt << endl;
   
   }
   
   return scatterDt;
}


Real System::specialOpsDt( const int a_step_number )
{
   Real specialOpsDt = DBL_MAX;
   if(m_use_specialOps) {
      specialOpsDt = m_specialOps->specialOpsDt();
      if(!procID()) cout << "special operator time step = " << specialOpsDt << endl;
   }
   return specialOpsDt;
}

void System::adaptDt( bool&  a_adapt_dt )
{
   if(m_advance_method==PIC_EM_EXPLICIT) a_adapt_dt = false;
   if(m_advance_method==PIC_EM_SEMI_IMPLICIT) a_adapt_dt = false;
}

void System::preRHSOp( const ODEVector<System>&  a_U,
                       const Real                a_dt )
{  
  CH_TIME("System::preRHSOp()");
  
  CH_assert(!m_electromagneticFields.isNull());
    
  // half dt advance of particle positions and velocities
  // and, if advancing E, compute current density
  //
  // xbar = xn + dt/2*vbar
  // vbar = vn + dt/2*q/m*(E(xbar) + vbar x B(xbar))
  // Jbar = Sum_sSum_p(qp*S(xbar-xg)*vbar)/dV
   
  for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
     auto this_picSpecies(m_pic_species_ptr_vect[s]);
     this_picSpecies->repositionOutcastsAndApplyForces( *m_electromagneticFields, 
                                                        a_dt, true );
     if(m_iter_max_particles>0) {
        this_picSpecies->advanceParticlesIteratively( *m_electromagneticFields, 
                                                      a_dt, m_iter_max_particles );
     }
     else {
        this_picSpecies->advancePositions(0.5*a_dt,true);
        this_picSpecies->interpolateFieldsToParticles( *m_electromagneticFields );
        this_picSpecies->advanceVelocities( a_dt, true );
     }
  }
  if (m_electromagneticFields->advanceE()) setCurrentDensity();

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
  m_electromagneticFields->zeroRHS();

  if (a_vec_type == b_only) {

    m_electromagneticFields->computeRHSMagneticField( a_time, cnormDt );
    copyBRHSToVec( a_F );

  } else if ( a_vec_type == e_only) {

    m_electromagneticFields->computeRHSElectricField( a_time, cnormDt );
    copyERHSToVec( a_F );

  } else if ( a_vec_type == e_and_b ) {

    m_electromagneticFields->computeRHSMagneticField( a_time, cnormDt );
    m_electromagneticFields->computeRHSElectricField( a_time, cnormDt );
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
  
  //m_electromagneticFields->zeroRHS();
  const Real cnormDt = a_dt*m_units->CvacNorm();

  if (m_advance_method == PIC_EM_THETA_IMPLICIT) {

    if (a_block == _EM_PICARD_B_FIELD_BLOCK_) {
      m_electromagneticFields->computeRHSMagneticField( a_time, cnormDt );
    } else if (a_block == _EM_PICARD_E_FIELD_BLOCK_) {
      m_electromagneticFields->computeRHSElectricField( a_time, cnormDt );
    }

  } else if (m_advance_method == PIC_EM_SEMI_IMPLICIT) {

    m_electromagneticFields->computeRHSElectricField( a_time, cnormDt );

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
    m_electromagneticFields->updateDataMagneticField( a_time );
    copyBToVec( a_U );
  } else if (a_vec_type == e_only ) {
    copyEFromVec( a_U );
    m_electromagneticFields->updateDataElectricField( a_time );
    copyEToVec( a_U );
  } else {
    copySolutionFromVec( a_U );
    m_electromagneticFields->updateDataMagneticField( a_time );
    m_electromagneticFields->updateDataElectricField( a_time );
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
      int offset = m_electromagneticFields->vecOffsetBField();
      copyBFromVec( a_U, offset );
      m_electromagneticFields->updateDataMagneticField( a_time );
      offset = m_electromagneticFields->vecOffsetBField();
      copyBToVec( a_U, offset );
    } else if (a_block == _EM_PICARD_E_FIELD_BLOCK_) {
      int offset = m_electromagneticFields->vecOffsetEField();
      copyEFromVec( a_U, offset );
      m_electromagneticFields->updateDataElectricField( a_time );
      offset = m_electromagneticFields->vecOffsetEField();
      copyEToVec( a_U, offset );
    }

  } else if (m_advance_method == PIC_EM_SEMI_IMPLICIT) {
    int offset = m_electromagneticFields->vecOffsetEField();
    copyEFromVec( a_U, offset );
    m_electromagneticFields->updateDataElectricField( a_time );
    offset = m_electromagneticFields->vecOffsetEField();
    copyEToVec( a_U, offset );
  }

}

void System::copyBFromVec(  const ODEVector<System>&  a_vec,
                            int&                      a_offset )
{
  LevelData<FluxBox>& Bfield( m_electromagneticFields->getMagneticField() );
  a_offset += SpaceUtils::copyToLevelData( Bfield, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<FArrayBox>& Bfield_virtual( m_electromagneticFields->getVirtualMagneticField() );
    a_offset += SpaceUtils::copyToLevelData( Bfield_virtual, a_vec.dataAt(a_offset) );
  }
  return;
}

void System::copyEFromVec( const ODEVector<System>&   a_vec,
                           int&                       a_offset)
{
  //int offset = m_electromagneticFields->vecOffsetEField();
  LevelData<EdgeDataBox>& Efield( m_electromagneticFields->getElectricField() );
  a_offset += SpaceUtils::copyToLevelData( Efield, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<NodeFArrayBox>& Efield_virtual(m_electromagneticFields->getVirtualElectricField());
    a_offset += SpaceUtils::copyToLevelData( Efield_virtual, a_vec.dataAt(a_offset) );
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
  const LevelData<FluxBox>& Bfield( m_electromagneticFields->getMagneticField() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Bfield );
  if (SpaceDim < 3) {
    const LevelData<FArrayBox>& Bfield_virtual( m_electromagneticFields->getVirtualMagneticField() );
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Bfield_virtual );
  }
  return;
}

void System::copyEToVec(  ODEVector<System>&  a_vec,
                          int&                a_offset ) const
{
  //int offset = m_electromagneticFields->vecOffsetEField();
  const LevelData<EdgeDataBox>& Efield( m_electromagneticFields->getElectricField() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Efield );
  if (SpaceDim < 3) {
    const LevelData<NodeFArrayBox>& Efield_virtual(m_electromagneticFields->getVirtualElectricField());
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), Efield_virtual );
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
  LevelData<FluxBox>& BfieldRHS( m_electromagneticFields->getMagneticFieldRHS() );
  a_offset += SpaceUtils::copyToLevelData( BfieldRHS, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<FArrayBox>& BfieldRHS_virtual( m_electromagneticFields->getVirtualMagneticFieldRHS() );
    a_offset += SpaceUtils::copyToLevelData( BfieldRHS_virtual, a_vec.dataAt(a_offset) );
  }
  return;
}

void System::copyERHSFromVec( const ODEVector<System>&  a_vec,
                              int&                      a_offset )
{
  LevelData<EdgeDataBox>& EfieldRHS( m_electromagneticFields->getElectricFieldRHS() );
  a_offset += SpaceUtils::copyToLevelData( EfieldRHS, a_vec.dataAt(a_offset) );
  if (SpaceDim < 3) {
    LevelData<NodeFArrayBox>& EfieldRHS_virtual(m_electromagneticFields->getVirtualElectricFieldRHS());
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
  const LevelData<FluxBox>& BfieldRHS( m_electromagneticFields->getMagneticFieldRHS() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), BfieldRHS );
  if (SpaceDim < 3) {
    const LevelData<FArrayBox>& BfieldRHS_virtual( m_electromagneticFields->getVirtualMagneticFieldRHS() );
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), BfieldRHS_virtual );
  }
  return;
}

void System::copyERHSToVec( ODEVector<System>&  a_vec,
                            int&                a_offset ) const
{
  const LevelData<EdgeDataBox>& EfieldRHS( m_electromagneticFields->getElectricFieldRHS() );
  a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), EfieldRHS );
  if (SpaceDim < 3) {
    const LevelData<NodeFArrayBox>& EfieldRHS_virtual(m_electromagneticFields->getVirtualElectricFieldRHS());
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
