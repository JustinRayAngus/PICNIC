#include "System.H"
#include "CH_Timer.H"

#include "BoxIterator.H"
#include "DomainGrid.H"

#include "dataFileIO.H"
#include <iostream>
#include <fstream>

#include "SpecialOperatorFactory.H"
#include "ScatteringFactory.H"

#include "NamespaceHeader.H"


System::System( ParmParse&  a_pp )
   :
     m_iter_max(0),
     m_theta(0.5),
     m_mesh(NULL),
     m_units(NULL),
     m_dataFile(NULL),
     m_meshInterp(NULL),
     m_verbosity(0),
     m_rtol(1e-6),
     m_atol(1e-12)
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

}

System::~System()
{
   delete m_mesh;
   delete m_units;
   delete m_dataFile;
   if(m_meshInterp!=NULL) {
      delete m_meshInterp;
      m_meshInterp = NULL;
   }
}

void System::initialize( const int     a_cur_step,
                         const double  a_cur_time )
{
   CH_TIME("System::initialize()");
   
   // initialize the pic species
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(a_cur_step==0) {
         this_picSpecies->initialize(*m_units); // set initial particle positions and velocities
      }
      else { // restart
         if(!procID()) cout << "System::initialize() called at step = " << a_cur_step << endl;
         if(!procID()) cout << "restarting not currently an option !!! " << endl;
         exit(EXIT_FAILURE);
      }
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
      m_electromagneticFields->initialize();
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         int this_charge = this_picSpecies->charge();
         if(this_charge != 0) this_picSpecies->setCurrentDensity();
      }
      m_electromagneticFields->setCurrentDensity(m_pic_species_ptr_vect);
   }
     
   // assert that m_electromagneticFields is not NULL if using forces 
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->forces()) {
         CH_assert(!m_electromagneticFields.isNull());
         break;
      }
   }

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
   
   bool use_fields = false;
   ppflds.query("use",use_fields);
   if(use_fields) {
      bool verbose = true;
      m_electromagneticFields = RefCountedPtr<ElectroMagneticFields>(new ElectroMagneticFields(ppflds,*m_mesh,*m_units,verbose)); 
      if(m_advance_method == DSMC ) {
         if(!procID()) cout << "advance_method DSMC cannot be used with electromagnetic fields on" << endl;
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
 
      species = m_pic_species_ptr_vect.size() + 1;

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
         PicSpecies* picSpecies = new PicSpecies( ppspc, name, *m_meshInterp, *m_mesh );

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
         s << "scattering." << m_scattering_ptr_vect.size() + 1; 
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
   int special_op_num = 0;
   while(more_ops) { // look for special operator...
 
      special_op_num = special_op_num + 1;

      stringstream s;
      s << "special_operator." << special_op_num; 
      ParmParse ppspop( s.str().c_str() );
     
      if(ppspop.contains("model")) {
         m_use_specialOps = true;
         m_specialOps = specialOpFactory.create( ppspop, 1 );
         m_specialOps->updateOp(*m_mesh,*m_units,0.0); // should make an initializeOp() function
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
                                                  s+1, a_cur_step, a_cur_time );
      }
      else {
         if(m_writeSpeciesChargeDensity) {
            this_picSpecies->setChargeDensity();      
            this_picSpecies->setChargeDensityOnFaces(); 
         }     
         if(m_writeSpeciesCurrentDensity) this_picSpecies->setCurrentDensity();      
         m_dataFile->writeChargedSpeciesDataFile( *this_picSpecies, 
                                                  s+1, a_cur_step, a_cur_time,
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
   
   if(a_startup) setupHistFile();

   std::ofstream histFile;
   if(!procID()) {
      histFile.open("history.txt", std::ios_base::app);
      if(histFile.is_open()) {
         histFile << a_cur_step << " ";
         histFile << std::setprecision(5) << std::scientific << a_cur_time << " ";
         //histFile << a_cur_time << " ";
         histFile << "\n";
         histFile.close();
      }
   }  

}

void System::setupHistFile() 
{
   if(!procID()) cout << "setting up the history file ..." << endl << endl;
 
   std::ofstream histFile("history.txt",std::ofstream::out);

   if(histFile.is_open()) {
      histFile << "#0 step number \n";
      histFile << "#1 time [code units] \n";
      histFile.close();
   }
   else {
      cout << "Problem creating/opening history.txt file" << endl;
   }

}

void System::parseParameters( ParmParse&  a_ppsys )
{

   //CH_assert(a_ppsys.contains("advance_method"));
   std::string advance_method_string = "DSMC";
   a_ppsys.query("advance_method", advance_method_string);
   if(advance_method_string=="DSMC") {
      m_advance_method = DSMC;
   }
   else if(advance_method_string=="PICMC_EXPLICIT") {
      m_advance_method = PICMC_EXPLICIT;
   }
   else if(advance_method_string=="PICMC_SEMI_IMPLICIT") {
      m_advance_method = PICMC_SEMI_IMPLICIT;
      a_ppsys.get("iter_max",m_iter_max);
      a_ppsys.query("abs_tol", m_atol);
      a_ppsys.query("rel_tol", m_rtol);
   }
   else if(advance_method_string=="PICMC_FULLY_IMPLICIT") {
      m_advance_method = PICMC_FULLY_IMPLICIT;
      a_ppsys.get("iter_max",m_iter_max);
      a_ppsys.query("abs_tol", m_atol);
      a_ppsys.query("rel_tol", m_rtol);
      a_ppsys.query("theta_parameter",m_theta);
      //CH_assert(m_theta>0.0 && m_theta<=1.0);
      CH_assert(m_theta>=0.5 && m_theta<=1.0);
   }
   else {
      if(!procID()) cout << "advance method = " << advance_method_string << " is not defined " << endl;
      exit(EXIT_FAILURE);
   }

   a_ppsys.query("write_species_charge_density", m_writeSpeciesChargeDensity);
   a_ppsys.query("write_species_current_density", m_writeSpeciesCurrentDensity);
 
   if(!procID()) {
      cout << "advance method  = " << advance_method_string << endl;
      if(m_advance_method==PICMC_SEMI_IMPLICIT) {
         cout << "iter_max = " << m_iter_max << endl;
         cout << "abs_tol  = " << m_atol << endl;
         cout << "rel_tol  = " << m_rtol << endl;
      }
      if(m_advance_method==PICMC_FULLY_IMPLICIT) {
         cout << "iter_max = " << m_iter_max << endl;
         cout << "abs_tol  = " << m_atol << endl;
         cout << "rel_tol  = " << m_rtol << endl;
         cout << "theta_parameter = " << m_theta << endl; 
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

   if(m_advance_method==PICMC_EXPLICIT) {  
   
      Real cnormHalfDt = 0.5*a_dt*m_units->CvacNorm();

      //
      // complete advance of E from t_{n} to t_{n+1/2}
      //

      if (!m_electromagneticFields.isNull() && m_electromagneticFields->advanceE()) {
         if (a_cur_time==0.0) { // initial advance of electric field by 1/2 time step
            m_electromagneticFields->setCurlB();
            //m_electromagneticFields->setCurrentDensity(m_pic_species_ptr_vect);
            m_electromagneticFields->advanceElectricField(cnormHalfDt);
         }
         else {
            m_electromagneticFields->advanceElectricField_2ndHalf();
         }
         m_electromagneticFields->updateOldElectricField();
      }
  
      //
      // complete advance of xp from t_{n} to t_{n+1/2}
      //

      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->motion()) {
            if (a_cur_time==0.0) { // initial advance of particle positions by 1/2 time step
               this_picSpecies->advancePositions(cnormHalfDt);
            }
            else {
               this_picSpecies->advancePositions_2ndHalf();
            }
            this_picSpecies->updateOldParticlePositions();
         }
      }

   }
   
   if(m_advance_method==PICMC_SEMI_IMPLICIT) {  

      //
      // complete advance of B from t_{n} to t_{n+1/2}
      //

      if (!m_electromagneticFields.isNull() && m_electromagneticFields->advanceB()) {
         if (a_cur_time==0.0) { // initial advance of magnetic field by 1/2 time step
            Real cnormHalfDt = 0.5*a_dt*m_units->CvacNorm();
            m_electromagneticFields->setCurlE();
            m_electromagneticFields->advanceMagneticField(cnormHalfDt);
         }
         else {
            m_electromagneticFields->advanceMagneticField_2ndHalf();
         }
         m_electromagneticFields->updateOldMagneticField();
         m_electromagneticFields->setCurlB(); // update curlB here
      }
  
   }

}

void System::advance( Real&  a_cur_time,
                      Real&  a_dt,
                      int&   a_step_number )
{  
   CH_TIME("System::advance()");

   switch(m_advance_method) {
      case DSMC : 
          advance_DSMC( a_dt );
          break;
      case PICMC_EXPLICIT : 
          advance_PICMC_EXPLICIT( a_dt );
          break;
      case PICMC_SEMI_IMPLICIT : 
          advance_PICMC_SEMI_IMPLICIT( a_dt );
          break;
      case PICMC_FULLY_IMPLICIT : 
          advance_PICMC_FULLY_IMPLICIT( a_dt );
          break;
      default :
          MayDay::Error("Invalid advance method");
   }
   
   // update current time and step number
   a_cur_time = a_cur_time + a_dt;
   a_step_number = a_step_number + 1;

}

void System::advance_DSMC( Real&  a_dt )
{  
   CH_TIME("System::advance_DSMC()");
   
   //
   // explicit advance of particle positions 
   // + velocity scatter
   // + special operators
   //
    
   Real cnormDt = a_dt*m_units->CvacNorm();

   // advance particle positions from t_{n} to t_{n+1} using vp_{n}
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->motion()) this_picSpecies->advancePositions(cnormDt);
   }
   
   // scatter the particles: vp_{n+1} ==> vp'_{n+1}
   if(m_use_scattering) scatterParticles( a_dt );
   
   // apply special operators
   if(m_use_specialOps) {
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) { // loop over pic species
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         m_specialOps->applyOp(*this_picSpecies,*m_mesh,cnormDt);
      }
      m_specialOps->updateOp(*m_mesh,*m_units,cnormDt);
   }

}

void System::advance_PICMC_EXPLICIT( Real&  a_dt )
{  
   CH_TIME("System::advance_PICMC_EXPLICIT()");
   
   //
   // explicit leap-frog time advance of particles and fields
   // B and vp are defined a whole time steps while E and xp are at half
   //
    
   Real cnormDt = a_dt*m_units->CvacNorm();
   Real cnormHalfDt = 0.5*cnormDt;
      
   if (!m_electromagneticFields.isNull()) {
   
      //
      // Step 1: advance B from t_{n} to t_{n+1/2} using E_{n+1/2}
      //
      if (m_electromagneticFields->advanceB()) {
         m_electromagneticFields->setCurlE();
         m_electromagneticFields->advanceMagneticField(cnormHalfDt);
      }

      //
      // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2} and advance vp from t_{n} to t_{n+1}
      //
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->forces()) {
            const LevelData<EdgeDataBox>& Efield = m_electromagneticFields->getElectricField();
            const LevelData<FluxBox>& Bfield = m_electromagneticFields->getMagneticField();
            this_picSpecies->interpolateElectricFieldToParticles( Efield );
            this_picSpecies->interpolateMagneticFieldToParticles( Bfield );
            if(SpaceDim<3) {
               const LevelData<NodeFArrayBox>& Efield_virt = m_electromagneticFields->getVirtualElectricField();
               const LevelData<FArrayBox>& Bfield_virt = m_electromagneticFields->getVirtualMagneticField();
               this_picSpecies->interpolateElectricFieldToParticles( Efield_virt );
               this_picSpecies->interpolateMagneticFieldToParticles( Bfield_virt );
            }
            this_picSpecies->advanceVelocities( cnormDt );
         }
      }

   }
   
   // scatter the particles: vp_{n+1} ==> vp'_{n+1}
   if(m_use_scattering) scatterParticles( a_dt );
   

   // advance particle positions from t_{n+1/2} to t_{n+1} using vp_{n+1}
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->motion()) this_picSpecies->advancePositions(cnormHalfDt);
   }
   
   if (!m_electromagneticFields.isNull()) {
   
      // complete advance of B from t_{n+1/2} to t_{n+1} and compute the curl
      if (m_electromagneticFields->advanceB()) {
         m_electromagneticFields->advanceMagneticField_2ndHalf();
         m_electromagneticFields->setCurlB();
      }

      // compute current density at t_{n+1} and advance E from t_{n+1/2} to t_{n+1} 
      // using B_{n+1} and J_{n+1}
      if (m_electromagneticFields->advanceE()) {
         for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
            PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
            int this_charge = this_picSpecies->charge();
            if(this_charge != 0) this_picSpecies->setCurrentDensity();
         }
         m_electromagneticFields->setCurrentDensity(m_pic_species_ptr_vect);
         m_electromagneticFields->advanceElectricField(cnormHalfDt);
      }

   }
   
   // apply special operators
   if(m_use_specialOps) {
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) { // loop over pic species
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         m_specialOps->applyOp(*this_picSpecies,*m_mesh,cnormDt);
      }
      m_specialOps->updateOp(*m_mesh,*m_units,cnormDt);
   }

}

void System::advance_PICMC_SEMI_IMPLICIT( Real&  a_dt )
{  
   CH_TIME("System::advance_PICMC_SEMI_IMPLICIT()");

   //
   // semi-implicit semi-energy-conservative scheme by Chen 2020
   // B is defined at half time steps while all other quantities
   // (E, xp, and vp) are defined at whole time steps
   //
    
   Real cnormDt = a_dt*m_units->CvacNorm();
   Real cnormHalfDt = 0.5*cnormDt;
  
   int iter = 0;
   Real norm_efield = 0, norm_efield0 = 0;
   while( 1 ) {

      //
      // Step 1: advance particle positions from t_{n} to t_{n+1/2} using vp_{n+1/2}
      //
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->motion()) this_picSpecies->advancePositions(cnormHalfDt);
      }
      
      //
      // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2} and advance vp from t_{n} to t_{n+1/2}
      //
      if (!m_electromagneticFields.isNull()) {
     
         for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
            PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
            if(this_picSpecies->forces()) {
               const LevelData<EdgeDataBox>& Efield = m_electromagneticFields->getElectricField();
               const LevelData<FluxBox>& Bfield = m_electromagneticFields->getMagneticField();
               this_picSpecies->interpolateElectricFieldToParticles( Efield );
               this_picSpecies->interpolateMagneticFieldToParticles( Bfield );
               if(SpaceDim<3) {
                  const LevelData<NodeFArrayBox>& Efield_virt = m_electromagneticFields->getVirtualElectricField();
                  const LevelData<FArrayBox>& Bfield_virt = m_electromagneticFields->getVirtualMagneticField();
                  this_picSpecies->interpolateElectricFieldToParticles( Efield_virt );
                  this_picSpecies->interpolateMagneticFieldToParticles( Bfield_virt );
               }
               this_picSpecies->advanceVelocities( cnormDt, true );
            }
         }
      
         // 
         // Step 3: compute current density at t_{n+1/2} and advance E from t_n to t_{n+1/2} 
         //
         if (m_electromagneticFields->advanceE()) {
            for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
               PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
               int this_charge = this_picSpecies->charge();
               if(this_charge != 0) this_picSpecies->setCurrentDensity();
            }
            m_electromagneticFields->setCurrentDensity(m_pic_species_ptr_vect);
            //m_electromagneticFields->setCurlB(); // done in preTimeStep()
            m_electromagneticFields->saveElectricField();
            m_electromagneticFields->advanceElectricField(cnormHalfDt);
            norm_efield = m_electromagneticFields->diffElectricField();
            if (iter == 0) norm_efield0 = norm_efield;

            if (!procID()) {
              printf("  iter = %3d,", iter);
              printf(" norm_efield = %1.4e (abs.), %1.4e (rel.)\n",
                     norm_efield, norm_efield/norm_efield0 );
            }
            if (norm_efield < m_atol) {
              if (!procID()) {
                printf("  exiting: satisfied absolute tolerance (%1.3e).\n", 
                       m_atol);
              }
              break;
            }
            if (norm_efield/norm_efield0 < m_rtol) {
              if (!procID()) {
                printf("  exiting: satisfied relative tolerance (%1.3e).\n", 
                       m_rtol);
              }
              break;
            }

         }

         if ( iter >= m_iter_max ) {
           if (!procID()) {
             printf("  exiting: iterations exceed max iterations (%d).\n", 
                    m_iter_max);
           }
           break;
         }
 
      }
      else {
         break;
      }
      
      iter = iter + 1;
  
   } // end iteration loop

   //
   // Step 4: 2nd half-advance of xp and vp from t_{n+1/2} to t_{n+1}
   //
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->motion()) this_picSpecies->advancePositions_2ndHalf();
      if(this_picSpecies->forces()) this_picSpecies->advanceVelocities_2ndHalf();
   }

   //
   // Step 5: 2nd half-advance of E from t_{n+1/2} to t_{n+1} and 
   //         1st half-advance of B from t_{n+1/2} to t_{n+1} using E_{n+1}
   //
   if (!m_electromagneticFields.isNull()) {
      if (m_electromagneticFields->advanceE()) m_electromagneticFields->advanceElectricField_2ndHalf();
      if (m_electromagneticFields->advanceB()) {
         m_electromagneticFields->setCurlE();
         m_electromagneticFields->advanceMagneticField(cnormHalfDt);
      }
   }

   //   
   // Step 6: scatter the particles: vp_{n+1} ==> vp'_{n+1}
   //
   if(m_use_scattering) scatterParticles( a_dt );

}

void System::advance_PICMC_FULLY_IMPLICIT( Real&  a_dt )
{  
   CH_TIME("System::advance_PICMC_FULLY_IMPLICIT()");

   //
   // fully-implicit energy-conservative scheme by Markidis and
   // Lapenta 2011. All quantities (E, B, xp, and vp) are defined 
   // at whole time steps
   //
    
   Real cnormDt = a_dt*m_units->CvacNorm();
   Real cnormHalfDt = 0.5*cnormDt;
   Real cnormThetaDt = m_theta*cnormDt;
  
   int iter = 0;
   Real  norm_efield = 0, norm_efield0 = 0, 
         norm_bfield = 0, norm_bfield0 = 0;
   while(1) {
 
      //
      // Step 1: advance particle positions from t_{n} to t_{n+1/2} using vp_{n+1/2}
      //
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->motion()) this_picSpecies->advancePositions(cnormHalfDt);
      }
      
      if (!m_electromagneticFields.isNull()) {
      
         //
         // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2} and advance 
         // vp from t_{n} to t_{n+1/2}
         //
         for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
            PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
            if(this_picSpecies->forces()) {
               const LevelData<EdgeDataBox>& Efield = m_electromagneticFields->getElectricField();
               const LevelData<FluxBox>& Bfield = m_electromagneticFields->getMagneticField();
               this_picSpecies->interpolateElectricFieldToParticles( Efield );
               this_picSpecies->interpolateMagneticFieldToParticles( Bfield );
               if(SpaceDim<3) {
                  const LevelData<NodeFArrayBox>& Efield_virt = m_electromagneticFields->getVirtualElectricField();
                  const LevelData<FArrayBox>& Bfield_virt = m_electromagneticFields->getVirtualMagneticField();
                  this_picSpecies->interpolateElectricFieldToParticles( Efield_virt );
                  this_picSpecies->interpolateMagneticFieldToParticles( Bfield_virt );
               }
               this_picSpecies->advanceVelocities( cnormDt, true );
            }
         }
      
         // 
         // Step 3: compute current density at t_{n+1/2} and advance B and E 
         // from t_n to t_{n+1/2} 
         //
         if (m_electromagneticFields->advance()) {
            
            // advance B first
            if (m_electromagneticFields->advanceB()) {
               m_electromagneticFields->setCurlE();
               m_electromagneticFields->saveMagneticField();
               m_electromagneticFields->advanceMagneticField(cnormThetaDt);
               norm_bfield = m_electromagneticFields->diffMagneticField();
            }

            // then advance E
            if (m_electromagneticFields->advanceE()) {
               for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
                  PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
                  int this_charge = this_picSpecies->charge();
                  if(this_charge != 0) this_picSpecies->setCurrentDensity();
               }
               m_electromagneticFields->setCurrentDensity(m_pic_species_ptr_vect);
               m_electromagneticFields->setCurlB();
               m_electromagneticFields->saveElectricField();
               m_electromagneticFields->advanceElectricField(cnormThetaDt);
               norm_efield = m_electromagneticFields->diffElectricField();
            }

            if (iter == 0) norm_bfield0 = norm_bfield;
            if (iter == 0) norm_efield0 = norm_efield;

            if (!procID()) {
              printf("  iter = %3d,", iter);
              printf(" norm_efield = %1.4e (abs.), %1.4e (rel.),",
                     norm_efield, norm_efield/norm_efield0 );
              printf(" norm_bfield = %1.4e (abs.), %1.4e (rel.)\n",
                     norm_bfield, norm_bfield/norm_bfield0 );
            }
            if ((norm_efield < m_atol) && (norm_bfield < m_atol)) {
              if (!procID()) {
                printf("  exiting: satisfied absolute tolerance (%1.3e).\n", 
                       m_atol);
              }
              break;
            }
            if (     (norm_efield/norm_efield0 < m_rtol)
                 &&  (norm_bfield/norm_bfield0 < m_rtol) ) {
              if (!procID()) {
                printf("  exiting: satisfied relative tolerance (%1.3e).\n", 
                       m_rtol);
              }
              break;
            }

         }
 
         if ( iter >= m_iter_max ) {
           if (!procID()) {
             printf("  exiting: iterations exceed max iterations (%d).\n", 
                    m_iter_max);
           }
           break;
         }

      }
      else {
         break;
      }
      
      iter = iter + 1;
  
   } // end iteration loop

   //
   // Step 4: 2nd half-advance of all quantities (E, B, xp, and vp)
   // from t_{n+1/2} to t_{n+1}
   //
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->motion()) this_picSpecies->advancePositions_2ndHalf();
      if(this_picSpecies->forces()) this_picSpecies->advanceVelocities_2ndHalf();
   }

   if (!m_electromagneticFields.isNull()) {
      if(m_electromagneticFields->advanceE()) m_electromagneticFields->advanceElectricField_2ndHalf(m_theta);
      if(m_electromagneticFields->advanceB()) m_electromagneticFields->advanceMagneticField_2ndHalf(m_theta);
   }

   //   
   // Step 6: scatter the particles: vp_{n+1} ==> vp'_{n+1}
   //
   if(m_use_scattering) scatterParticles( a_dt );

}

void System::scatterParticles( const Real&  a_dt )
{  
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
      this_picSpecies1->setEnergyDensity();
      if(this_species1==this_species2) { // self-species scattering
         this_scattering->applySelfScattering( *this_picSpecies1, *m_mesh, dt_sec );
      }
      else { // inter-species scattering
         PicSpeciesPtr this_picSpecies2(m_pic_species_ptr_vect[this_species2]);
         this_picSpecies2->setEnergyDensity();
         this_scattering->applyInterScattering( *this_picSpecies1, *this_picSpecies2, *m_mesh, dt_sec );
      }

   }

}   

void System::postTimeStep( const Real&  a_cur_time,
                           const Real&  a_dt,
                           const int&   a_step_number )
{  
   CH_TIME("System::postTimeStep()");
   
   if(m_advance_method==DSMC) {  
   
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->motion()) this_picSpecies->updateOldParticlePositions();
         if(this_picSpecies->forces()) this_picSpecies->updateOldParticleVelocities();
      }

   }
   
   if(m_advance_method==PICMC_EXPLICIT) {  
   
      if (!m_electromagneticFields.isNull() && m_electromagneticFields->advanceB()) {
         m_electromagneticFields->updateOldMagneticField();
      }
   
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->forces()) this_picSpecies->updateOldParticleVelocities();
      }

   }
   
   if(m_advance_method==PICMC_SEMI_IMPLICIT) {  
   
      if (!m_electromagneticFields.isNull() && m_electromagneticFields->advanceE()) {
         m_electromagneticFields->updateOldElectricField();
      }
   
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->motion()) this_picSpecies->updateOldParticlePositions();
         if(this_picSpecies->forces()) this_picSpecies->updateOldParticleVelocities();
      }

   }
   
   if(m_advance_method==PICMC_FULLY_IMPLICIT) {  
   
      if (!m_electromagneticFields.isNull()) {
         if(m_electromagneticFields->advanceE()) m_electromagneticFields->updateOldElectricField();
         if(m_electromagneticFields->advanceB()) m_electromagneticFields->updateOldMagneticField();
      }
   
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->motion()) this_picSpecies->updateOldParticlePositions();
         if(this_picSpecies->forces()) this_picSpecies->updateOldParticleVelocities();
      }

   }


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
   if(m_advance_method==PICMC_EXPLICIT) a_adapt_dt = false;
   if(m_advance_method==PICMC_SEMI_IMPLICIT) a_adapt_dt = false;
}


#include "NamespaceFooter.H"
