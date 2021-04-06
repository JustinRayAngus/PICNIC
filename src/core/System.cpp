#include "System.H"
#include "CH_Timer.H"

#include "BoxIterator.H"
#include "DomainGrid.H"

#include "dataFileIO.H"

#include "SpecialOperatorFactory.H"
#include "ScatteringFactory.H"

#include "NamespaceHeader.H"


System::System( ParmParse& a_pp )
   :
     m_mesh(NULL),
     m_units(NULL),
     m_dataFile(NULL),
     m_meshInterp(NULL),
     m_verbosity(0)
{
   ParmParse ppsys("system");
   parseParameters(ppsys);

   m_units = new CodeUnits( ppsys );
   if(!procID()) m_units->printParameters( procID() );
 
   createProblemDomain();          // create the problem domain
  
   DisjointBoxLayout grids;
   getDisjointBoxLayout( grids );  // define the disjointBoxLayout
   
   // initialize the coordinates and grid
   //
   ParmParse ppgrid( "grid" );
   m_mesh = new DomainGrid( ppgrid, m_num_ghosts, m_domain, grids ); 

   m_dataFile = new dataFileIO( a_pp, *m_mesh );
   //m_dataFile = RefCountedPtr<dataFileIO>(new dataFileIO( a_pp, *mesh));
     
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
         this_picSpecies->initialize(); // set initial particle positions and velocities
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
   }

}

void System::createProblemDomain()
{
   CH_TIME("System::createProblemDomain()");

   ParmParse ppgrid( "grid" );
   ppgrid.get( "geometry", m_geom_type );
   ppgrid.query( "num_ghosts", m_num_ghosts );
   int DIM = SpaceDim;

   if ( m_geom_type == "cylindrical" || m_geom_type == "cartesian") {
     
      // Set the grid size
      m_num_cells.resize( DIM );
      for (int i=0; i<DIM; ++i) m_num_cells[i] = 0;
      ppgrid.getarr( "num_cells", m_num_cells, 0, DIM );
      for (int i=0; i<DIM; ++i) CH_assert( m_num_cells[i]>0 );

      // Determine which spatial directions are periodic
      m_is_periodic.resize(DIM);
      vector<int> isPeriodic( DIM ); // why should I have to do this?
      ppgrid.getarr( "is_periodic", isPeriodic, 0, DIM );
      for (int dim=0; dim<SpaceDim; dim++)  {
         m_is_periodic[dim] = (isPeriodic[dim] == 1);
      }

      // Get the domain box decomposition parameters
      if (ppgrid.contains("config_decomp")) {
         m_config_decomp.resize( DIM );
         for (int i=0; i<DIM; ++i) m_config_decomp[i] = 0;
         ppgrid.getarr( "config_decomp", m_config_decomp, 0, DIM );
         for (int i=0; i<DIM; ++i) CH_assert( m_config_decomp[i]>0 );
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
   //
   IntVect boxSize;
   boxSize[0] = domain_box.size(0)/m_config_decomp[0];
   for (int dir=1; dir<SpaceDim; ++dir) {
      boxSize[dir] = domain_box.size(dir)/m_config_decomp[dir];
      CH_assert(boxSize[dir]==boxSize[0]);
      //if(!procID()) cout << "JRA: boxSize[dir] = " << boxSize[dir] << endl;
   }
   
   // Chop up the configuration space domain box over the number of processors specified
   // for this block (single block here).  At this point, we insist that the box 
   // decomposes uniformly, or an error is thrown.
   //
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
   // create/set the meshInterp object pointer
   //m_meshInterp = new MeshInterp( domain.domainBox(),
   //                               meshSpacing,
   //                               meshOrigin );
   m_meshInterp = static_cast<MeshInterp*> (new MeshInterp( domain.domainBox(),
                                                            meshSpacing,
                                                            meshOrigin  ));
   //m_meshInterp = static_cast<RefCountedPtr<MeshInterp>> (new MeshInterp( domain.domainBox(),
   //                                                                       meshSpacing,
   //                                                                       meshOrigin  ));

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
      m_electromagneticFields = RefCountedPtr<ElectroMagneticFields>(new ElectroMagneticFields(ppflds,*m_mesh,verbose)); 
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
         m_specialOps->updateOp(*m_mesh,0.0); // should make an initializeOp() function
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
   
   //m_dataFile->writeFieldDataFile();
   
}

void System::writeHistFile( const int     a_cur_step,
                            const double  a_cur_time,
                            const bool    a_startup )
{
   CH_TIME("System::writeHistFile()");
   if(!procID()) cout << "System::writeHistFile() a_startup " << a_startup << endl;
   
   if(a_startup) setupHistFile();

}

void System::setupHistFile() 
{
   if(!procID()) cout << "System::setupHistFile() " << endl;
 
   // need to create the hdf5 file and set it up for being appended to 

}

void System::parseParameters( ParmParse&  a_ppsys )
{
   
   a_ppsys.query("write_species_charge_density", m_writeSpeciesChargeDensity);
   a_ppsys.query("write_species_current_density", m_writeSpeciesCurrentDensity);
 
   if(!procID()) {
      cout << "write species charge density  = " << m_writeSpeciesChargeDensity << endl;
      cout << "write species current density = " << m_writeSpeciesCurrentDensity << endl << endl;
   }
 
}

void System::advance( Real&  a_cur_time,
                      Real&  a_dt,
                      int&   a_step_number )
{  
   CH_TIME("System::advance()");
   

   // advance the electromagnetic fields
   if (!m_electromagneticFields.isNull()) {
      if(m_electromagneticFields->advance()) {
         m_electromagneticFields->setCurlB();
         m_electromagneticFields->advanceElectricField(a_dt);
         m_electromagneticFields->setCurlE();
         m_electromagneticFields->advanceMagneticField(a_dt);
         m_electromagneticFields->updateOldFieldValues();
      }
   }
   
   // advance particle positions
   for (int s=0; s<m_pic_species_ptr_vect.size(); s++) { // loop over pic species
      PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
      if(this_picSpecies->motion()) this_picSpecies->advancePositions(a_dt);
   }
   
   // scatter the particles
   if(m_use_scattering) {
      
      // prepare all species to be scattered for scattering
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) {
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         if(this_picSpecies->scatter()) {
            this_picSpecies->binTheParticles();
            this_picSpecies->setNumberDensity();
            this_picSpecies->setEnergyDensity(); // may want to do this after or before each scatter
         }
      }

      // Should we shuffle the m_scattering_ptr_vect each time ?
      for (int sct=0; sct<m_scattering_ptr_vect.size(); sct++) { // loop over scattering objects
         
         ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
         const int this_species1 = this_scattering->species1();
         const int this_species2 = this_scattering->species2();
         
         PicSpeciesPtr this_picSpecies1(m_pic_species_ptr_vect[this_species1]);
         if(this_species1==this_species2) { // self-species scattering
            this_scattering->applySelfScattering( *this_picSpecies1, *m_mesh, a_dt );
         }
         else { // inter-species scattering
            PicSpeciesPtr this_picSpecies2(m_pic_species_ptr_vect[this_species2]);
            this_scattering->applyInterScattering( *this_picSpecies1, *this_picSpecies2, *m_mesh, a_dt );
         }

      }

   }
   
   // apply special operators
   if(m_use_specialOps) {
      for (int s=0; s<m_pic_species_ptr_vect.size(); s++) { // loop over pic species
         PicSpeciesPtr this_picSpecies(m_pic_species_ptr_vect[s]);
         m_specialOps->applyOp(*this_picSpecies,*m_mesh,a_dt);
      }
      m_specialOps->updateOp(*m_mesh,a_dt);
   }

   // update current time and step number
   a_cur_time = a_cur_time + a_dt;
   a_step_number = a_step_number + 1;

}

Real System::fieldsDt( const int a_step_number )
{
   Real fldsDt = DBL_MAX;
   if (!m_electromagneticFields.isNull()) {
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
      this_picSpecies->setStableDt();
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
   if(m_use_scattering) {
      for (int sct=0; sct<m_scattering_ptr_vect.size(); sct++) {
         ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
         scatterDt = min(scatterDt,this_scattering->scatterDt());
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


#include "NamespaceFooter.H"
