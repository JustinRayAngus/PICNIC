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
     m_picSpecies(NULL),
     m_verbosity(0)
{
   ParmParse ppsys("system");

   m_units = new CodeUnits( ppsys );
   if(!procID()) m_units->printParameters( procID() );
 
   createProblemDomain();          // create the problem domain
  
   DisjointBoxLayout grids;
   getDisjointBoxLayout( grids );  // define the disjointBoxLayout
   
   // initialize the coordinates and grid
   //
   ParmParse ppgrid( "grid" );
   m_mesh = new DomainGrid( ppgrid, m_num_ghosts, m_domain, grids ); 
   //DomainGrid domainGrid( ppgrid, m_num_ghosts, m_domain, grids ); 
   //DomainGrid* mesh = DomainGrid::mesh;

   m_dataFile = new dataFileIO( a_pp, *m_mesh );
   //m_dataFile = RefCountedPtr<dataFileIO>(new dataFileIO( a_pp, *mesh));
     
   createMeshInterp();

   createState( a_pp );

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
   delete m_picSpecies;
 
}

void System::initialize( const int     a_cur_step,
                         const double  a_cur_time )
{
   CH_TIME("System::initialize()");
   
   //DomainGrid* mesh = DomainGrid::mesh;
   
   // initialize the state_variable data members (ICs) and the operators
   // Will eventually loop over all state_variables objects

   //ParticleData<Particle>& m_Ptest = m_picSpecies->partData(); //ref, so can change
   //const ParticleData<Particle>& m_Ptest = m_picSpecies->partData(); // const ref, so can't change
   if(a_cur_step==0) {
      m_picSpecies->initialize(); // set initial particle positions and velocities
      //m_picSpecies->initialize( m_mesh ); // set initial particle positions and velocities
   }
   else { // restart
      if(!procID()) cout << "System::initialize() called at step = " << a_cur_step << endl;
      if(!procID()) cout << "restarting not currently an option !!! " << endl;
      exit(EXIT_FAILURE);
   }

   // check for scattering and initialize scattering models
   //
   if(m_picSpecies->scatter()) {

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
      //
      m_num_cells.resize( DIM );
      for (int i=0; i<DIM; ++i) m_num_cells[i] = 0;
      ppgrid.getarr( "num_cells", m_num_cells, 0, DIM );
      for (int i=0; i<DIM; ++i) CH_assert( m_num_cells[i]>0 );

      // Determine which spatial directions are periodic
      //
      m_is_periodic.resize(DIM);
      vector<int> isPeriodic( DIM ); // why should I have to do this?
      ppgrid.getarr( "is_periodic", isPeriodic, 0, DIM );
      for (int dim=0; dim<SpaceDim; dim++)  {
         m_is_periodic[dim] = (isPeriodic[dim] == 1);
      }

      // Get the domain box decomposition parameters
      //
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
   //
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
   CH_TIME("System::createState()");
   
   // create kinetic species
   //
   //PICspeciesPtrVect pic_species;
   createPICspecies();
   
   createScattering();

   // create special operators
   //
   createSpecialOperators();

   /*

   FieldPtrVect fields;
   createFields( fields );

   */

}


void System::createPICspecies()
{
   //
   // Create the vector of species model (pointers), and when a kinetic
   // species is created, add it to the kinetic species vector
   // number species.
   //
   
   if(!procID()) {
      cout << "Adding PIC species..." << endl;
   }
   //DomainGrid* mesh = DomainGrid::mesh;

   bool more_vars(true);
   string name0;
   int species = 0;
   while(more_vars) { // look for pic species...
 
      species = species + 1;

      stringstream s;
      s << "pic_species." << species; 

      //ParmParse ppspc( "pic_species.1" );
      ParmParse ppspc( s.str().c_str() );
     
      string name;
      if(ppspc.contains("name")) {
         ppspc.get("name",name);
         name0 = name;
      } 
      else {
         more_vars = false;
      }
   
      if(more_vars) {
         //m_picSpecies = new PicSpecies( ppspc, name,  *mesh );
         //m_picSpecies = new PicSpecies( ppspc, name, m_meshInterp, *mesh );
         m_picSpecies = new PicSpecies( ppspc, name, *m_meshInterp, *m_mesh );
         //m_picSpecies = RefCountedPtr<PicSpecies>(new PicSpecies( ppspc, name, *m_meshInterp, *m_mesh ));
         //m_picSpecies = RefCountedPtr<PicSpecies>(new PicSpecies( ppspc, name, *m_mesh ));
         //m_picSpecies = RefCountedPtr<PicSpecies>(new PicSpecies( ppspc, name, *m_mesh ));
      }

   }

   if(!procID()) {
      cout << "Done adding PIC species" << endl << endl;
   }

}


void System::createScattering()
{
   if(m_picSpecies->scatter()) {  

      ScatteringFactory  scatteringFactory;
   
      if(!procID()) {
         cout << "Adding scattering operators... " << endl;
      }
      
      int species_num = 1;
      stringstream s;
      s << "scattering." << m_picSpecies->name(); 
      
      ParmParse pp_scatter( s.str().c_str() );
      if(pp_scatter.contains("model")) {
         m_use_scattering = true;
         m_scattering = scatteringFactory.create( pp_scatter, *m_picSpecies, 1 );
      } 
      
      if(!procID()) {
         cout << "Done adding scattering models" << endl << endl;
      }

   }

}

void System::createSpecialOperators()
{
   //
   // Create the vector of special operator (pointers)
   //
   
   SpecialOperatorFactory  specialOpFactory;
   
   if(!procID()) {
      cout << "Adding special operators..." << endl;
   }
   //DomainGrid* mesh = DomainGrid::mesh;

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

   // eventually will loop over species/vars here
   //
   //ParticleData<Particle>& Ptest = m_picSpecies->partData(); //ref, so can change
   //const ParticleData<Particle>& Pdata = m_picSpecies->partData(); // const ref, so can't change
   const ParticleData<JustinsParticle>& Pdata = m_picSpecies->partData(); // const ref, so can't change
   //const LevelData<BinFab<JustinsParticle>>& Pdata_binfab = m_picSpecies->partData_binfab(); // const ref, so can't change
   //m_picSpecies->setNumberDensity(); 
   const bool setDensity = true;
   const LevelData<FArrayBox>& density = m_picSpecies->getNumberDensity(setDensity);
   /*
   const DisjointBoxLayout& grids = density.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      const Box gridBox = grids.get(dit);
      BoxIterator gbit(gridBox); 
      for(gbit.begin(); gbit.ok(); ++gbit) {
         if(!procID()) {
            const IntVect ig = gbit();
            const Real thisDen = density[dit].get(ig,0); 
            cout << "JRA: ig = " << ig << endl;
            cout << "JRA: thisDen = " << thisDen << endl;
         } 
      }
   }
   */
   /*
   const DisjointBoxLayout& grids = density.disjointBoxLayout();
   LevelData<FArrayBox> rho;
   rho.define(grids,density.nComp(),density.ghostVect());
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      rho[dit].copy(density[dit]);
   }
   rho.exchange();
   */
   const bool setMomentum = true;
   const LevelData<FArrayBox>& momentum = m_picSpecies->getMomentumDensity(setMomentum);
   
   const bool setEnergy = true;
   const LevelData<FArrayBox>& energy = m_picSpecies->getEnergyDensity(setEnergy);

   //m_picSpecies->inspectBinFab();

   m_dataFile->writeParticleDataFile( Pdata, density, momentum, energy, m_picSpecies->mass(),
                                      a_cur_step, a_cur_time );
   
   //m_dataFile->writeBinFabDataFile( Pdata_binfab, a_cur_step, a_cur_time );

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
}

void System::advance( Real&  a_cur_time,
                      Real&  a_dt,
                      int&   a_step_number )
{  
   CH_TIME("System::advance()");
   
   //m_integrator->setTimeStepSize( a_dt ); // pass time step to integrator
  
   // advance particle positions
   //
   if(m_picSpecies->motion()) {
      //ParticleData<JustinsParticle>& Pdata = m_picSpecies->partData();
      m_picSpecies->advancePositions(a_dt);
   }
   
   // scatter the particles
   //
   if(m_use_scattering) {
      //m_picSpecies->testParticleShuffling(a_dt);
      if(m_picSpecies->scatter()) {
         m_picSpecies->binTheParticles(); // bin particles up for scattering call below
         const bool setMoments = true;
         const LevelData<FArrayBox>& numberDensity = m_picSpecies->getNumberDensity(setMoments);
         const LevelData<FArrayBox>& energyDensity = m_picSpecies->getEnergyDensity(setMoments);
         m_scattering->applySelfScattering( *m_picSpecies, *m_mesh, numberDensity, energyDensity, a_dt );
      }
   }
   
   // apply special operators
   //
   if(m_use_specialOps) {
      m_specialOps->applyOp(*m_picSpecies,*m_mesh,a_dt);
      m_specialOps->updateOp(*m_mesh,a_dt);
   }

   a_cur_time = a_cur_time + a_dt;
   a_step_number = a_step_number + 1;

}


Real System::stableDt( const int a_step_number )
{
   m_picSpecies->setStableDt();
   Real stableDt = m_picSpecies->stableDt();
   if(!procID()) cout << "stable particle time step = " << stableDt << endl;
   return stableDt;
}


Real System::scatterDt( const int a_step_number )
{
   Real scatterDt = DBL_MAX;
   if(m_use_scattering) scatterDt = m_scattering->scatterDt();
   if(!procID()) cout << "mean free scattering time = " << scatterDt << endl;
   return scatterDt;
}


#include "NamespaceFooter.H"
