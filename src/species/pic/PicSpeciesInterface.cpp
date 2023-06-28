
#include "PicSpeciesInterface.H"
#include "ParmParse.H"
#include "CH_HDF5.H"

#include "FieldsF_F.H"

#include "SpaceUtils.H"
#include "MathUtils.H"
#include "PicnicConstants.H"

#include "NamespaceHeader.H"

PicSpeciesInterface::PicSpeciesInterface( const CodeUnits&   a_units,
                                          const DomainGrid&  a_mesh )
   : m_verbosity(true),
     m_writeSpeciesChargeDensity(false),
     m_writeSpeciesCurrentDensity(false),
     m_writeSpeciesNppc(false),
     m_part_order_swap(false),
     m_iter_min_two(false),
     m_iter_max_particles(0),
     m_rtol_particles(1.0e-12),
     m_newton_num_guess(1),
     m_freeze_particles_jacobian(false),
     m_quasi_freeze_particles_jacobian(false),
     m_use_mass_matrices(false),
     m_mod_init_advance(true),
     m_num_species(0),
     m_mesh(a_mesh),
     m_courant_dt(DBL_MAX)
{

   // parse parameters from input file
   ParmParse pp("pic_species");
   pp.query("write_species_charge_density", m_writeSpeciesChargeDensity);
   pp.query("write_species_current_density", m_writeSpeciesCurrentDensity);
   pp.query("write_species_nppc", m_writeSpeciesNppc);
   //
   pp.query("part_order_swap",m_part_order_swap);
   pp.query("iter_min_two",m_iter_min_two);
   pp.query("iter_max_particles",m_iter_max_particles);
   pp.query("rtol_particles",m_rtol_particles);
   pp.query("newton_num_guess",m_newton_num_guess);
   //
   pp.query("freeze_particles_jacobian",m_freeze_particles_jacobian);
   pp.query("quasi_freeze_particles_jacobian",m_quasi_freeze_particles_jacobian);
   pp.query("use_mass_matrices",m_use_mass_matrices);
   pp.query("mod_init_advance",m_mod_init_advance);
      
   if(m_use_mass_matrices) m_quasi_freeze_particles_jacobian = true;
   if(m_quasi_freeze_particles_jacobian) m_freeze_particles_jacobian=false;
   
   createAllPicSpecies( a_mesh );
      
   m_num_species = m_pic_species_ptr_vect.size();
   
   // define the LevelData members for moments
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   m_chargeDensity.define(grids,1,ghostVect);
   m_chargeDensity_faces.define(grids,1,ghostVect);
   m_chargeDensity_nodes.define(grids,1,ghostVect);
   
   m_currentDensity.define(grids,1,ghostVect);
   m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect);
   
   m_DebyeLength.define(grids,1,ghostVect);
   
}

void
PicSpeciesInterface::createAllPicSpecies( const DomainGrid&  a_mesh )
{

   if(!procID()) cout << "Creating PIC species..." << endl << endl;

   bool more_vars(true);
   int species;
   while(more_vars) { // look for pic species...
 
      species = m_pic_species_ptr_vect.size();
      stringstream s;
      s << "pic_species." << species; 
      ParmParse pp_spc( s.str().c_str() );
     
      string name;
      if(pp_spc.contains("name")) pp_spc.get("name",name);
      else more_vars = false;
   
      if(more_vars) {
         PicSpecies* picSpecies = new PicSpecies( pp_spc, species, name, a_mesh );
         m_pic_species_ptr_vect.push_back(PicSpeciesPtr(picSpecies));
      }

   }

   if(!procID()) {
      cout << "Finished creating " << m_pic_species_ptr_vect.size() << " PIC species" << endl << endl;
   }

}
  
void 
PicSpeciesInterface::initialize( const CodeUnits&    a_units,
                                 const bool          a_implicit_advance,
                                 const Real          a_cur_time,
                                 const std::string&  a_restart_file_name )
{
   if(!procID()) cout << "Initializing pic species..." << endl << endl;
      
   // initialize the pic species
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
      
      PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      species->initialize( a_units, a_cur_time, a_restart_file_name );
      
      if(!a_implicit_advance) {
         bool planar_push = species->planarPush();
         bool boris_inertia = species->borisInertia();
         if(!planar_push && !boris_inertia) {
            if(!procID()) cout << "EXIT_FAILURE: Using non-planar and non-boris push for sp = " << sp << endl;
            if(!procID()) cout << "not setup to work with EXPLICIT advance" << endl;
            exit(EXIT_FAILURE);
         }
      }

      if(a_implicit_advance) {
      
         species->setParticleSolverParams( m_part_order_swap,
                                           m_iter_min_two,
                                           m_iter_max_particles,
                                           m_rtol_particles,
                                           m_newton_num_guess );
      }
  
   }
    
   // additional species-pair initialization
   speciesPairingInit();

   if(a_implicit_advance) {
      bool init_mass_matrices = false;
      if(m_use_mass_matrices) init_mass_matrices = true;
#ifdef MASS_MATRIX_TEST
      CH_assert(!m_use_mass_matrices);
      init_mass_matrices = true;
#endif
      if(init_mass_matrices) initializeMassMatrices( a_restart_file_name );
   }
   else m_use_mass_matrices = false;

   if(!procID()) {
      cout << "Finished initializing " << m_pic_species_ptr_vect.size();
      cout << " pic species" << endl << endl;
      cout << " write species charge density  = " << m_writeSpeciesChargeDensity << endl;
      cout << " write species current density = " << m_writeSpeciesCurrentDensity << endl;
      cout << " write species nppc = " << m_writeSpeciesNppc << endl;
      if (a_implicit_advance) {
         cout << " implicit advance parameters:" << endl;
         cout << "  part_order_swap = " << m_part_order_swap << endl;
         cout << "  iter_min_two = " << m_iter_min_two << endl;
         cout << "  iter_max_particles = " << m_iter_max_particles << endl;
         cout << "  rtol_particles = " << m_rtol_particles << endl;
         cout << "  newton_num_guess = " << m_newton_num_guess << endl;
         cout << "  freeze_particles_jacobian = " << m_freeze_particles_jacobian << endl;
         cout << "  quasi_freeze_particles_jacobian = " << m_quasi_freeze_particles_jacobian << endl;
         cout << "  use_mass_matrices = " << m_use_mass_matrices << endl;
         cout << "  mod_init_advance = " << m_mod_init_advance << endl << endl;
      }
   }

}

void PicSpeciesInterface::initializeMassMatrices( const std::string&  a_restart_file_name )
{
   CH_TIME("PicSpeciesInterface::initializeMassMatrices()");

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit; 

   InterpType interp_type = getInterpForMassMatrix();

   // initialize number of components for CIC/TSC
   int tsc = 0; 
   if(interp_type==TSC) { 
      tsc = 2;
      CH_assert(ghosts>=3);
   }
   for (int dir=0; dir<SpaceDim; ++dir) {
      if(dir==0) {
         m_num_xx_comps[dir] = 3 + tsc;
         m_num_xy_comps[dir] = 4 + tsc;
         m_num_xz_comps[dir] = 4 + tsc;
         //
         m_num_yx_comps[dir] = 4 + tsc;
         m_num_yy_comps[dir] = 3 + tsc;
         m_num_yz_comps[dir] = 3 + tsc;
         //
         m_num_zx_comps[dir] = 4 + tsc;
         m_num_zy_comps[dir] = 3 + tsc;
         m_num_zz_comps[dir] = 3 + tsc;
      }
      if(dir==1) {
         m_num_xx_comps[dir] = 3 + tsc;
         m_num_xy_comps[dir] = 4 + tsc;
         m_num_xz_comps[dir] = 3 + tsc;
         //
         m_num_yx_comps[dir] = 4 + tsc;
         m_num_yy_comps[dir] = 3 + tsc;
         m_num_yz_comps[dir] = 4 + tsc;
         //
         m_num_zx_comps[dir] = 3 + tsc;
         m_num_zy_comps[dir] = 4 + tsc;
         m_num_zz_comps[dir] = 3 + tsc;
      }
      if(dir==2) {
         m_num_xx_comps[dir] = 3 + tsc;
         m_num_xy_comps[dir] = 3 + tsc;
         m_num_xz_comps[dir] = 4 + tsc;
         //
         m_num_yx_comps[dir] = 3 + tsc;
         m_num_yy_comps[dir] = 3 + tsc;
         m_num_yz_comps[dir] = 4 + tsc;
         //
         m_num_zx_comps[dir] = 4 + tsc;
         m_num_zy_comps[dir] = 4 + tsc;
         m_num_zz_comps[dir] = 3 + tsc;
      }
   }
   
   if(interp_type==CIC) {
      CH_assert(ghosts>=2);
   }
   if(interp_type==CC0) {
      m_num_xx_comps[0] = 5; // 2 cell crossings permitted
      CH_assert(ghosts>=2);
      CH_assert(SpaceDim==1);
   }
   if(interp_type==CC1 && SpaceDim==1) {
      CH_assert(ghosts>=2);
      int maxXings = 2;
      if(ghosts==2) maxXings = 1; // use ghost to limit crossings
      m_num_xx_comps[0] = 3 + 2*maxXings;
      m_num_xy_comps[0] = 2 + 2*maxXings;
      m_num_xz_comps[0] = 2 + 2*maxXings;
      m_num_yx_comps[0] = 2 + 2*maxXings;
      m_num_zx_comps[0] = 2 + 2*maxXings;
   }
   if(interp_type==CC1 && SpaceDim==2) {
      CH_assert(ghosts>=3);
      int maxXings = 2;
      if(ghosts==3) maxXings = 1; // use ghost to limit crossings
      for (int dir=0; dir<SpaceDim; ++dir) {
         if(dir==0) {
            m_num_xx_comps[dir] = 3 + 2*maxXings; 
            m_num_xy_comps[dir] = 4 + 2*maxXings; 
            m_num_xz_comps[dir] = 2 + 2*maxXings;
            //
            m_num_yx_comps[dir] = 4 + 2*maxXings; 
            m_num_yy_comps[dir] = 5 + 2*maxXings; 
            m_num_yz_comps[dir] = 3 + 2*maxXings;
            //
            m_num_zx_comps[dir] = 2 + 2*maxXings; 
            m_num_zy_comps[dir] = 3 + 2*maxXings; 
            m_num_zz_comps[dir] = 3;
         } 
         if(dir==1) {
            m_num_xx_comps[dir] = 5 + 2*maxXings; 
            m_num_xy_comps[dir] = 4 + 2*maxXings; 
            m_num_xz_comps[dir] = 3 + 2*maxXings;
            //
            m_num_yx_comps[dir] = 4 + 2*maxXings; 
            m_num_yy_comps[dir] = 3 + 2*maxXings; 
            m_num_yz_comps[dir] = 2 + 2*maxXings;
            //
            m_num_zx_comps[dir] = 3 + 2*maxXings; 
            m_num_zy_comps[dir] = 2 + 2*maxXings; 
            m_num_zz_comps[dir] = 3;
         } 
      }
   }
   
   m_J0.define(grids,1,ghostVect);    // EdgeDataBox with 1 comp
   m_J0_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
   
   m_E0.define(grids,1,ghostVect);    // EdgeDataBox with 1 comp
   m_E0_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
   
   // poor naming convention. m_sigma_x* containers are EdgeBox types.
   // In 2D, m_sigma_xj[0] is sigma_xj for j = x,y,z, but
   //        m_sigma_xx[1] is sigma_yy 
   //        m_sigma_xy[1] is sigma_yx 
   //        m_sigma_xz[1] is sigma_yz 
   m_sigma_xx.define(grids,m_num_xx_comps.product(),ghostVect);
   m_sigma_xy.define(grids,m_num_xy_comps.product(),ghostVect);
   m_sigma_xz.define(grids,m_num_xz_comps.product(),ghostVect);   
#if CH_SPACEDIM==1
   m_sigma_yx.define(grids,m_num_yx_comps.product(),ghostVect);
   m_sigma_yy.define(grids,m_num_yy_comps.product(),ghostVect);
   m_sigma_yz.define(grids,m_num_yz_comps.product(),ghostVect);
#endif
   m_sigma_zx.define(grids,m_num_zx_comps.product(),ghostVect);
   m_sigma_zy.define(grids,m_num_zy_comps.product(),ghostVect);
   m_sigma_zz.define(grids,m_num_zz_comps.product(),ghostVect);
   
   // set all to zero (effects PC for initial step... bug?)
   SpaceUtils::zero( m_J0 );
   SpaceUtils::zero( m_J0_virtual );
   
   SpaceUtils::zero( m_E0 );
   SpaceUtils::zero( m_E0_virtual );
   
   SpaceUtils::zero( m_sigma_xx );
   SpaceUtils::zero( m_sigma_xy );
   SpaceUtils::zero( m_sigma_xz );
#if CH_SPACEDIM==1
   SpaceUtils::zero( m_sigma_yx );
   SpaceUtils::zero( m_sigma_yy );
   SpaceUtils::zero( m_sigma_yz );
#endif
   SpaceUtils::zero( m_sigma_zx );
   SpaceUtils::zero( m_sigma_zy );
   SpaceUtils::zero( m_sigma_zz );

#ifdef MASS_MATRIX_TEST
   m_currentDensity_TEST.define(grids,1,ghostVect);    // EdgeDataBox with 1 comp
   if(SpaceDim<3) m_currentDensity_virtual_TEST.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
   SpaceUtils::zero( m_currentDensity_TEST );
   SpaceUtils::zero( m_currentDensity_virtual_TEST );
#endif

   if(!a_restart_file_name.empty()) {

      if(!procID()) cout << "reading mass matrices from restart file..." << endl;
      HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );
      int group_id;

      //
      // read in the mass matrices (needed for PC on restart)
      //
      
      group_id = handle.setGroup("sigmaxx");
      if(group_id<0) cout << "Warning: sigmaxx not found in restart file" << endl;
      else read(handle, m_sigma_xx, "data", grids);
      
      group_id = handle.setGroup("sigmaxy");
      if(group_id<0) cout << "Warning: sigmaxy not found in restart file" << endl;
      else read(handle, m_sigma_xy, "data", grids);
      
      group_id = handle.setGroup("sigmaxz");
      if(group_id<0) cout << "Warning: sigmaxz not found in restart file" << endl;
      else read(handle, m_sigma_xz, "data", grids);

#if CH_SPACEDIM==1
      group_id = handle.setGroup("sigmayx");
      if(group_id<0) cout << "Warning: sigmayx not found in restart file" << endl;
      else read(handle, m_sigma_yx, "data", grids);
      
      group_id = handle.setGroup("sigmayy");
      if(group_id<0) cout << "Warning: sigmayy not found in restart file" << endl;
      else read(handle, m_sigma_yy, "data", grids);
      
      group_id = handle.setGroup("sigmayz");
      if(group_id<0) cout << "Warning: sigmayz not found in restart file" << endl;
      else read(handle, m_sigma_yz, "data", grids);
#endif

#if CH_SPACEDIM<3
      group_id = handle.setGroup("sigmazx");
      if(group_id<0) cout << "Warning: sigmazx not found in restart file" << endl;
      else read(handle, m_sigma_zx, "data", grids);
      
      group_id = handle.setGroup("sigmazy");
      if(group_id<0) cout << "Warning: sigmazy not found in restart file" << endl;
      else read(handle, m_sigma_zy, "data", grids);
      
      group_id = handle.setGroup("sigmazz");
      if(group_id<0) cout << "Warning: sigmazz not found in restart file" << endl;
      else read(handle, m_sigma_zz, "data", grids);
#endif
 
      handle.close();

   }
   
}

void
PicSpeciesInterface::speciesPairingInit()
{
   ParmParse pp("pic_species.pairing_init");

   bool more_pairings(true);
   int pairing_number = -1;
   while(more_pairings) { // look for pic species...
 
      pairing_number += 1;
      stringstream s;
      s << "pic_species.pairing_init." << pairing_number; 
      ParmParse pp_spi( s.str().c_str() );
     
      string name;
      if(pp_spi.contains("name")) pp_spi.get("name",name);
      else more_pairings = false;

      if(more_pairings) {

         if(!procID()) cout << "found pairing_init with name = " << name << endl;

         int speciesA, speciesB;
         pp_spi.get("speciesA",speciesA);
         pp_spi.get("speciesB",speciesB);
         CH_assert(speciesA<m_num_species);
         CH_assert(speciesB<m_num_species);
         const PicSpeciesPtr spA(m_pic_species_ptr_vect[speciesA]);
         const PicSpeciesPtr spB(m_pic_species_ptr_vect[speciesB]);
         if(!procID()) cout << "species A = " << spA->name() << endl;
         int numP_A = spA->numParticles();
         int numP_B = spB->numParticles();
         if(!procID()) cout << "num A parts = " << numP_A << endl;
         if(!procID()) cout << "species B = " << spB->name() << endl;
         if(!procID()) cout << "num B parts = " << numP_B << endl;
         CH_assert(numP_A==numP_B);
        
         // get the mode amplitude and number 
         Real amplitude;
         pp_spi.get("amplitude",amplitude);
         std::vector<int> temp(SpaceDim,0);
         pp_spi.getarr( "mode", temp, 0, SpaceDim );
         IntVect mode = IntVect( temp );
         if(!procID()) cout << "amplitude = " << amplitude << endl;
         if(!procID()) cout << "mode = " << mode << endl;
   
         const RealVect Xmin = m_mesh.getXmin();
         const RealVect Xmax = m_mesh.getXmax();
         const RealVect L = Xmax-Xmin;
         const Real twoPi = Constants::TWOPI;
         
         // bin the particles up by grid cells
         spA->binTheParticles();
         spB->binTheParticles();
   
         // predefine some pointers to be used below
         JustinsParticle* partA_ptr = NULL;  
         JustinsParticle* partB_ptr = NULL;  
 
         LevelData<BinFab<JustinsParticlePtr>>& dataA_ptr = spA->partData_binfab();
         LevelData<BinFab<JustinsParticlePtr>>& dataB_ptr = spB->partData_binfab();
   
         const DisjointBoxLayout& grids = dataA_ptr.disjointBoxLayout();
         DataIterator dit(grids);
         for (dit.begin(); dit.ok(); ++dit) {

            BinFab<JustinsParticlePtr>& binFabA_ptr = dataA_ptr[dit];
            BinFab<JustinsParticlePtr>& binFabB_ptr = dataB_ptr[dit];
       
            std::vector<JustinsParticlePtr> vector_partA_ptrs;
            std::vector<JustinsParticlePtr> vector_partB_ptrs;
      
            const Box gridBox = grids.get(dit);
            BoxIterator bit(gridBox);
            for (bit.begin(); bit.ok(); ++bit) { // loop over cells
         
               const IntVect ig = bit(); // grid index

               List<JustinsParticlePtr>& cell_pListA = binFabA_ptr(ig,0);
               List<JustinsParticlePtr>& cell_pListB = binFabB_ptr(ig,0);
               int cell_numA = cell_pListA.length();
               int cell_numB = cell_pListB.length();
               CH_assert(cell_numA==cell_numB);
         
               // copy the iterators to a vector in order to loop
               vector_partA_ptrs.clear();
               vector_partA_ptrs.reserve(cell_numA);
               ListIterator<JustinsParticlePtr> litA(cell_pListA);
               for (litA.begin(); litA.ok(); ++litA) vector_partA_ptrs.push_back(litA());
 
               vector_partB_ptrs.clear();
               vector_partB_ptrs.reserve(cell_numB);
               ListIterator<JustinsParticlePtr> litB(cell_pListB);
               for (litB.begin(); litB.ok(); ++litB) vector_partB_ptrs.push_back(litB());
         
               for (int p=0; p<cell_numA; p++) { // loop over particle

                  // get particle data for particle A
                  JustinsParticlePtr& partA = vector_partA_ptrs[p];
                  partA_ptr = partA.getPointer();
                  RealVect& xpA = partA_ptr->position();
                  RealVect& xpoldA = partA_ptr->position_old();
                  
                  // get particle data for particle A
                  JustinsParticlePtr& partB = vector_partB_ptrs[p];
                  partB_ptr = partB.getPointer();
                  RealVect& xpB = partB_ptr->position();
                  RealVect& xpoldB = partB_ptr->position_old();
         
                  // reposition particle A, then put particle B at same place
                  RealVect xp0;
                  for (int dir=0; dir<SpaceDim; dir++) {
                     xp0[dir] = Xmin[dir] + MathUtils::rand()*L[dir];
                     Real arg = twoPi*mode[dir]*xp0[dir]/L[dir];
                     xpA[dir] = xp0[dir] + amplitude*cos(arg);
                     xpoldA[dir] = xpA[dir];
                     xpB[dir] = xpA[dir];
                     xpoldB[dir] = xpB[dir];
                  }
            
               }

            }

         }
         spA->applyBCs(false);
         spB->applyBCs(false);

         // don't forget to set pointers back to NULL and delete
         partA_ptr = NULL;
         partB_ptr = NULL;
         delete partA_ptr;
         delete partB_ptr;

         if(!procID()) cout << endl;
      }

   }

}

void PicSpeciesInterface::computeJfromMassMatrices( const ElectroMagneticFields&  a_emfields )
{
   CH_TIME("PicSpeciesInterface::computeJfromMassMatrices()");
   
   if (!a_emfields.advanceE()) return;
   
   const LevelData<EdgeDataBox>& electricField = a_emfields.getElectricField();
   const LevelData<NodeFArrayBox>& electricField_virtual = a_emfields.getVirtualElectricField();

#ifdef MASS_MATRIX_TEST
   LevelData<EdgeDataBox>& currentDensity = m_currentDensity_TEST;
   LevelData<NodeFArrayBox>& currentDensity_virtual = m_currentDensity_virtual_TEST;
#else
   LevelData<EdgeDataBox>& currentDensity = m_currentDensity;
   LevelData<NodeFArrayBox>& currentDensity_virtual = m_currentDensity_virtual;
#endif      

   // this calculation requires 2 ghost cells for cic
   const int nG = m_mesh.ghosts();
   CH_assert(nG>=2);

   const DisjointBoxLayout& grids(m_mesh.getDBL());

   for(DataIterator dit(grids); dit.ok(); ++dit) {
      
      Box node_box = surroundingNodes(grids[dit]);
      
      FArrayBox& Jx = currentDensity[dit][0];
      const FArrayBox& Ex = electricField[dit][0];
      const FArrayBox& Jx0 = m_J0[dit][0];
      const FArrayBox& Ex0 = m_E0[dit][0];
#if CH_SPACEDIM<3
      FArrayBox& Jv = currentDensity_virtual[dit].getFab();      
      const FArrayBox& Ev = electricField_virtual[dit].getFab();
      const FArrayBox& Jv0 = m_J0_virtual[dit].getFab();
      const FArrayBox& Ev0 = m_E0_virtual[dit].getFab();
#endif
#if CH_SPACEDIM>=2
      FArrayBox& Jy = currentDensity[dit][1];
      const FArrayBox& Ey = electricField[dit][1];
      const FArrayBox& Jy0 = m_J0[dit][1];
      const FArrayBox& Ey0 = m_E0[dit][1];
#endif
#if CH_SPACEDIM==3
      FArrayBox& Jz = currentDensity[dit][2];
      const FArrayBox& Ez = electricField[dit][2];
      const FArrayBox& Jz0 = m_J0[dit][2];
      const FArrayBox& Ez0 = m_E0[dit][2];
#endif

      //
      // compute Jx = Jx0 + sigma*E
      //

      const FArrayBox& sigxx = m_sigma_xx[dit][0];      
      const FArrayBox& sigxy = m_sigma_xy[dit][0];      
      const FArrayBox& sigxz = m_sigma_xz[dit][0];

      Box edge_box0 = enclosedCells(node_box,0);
      FORT_COMPUTE_JX_FROM_MASS_MATRIX( CHF_BOX(edge_box0),
                                        CHF_CONST_INTVECT(m_num_xx_comps),
                                        CHF_CONST_INTVECT(m_num_xy_comps),
                                        CHF_CONST_INTVECT(m_num_xz_comps),
                                        CHF_CONST_FRA(sigxx),
                                        CHF_CONST_FRA(sigxy),
                                        CHF_CONST_FRA(sigxz),
                                        CHF_CONST_FRA1(Ex0,0),
                                        CHF_CONST_FRA1(Ex,0),
#if CH_SPACEDIM==1
                                        CHF_CONST_FRA1(Ev0,0),
                                        CHF_CONST_FRA1(Ev,0),
                                        CHF_CONST_FRA1(Ev0,1),
                                        CHF_CONST_FRA1(Ev,1),
#elif CH_SPACEDIM==2
                                        CHF_CONST_FRA1(Ey0,0),
                                        CHF_CONST_FRA1(Ey,0),
                                        CHF_CONST_FRA1(Ev0,0),
                                        CHF_CONST_FRA1(Ev,0),
#elif CH_SPACEDIM==3
                                        CHF_CONST_FRA1(Ey0,0),
                                        CHF_CONST_FRA1(Ey,0),
                                        CHF_CONST_FRA1(Ez0,0),
                                        CHF_CONST_FRA1(Ez,0),
#endif
                                        CHF_CONST_FRA1(Jx0,0),
                                        CHF_FRA1(Jx,0) );

      //
      // compute Jy = Jy0 + sigma*E
      //

#if CH_SPACEDIM==1
      const FArrayBox& sigyx = m_sigma_yx[dit].getFab();
      const FArrayBox& sigyy = m_sigma_yy[dit].getFab();
      const FArrayBox& sigyz = m_sigma_yz[dit].getFab();
#else
      const FArrayBox& sigyx = m_sigma_xy[dit][1];      
      const FArrayBox& sigyy = m_sigma_xx[dit][1];      
      const FArrayBox& sigyz = m_sigma_xz[dit][1];
#endif

#if CH_SPACEDIM>1      
      Box edge_box1 = enclosedCells(node_box,1);
#else 
      Box edge_box1 = node_box;
#endif
      FORT_COMPUTE_JY_FROM_MASS_MATRIX( CHF_BOX(edge_box1),
                                        CHF_CONST_INTVECT(m_num_yx_comps),
                                        CHF_CONST_INTVECT(m_num_yy_comps),
                                        CHF_CONST_INTVECT(m_num_yz_comps),
                                        CHF_CONST_FRA(sigyx),
                                        CHF_CONST_FRA(sigyy),
                                        CHF_CONST_FRA(sigyz),
                                        CHF_CONST_FRA1(Ex0,0),
                                        CHF_CONST_FRA1(Ex,0),
#if CH_SPACEDIM==1
                                        CHF_CONST_FRA1(Ev0,0),
                                        CHF_CONST_FRA1(Ev,0),
                                        CHF_CONST_FRA1(Ev0,1),
                                        CHF_CONST_FRA1(Ev,1),
                                        CHF_CONST_FRA1(Jv0,0),
                                        CHF_FRA1(Jv,0) );
#elif CH_SPACEDIM==2
                                        CHF_CONST_FRA1(Ey0,0),
                                        CHF_CONST_FRA1(Ey,0),
                                        CHF_CONST_FRA1(Ev0,0),
                                        CHF_CONST_FRA1(Ev,0),
                                        CHF_CONST_FRA1(Jy0,0),
                                        CHF_FRA1(Jy,0) );
#elif CH_SPACEDIM==3
                                        CHF_CONST_FRA1(Ey0,0),
                                        CHF_CONST_FRA1(Ey,0),
                                        CHF_CONST_FRA1(Ez0,0),
                                        CHF_CONST_FRA1(Ez,0),
                                        CHF_CONST_FRA1(Jy0,0),
                                        CHF_FRA1(Jy,0) );
#endif

      //
      // compute Jz = Jz0 + sigma*E
      //

      const FArrayBox& sigzx = m_sigma_zx[dit].getFab();
      const FArrayBox& sigzy = m_sigma_zy[dit].getFab();
      const FArrayBox& sigzz = m_sigma_zz[dit].getFab();
      
#if CH_SPACEDIM==3     
      Box edge_box2 = enclosedCells(node_box,2);
#else 
      Box edge_box2 = node_box;
#endif
      FORT_COMPUTE_JZ_FROM_MASS_MATRIX( CHF_BOX(edge_box2),
                                        CHF_CONST_INTVECT(m_num_zx_comps),
                                        CHF_CONST_INTVECT(m_num_zy_comps),
                                        CHF_CONST_INTVECT(m_num_zz_comps),
                                        CHF_CONST_FRA(sigzx),
                                        CHF_CONST_FRA(sigzy),
                                        CHF_CONST_FRA(sigzz),
                                        CHF_CONST_FRA1(Ex0,0),
                                        CHF_CONST_FRA1(Ex,0),
#if CH_SPACEDIM==1
                                        CHF_CONST_FRA1(Ev0,0),
                                        CHF_CONST_FRA1(Ev,0),
                                        CHF_CONST_FRA1(Ev0,1),
                                        CHF_CONST_FRA1(Ev,1),
                                        CHF_CONST_FRA1(Jv0,1),
                                        CHF_FRA1(Jv,1) );
#elif CH_SPACEDIM==2
                                        CHF_CONST_FRA1(Ey0,0),
                                        CHF_CONST_FRA1(Ey,0),
                                        CHF_CONST_FRA1(Ev0,0),
                                        CHF_CONST_FRA1(Ev,0),
                                        CHF_CONST_FRA1(Jv0,0),
                                        CHF_FRA1(Jv,0) );
#elif CH_SPACEDIM==3
                                        CHF_CONST_FRA1(Ey0,0),
                                        CHF_CONST_FRA1(Ey,0),
                                        CHF_CONST_FRA1(Ez0,0),
                                        CHF_CONST_FRA1(Ez,0),
                                        CHF_CONST_FRA1(Jz0,0),
                                        CHF_FRA1(Jz,0) );
#endif

   }
   
#if CH_SPACEDIM==1
   // apply BC to J (needed for symmetry/axis BCs)
   // will eventually remove the need to do this through
   // the field_bc. This is legacy.
   const FieldBC& field_bc = a_emfields.getFieldBC();
   field_bc.prepJforBC( currentDensity, currentDensity_virtual,
                        m_J0, m_J0_virtual,
                        m_E0, m_E0_virtual,
                        electricField, electricField_virtual,
                        m_sigma_xx, m_sigma_xy, m_sigma_xz,
                        m_sigma_yx, m_sigma_yy, m_sigma_yz,
                        m_sigma_zx, m_sigma_zy, m_sigma_zz );
   field_bc.applyToJ( currentDensity, currentDensity_virtual );
#endif
   
   // divide by Jacobian after doing exchange (Corrected Jacobian does not have ghosts) 
   const LevelData<EdgeDataBox>& Jec = m_mesh.getCorrectedJec();  
   const LevelData<NodeFArrayBox>& Jnc = m_mesh.getCorrectedJnc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         currentDensity[dit][dir].divide(Jec[dit][dir],0,0,1);
      }
#if CH_SPACEDIM<3
      for (int comp=0; comp<currentDensity_virtual.nComp(); ++comp) {
         currentDensity_virtual[dit].getFab().divide(Jnc[dit].getFab(),0,comp,1);
      }
#endif
   }

   // numerical energy test will heat after 600 wpedt when using PC and
   // not calling exchange here (for virtual J only actually)... 
   // Not sure why exchange is needed here for PC to work correctly...
   SpaceUtils::exchangeEdgeDataBox(currentDensity);
#if CH_SPACEDIM<3
   SpaceUtils::exchangeNodeFArrayBox(currentDensity_virtual);
#endif
   
}


void PicSpeciesInterface::setChargeDensityOnNodes()
{
   CH_TIME("PicSpeciesInterface::setChargeDensityOnNodes()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   // set the charge density to zero
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& this_rho = m_chargeDensity_nodes[dit].getFab();
      this_rho.setVal(0.0);
   }

   // loop over all pic species and add the charge density 
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); ++sp) {

      const PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if(species->charge() == 0) continue;

      species->setChargeDensityOnNodes();
      const LevelData<NodeFArrayBox>& species_rho = species->getChargeDensityOnNodes();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& this_species_rho = species_rho[dit].getFab();
         FArrayBox& this_rho = m_chargeDensity_nodes[dit].getFab();
         this_rho.plus(this_species_rho,0,0,m_chargeDensity_nodes.nComp());
      }

   }

}

void PicSpeciesInterface::setCurrentDensity( const bool  a_from_explicit_solver )
{
   CH_TIME("PicSpeciesInterface::setCurrentDensity()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
 
   // set the current density member to zero
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         m_currentDensity[dit][dir].setVal(0.0);
      }
#if CH_SPACEDIM<3
      FArrayBox& this_Jv_nodes = m_currentDensity_virtual[dit].getFab();
      this_Jv_nodes.setVal(0.0);
#endif
   }

   // loop over all pic species and add the current density 
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); ++sp) {
      const PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if(species->charge()==0) continue;

      species->setCurrentDensity( a_from_explicit_solver );
      const LevelData<EdgeDataBox>& species_J = species->getCurrentDensity();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
            const FArrayBox& this_species_Jdir = species_J[dit][dir];
            m_currentDensity[dit][dir].plus(this_species_Jdir,0,0,1);
         }
      }
#if CH_SPACEDIM<3
      const LevelData<NodeFArrayBox>& species_Jv = species->getCurrentDensity_virtual();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& this_species_Jv = species_Jv[dit].getFab();
         FArrayBox& this_Jv = m_currentDensity_virtual[dit].getFab();
         this_Jv.plus(this_species_Jv,0,0,m_currentDensity_virtual.nComp());
         //if(!procID()) SpaceUtils::inspectFArrayBox( this_Jv, this_Jv.box(), 1);
      }
#endif
   }

   // call exchange 
   // JRA 1-13-22
   // When I call exchange for virtual J here, the numerical energy regression test 
   // that uses jfnk fails. However, when using jfnk it is possible that a machine-level
   // difference can grow into a tolerance-level difference
   //SpaceUtils::exchangeEdgeDataBox(m_currentDensity);
   //if(SpaceDim<3) SpaceUtils::exchangeNodeFArrayBox(m_currentDensity_virtual);

}

void PicSpeciesInterface::preRHSOp( const bool                    a_from_emjacobian, 
                                    const ElectroMagneticFields&  a_emfields,
                                    const Real                    a_dt,
                                    const int                     a_nonlinear_iter )
{
   CH_TIME("PicSpeciesInterface::preRHSOp()");
  
   if(a_from_emjacobian && m_freeze_particles_jacobian) return;
  
   // half dt advance of particle positions and velocities
   // and then compute current density at half time step
   //
   // xbar = xn + dt/2*vbar
   // upbar = upn + dt/2*q/m*(E(xbar) + upbar/gammapbar x B(xbar))
   // Jbar = Sum_sSum_p(qp*S(xbar-xg)*upbar/gammap)/dV

   if(a_from_emjacobian) { // called from linear stage of jfnk

      if (m_use_mass_matrices) {

         computeJfromMassMatrices( a_emfields );

      }
      else {

         for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
            auto species(m_pic_species_ptr_vect[sp]);
            if (m_quasi_freeze_particles_jacobian) {
               species->interpolateFieldsToParticles( a_emfields );
               species->addExternalFieldsToParticles( a_emfields ); 
               species->advanceVelocities( a_dt, true );
            }
            else {
               species->advanceParticlesIteratively( a_emfields, a_dt );
            }
         }
         setCurrentDensity();
#ifdef MASS_MATRIX_TEST
         setMassMatrices( a_emfields, a_dt );
         computeJfromMassMatrices( a_emfields );
#endif

      }

   }
   else { // called from nonlinear stage of jfnk or from picard solver
  
      for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
         auto species(m_pic_species_ptr_vect[sp]);
         if(a_nonlinear_iter==0 && m_mod_init_advance) {
            species->advancePositionsImplicit( a_dt, true );
            species->interpolateFieldsToParticles( a_emfields );
            species->addExternalFieldsToParticles( a_emfields ); 
            species->advanceVelocities( a_dt, true );
         }
         else species->advanceParticlesIteratively( a_emfields, a_dt );
      }
     
      if(m_use_mass_matrices) {
#ifdef MASS_MATRIX_COST_TEST
         for (int iter=0; iter<2; iter++) {
            for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
               auto species(m_pic_species_ptr_vect[sp]);
               species->interpolateFieldsToParticles( a_emfields );
               species->advanceVelocities( a_dt, true );
            }
         }
         setCurrentDensity();
#endif
         setMassMatrices( a_emfields, a_dt );
         computeJfromMassMatrices( a_emfields );
      }
      else {
         setCurrentDensity();
#ifdef MASS_MATRIX_TEST
         setMassMatrices( a_emfields, a_dt );
         computeJfromMassMatrices( a_emfields );
#endif
      }

   }
   
   LevelData<EdgeDataBox>& J = m_currentDensity;
   LevelData<NodeFArrayBox>& Jv = m_currentDensity_virtual;
   addInflowJ( J, Jv, a_emfields, a_dt );   

}
 
void PicSpeciesInterface::postNewtonUpdate( const ElectroMagneticFields&  a_emfields,
                                            const Real                    a_time,
                                            const Real                    a_dt )
{
   CH_TIME("PicSpeciesInterface::postNewtonUpdate()");
  
   // this is called on last Newton update for energy-conserving
   // semi-implicit option from theta implicit integrator. To get
   // exact energy conservation, need to update particle 
   // Epbar only and then re-compute particle vpbar immediately 
   // following the newton correction update to the fields. 
   // This is done by terminating the newton solve with iter_max 
   // and calling this function rather than preRHSOp() upon the 
   // final iteration in the newton solve(). 
        
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
      auto species(m_pic_species_ptr_vect[sp]);
      species->interpolateEfieldToParticles( a_emfields );
      species->advanceVelocities( a_dt, true );
      species->advancePositionsImplicit( a_dt, true );
   }
 
}

void PicSpeciesInterface::setMassMatrices( const ElectroMagneticFields&  a_emfields,
                                           const Real  a_dt )
{
   CH_TIME("PicSpeciesInterface::setMassMatrices()");
   
   // set all to zero
   SpaceUtils::zero( m_J0 );
   SpaceUtils::zero( m_J0_virtual );
   
   SpaceUtils::zero( m_sigma_xx );
   SpaceUtils::zero( m_sigma_xy );
   SpaceUtils::zero( m_sigma_xz );
#if CH_SPACEDIM==1
   SpaceUtils::zero( m_sigma_yx );
   SpaceUtils::zero( m_sigma_yy );
   SpaceUtils::zero( m_sigma_yz );
#endif
   SpaceUtils::zero( m_sigma_zx );
   SpaceUtils::zero( m_sigma_zy );
   SpaceUtils::zero( m_sigma_zz );

   // loop over all pic species and add their contribution 
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); ++sp) {
      const PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      species->accumulateMassMatrices( m_sigma_xx, m_sigma_xy, m_sigma_xz,
#if CH_SPACEDIM==1
                                       m_sigma_yx, m_sigma_yy, m_sigma_yz,
#endif
                                       m_sigma_zx, m_sigma_zy, m_sigma_zz, 
                                       m_J0, m_J0_virtual, a_emfields, a_dt );
   }
   
   // add ghost cells to valid cells
   LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   m_J0.exchange( m_J0.interval(), m_mesh.reverseCopier(), addEdgeOp );
   m_sigma_xx.exchange( m_sigma_xx.interval(), m_mesh.reverseCopier(), addEdgeOp );
   m_sigma_xy.exchange( m_sigma_xy.interval(), m_mesh.reverseCopier(), addEdgeOp );
   m_sigma_xz.exchange( m_sigma_xz.interval(), m_mesh.reverseCopier(), addEdgeOp );

   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   m_J0_virtual.exchange( m_J0_virtual.interval(), m_mesh.reverseCopier(), addNodeOp );
#if CH_SPACEDIM==1
   m_sigma_yx.exchange( m_sigma_yx.interval(), m_mesh.reverseCopier(), addNodeOp );
   m_sigma_yy.exchange( m_sigma_yy.interval(), m_mesh.reverseCopier(), addNodeOp );
   m_sigma_yz.exchange( m_sigma_yz.interval(), m_mesh.reverseCopier(), addNodeOp );
#endif
   m_sigma_zx.exchange( m_sigma_zx.interval(), m_mesh.reverseCopier(), addNodeOp );
   m_sigma_zy.exchange( m_sigma_zy.interval(), m_mesh.reverseCopier(), addNodeOp );
   m_sigma_zz.exchange( m_sigma_zz.interval(), m_mesh.reverseCopier(), addNodeOp );
   
   // save electric field to compute perturbed E during GMRES
   const LevelData<EdgeDataBox>& electricField = a_emfields.getElectricField();
   const LevelData<NodeFArrayBox>& electricField_virtual = a_emfields.getVirtualElectricField();
   
#ifdef MASS_MATRIX_TEST
   SpaceUtils::zero( m_E0 );
   SpaceUtils::zero( m_E0_virtual );
#else
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         m_E0[dit][dir].copy(electricField[dit][dir]);
      }
      m_E0_virtual[dit].getFab().copy(electricField_virtual[dit].getFab());
   } 
#endif
   
}

void PicSpeciesInterface::addInflowJ( LevelData<EdgeDataBox>&    a_J,
                                      LevelData<NodeFArrayBox>&  a_Jv,
                                const ElectroMagneticFields&     a_emfields,
                                const Real                       a_dt )
{
   CH_TIME("PicSpeciesInterface::addInflowJ()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   // if using the sub-orbit model for inflow particles, update the inflow particle
   // quantities, compute their current, and add it to the total J

   // loop over all pic species and add the current density from inflow particles 
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); ++sp) {
      const PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if(species->charge()==0) continue;
      if(!species->suborbit_inflowJ()) continue;

      // update the inflow particle quantities     
      species->advanceInflowParticlesIteratively( a_emfields, a_dt );
      
      // compute inflow J and add it to the total J     
      species->setInflowJ( a_dt );
      const LevelData<EdgeDataBox>& this_inflowJ = species->getInflowJ();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
            const FArrayBox& inflowJ_dir = this_inflowJ[dit][dir];
            a_J[dit][dir].plus(inflowJ_dir,0,0,1);
         }
      }
#if CH_SPACEDIM<3
      const LevelData<NodeFArrayBox>& this_inflowJv = species->getInflowJ_virtual();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& inflowJv = this_inflowJv[dit].getFab();
         a_Jv[dit].getFab().plus(inflowJv,0,0,a_Jv.nComp());
      }
#endif
   }

}

void PicSpeciesInterface::prepForScatter( const int   a_num_coulomb,
	                                  const bool  a_from_scatterDt )
{
   CH_TIME("PicSpeciesInterface::prepForScatter()");
   
   // prepare all species to be scattered for scattering
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
      PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      species->binTheParticles();
      species->setNumberDensityFromBinFab();
      if(a_from_scatterDt) {
         species->setMomentumDensityFromBinFab();
         species->setEnergyDensityFromBinFab();
      }
      else if(a_num_coulomb>0 && species->charge()!=0) {
         species->setMomentumDensityFromBinFab();
         species->setEnergyDensityFromBinFab();
      }
   }
   
   if(a_num_coulomb>0) setDebyeLength();

}

void PicSpeciesInterface::setDebyeLength() 
{
   CH_TIME("PicSpeciesInterface::setDebyeLength()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());

   Real mcSq_eV = Constants::ME*Constants::CVAC*Constants::CVAC*Constants::EV_PER_JOULE;
   SpaceUtils::zero( m_DebyeLength );
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {

      PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if(species->charge()==0) continue;

      // job of caller to make sure the momentus are precomputed
      const LevelData<FArrayBox>& numDen = species->getNumberDensityFromBinFab();
      const LevelData<FArrayBox>& momDen = species->getMomentumDensityFromBinFab();
      const LevelData<FArrayBox>& eneDen = species->getEnergyDensityFromBinFab();
      const Real Aconst = Constants::EP0/Constants::QE/(species->charge()*species->charge());
   
      for(DataIterator dit(grids); dit.ok(); ++dit) {
   
               FArrayBox& this_LDe = m_DebyeLength[dit];
         const FArrayBox& this_numDen = numDen[dit];
         const FArrayBox& this_momDen = momDen[dit];
         const FArrayBox& this_eneDen = eneDen[dit];
     
         const Box gridBox = grids.get(dit);
         BoxIterator gbit(gridBox);
         for (gbit.begin(); gbit.ok(); ++gbit) {

            const IntVect ig = gbit();
       
            Real N = this_numDen.get(ig,0); // [1/m^3]
            if(N == 0.0) continue;
            Real rho = N*species->mass();

	    Real rhoUx = this_momDen.get(ig,0); // [m/s*1/m^3]
	    Real rhoUy = this_momDen.get(ig,1); // [m/s*1/m^3]
	    Real rhoUz = this_momDen.get(ig,2); // [m/s*1/m^3]
	    const Real meanE = (rhoUx*rhoUx + rhoUy*rhoUy + rhoUz*rhoUz)/rho/2.0;
         
            Real rhoE = 0.0;
            for( int dir=0; dir<3; dir++) rhoE += this_eneDen.get(ig,dir); 
            
            Real T_eV = 2.0/3.0*(rhoE - meanE)/N*mcSq_eV;
	    T_eV = std::max(T_eV, 0.01);
            //if(!procID()) cout << "JRA: species = " << species->name() << endl;
            //if(!procID()) cout << "JRA: numDen = " << N << endl;
            //if(!procID()) cout << "JRA: meanE[ig="<<ig<<"] = " << meanE << endl;
            //if(!procID()) cout << "JRA: T_eV[ig="<<ig<<"]  = " << T_eV << endl;
	    const Real LDe_sq = Aconst*T_eV/N; // [m^2] 

	    // store sum of 1/LDe^2 in LDe container
	    Real sumLDe_sq_inv = this_LDe.get(ig,0);
	    sumLDe_sq_inv += 1.0/LDe_sq;
	    this_LDe.set(ig,0,sumLDe_sq_inv);

	 }

      }

   }

   // invert sum of species inverse Debye lengths squared and take square root
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& this_LDe = m_DebyeLength[dit];
      const Box gridBox = grids.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) {
         const IntVect ig = gbit();
         Real LDe_sq_inv = this_LDe.get(ig,0);
	 Real LDe = 1.0/std::sqrt(LDe_sq_inv);
	 this_LDe.set(ig,0,LDe);
         //if(!procID()) cout << "JRA: LDe[ig="<<ig<<"] = " << LDe << " m" << endl;
      }
   }

} 

Real 
PicSpeciesInterface::courantDt()
{

   m_courant_dt = DBL_MAX;
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
      PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      species->setStableDt();
      m_courant_dt = std::min(m_courant_dt,species->stableDt());
   }
   
//   Real m_courant_dt = courantDt_local;
//#ifdef CH_MPI
//   MPI_Allreduce( &courantDt_local, &m_courant_dt, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD ); 
//#endif
  
   return m_courant_dt;

}

#include "NamespaceFooter.H"

