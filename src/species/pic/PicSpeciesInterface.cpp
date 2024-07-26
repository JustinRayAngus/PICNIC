
#include "PicSpeciesInterface.H"
#include "ParmParse.H"
#include "CH_HDF5.H"

#include "FieldsF_F.H"

#include "MathUtils.H"
#include "PicnicConstants.H"
#include "SpaceUtils.H"

#include "NamespaceHeader.H"

PicSpeciesInterface::PicSpeciesInterface( const CodeUnits&   a_units,
                                          const DomainGrid&  a_mesh )
   : m_verbosity(true),
     m_writeSpeciesChargeDensity(false),
     m_writeSpeciesSurfaceCharge(false),
     m_writeSpeciesCurrentDensity(false),
     m_writeSpeciesNppc(false),
     m_writeSpeciesEnergyOffDiag(false),
     m_writeSpeciesEnergyFlux(false),
     m_verbose_particles(false),
     m_part_order_swap(false),
     m_iter_max_particles(0),
     m_rtol_particles(1.0e-12),
     m_newton_maxits(20),
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
   pp.query("write_species_surface_charge", m_writeSpeciesSurfaceCharge);
   pp.query("write_species_current_density", m_writeSpeciesCurrentDensity);
   pp.query("write_species_nppc", m_writeSpeciesNppc);
   pp.query("write_species_energy_off_diagonal", m_writeSpeciesEnergyOffDiag);
   pp.query("write_species_energy_flux", m_writeSpeciesEnergyFlux);
   if(m_writeSpeciesEnergyFlux) m_writeSpeciesEnergyOffDiag = true;
   //
   pp.query("part_order_swap",m_part_order_swap);
   pp.query("iter_max_particles",m_iter_max_particles);
   pp.query("verbose_particles",m_verbose_particles);
   pp.query("rtol_particles",m_rtol_particles);
   pp.query("newton_maxits",m_newton_maxits);
   pp.query("newton_num_guess",m_newton_num_guess);
   //
   pp.query("freeze_particles_jacobian",m_freeze_particles_jacobian);
   pp.query("quasi_freeze_particles_jacobian",m_quasi_freeze_particles_jacobian);
   pp.query("use_mass_matrices",m_use_mass_matrices);
   pp.query("mod_init_advance",m_mod_init_advance);
      
   if(m_use_mass_matrices) m_quasi_freeze_particles_jacobian = true;
   if(m_quasi_freeze_particles_jacobian) m_freeze_particles_jacobian=false;
   
   createAllPicSpecies( a_mesh );
  
#ifdef PRINT_COMPS
   m_print_comps = false;
   if(!procID()) m_print_comps = true;
#endif
      
   m_num_species = m_pic_species_ptr_vect.size();
   
   // define the LevelData members for moments
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const int ghosts(m_mesh.ghosts());
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   m_chargeDensity.define(grids,1,ghostVect);
   m_chargeDensity_faces.define(grids,1,ghostVect);
   m_chargeDensity_nodes.define(grids,1,ghostVect);
   m_surfaceCharge_nodes.define(grids,1,ghostVect);
   
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
      
         species->setParticleSolverParams( m_verbose_particles,
			                   m_part_order_swap,
                                           m_iter_max_particles,
                                           m_rtol_particles,
                                           m_newton_maxits,
                                           m_newton_num_guess );

      }
  
   }
    
   // additional species-pair initialization
   speciesPairingInit();

   if(!a_implicit_advance) { m_use_mass_matrices = false; }

   if(!procID()) {
      cout << "Finished initializing " << m_pic_species_ptr_vect.size();
      cout << " pic species" << endl << endl;
      cout << " write species charge density  = " << (m_writeSpeciesChargeDensity?"true":"false") << endl;
      cout << " write species surface charge  = " << (m_writeSpeciesSurfaceCharge?"true":"false") << endl;
      cout << " write species current density = " << (m_writeSpeciesCurrentDensity?"true":"false") << endl;
      cout << " write species nppc = " << (m_writeSpeciesNppc?"true":"false") << endl;
      cout << " write species energy off diagonal = " << (m_writeSpeciesEnergyOffDiag?"true":"false") << endl;
      cout << " write species energy flux = " << (m_writeSpeciesEnergyFlux?"true":"false") << endl;
      if (a_implicit_advance) {
         cout << " implicit advance parameters:" << endl;
         cout << "  verbose_particles = " << m_verbose_particles << endl;
         cout << "  part_order_swap = " << m_part_order_swap << endl;
         cout << "  iter_max_particles = " << m_iter_max_particles << endl;
         cout << "  rtol_particles = " << m_rtol_particles << endl;
         cout << "  newton_maxits = " << m_newton_maxits << endl;
         cout << "  newton_num_guess = " << m_newton_num_guess << endl;
         cout << "  freeze_particles_jacobian = " << (m_freeze_particles_jacobian?"true":"false") << endl;
         cout << "  quasi_freeze_particles_jacobian = " << (m_quasi_freeze_particles_jacobian?"true":"false") << endl;
         cout << "  use_mass_matrices = " << (m_use_mass_matrices?"true":"false") << endl;
         cout << "  mod_init_advance = " << (m_mod_init_advance?"true":"false") << endl << endl;
      }
   }

}

void PicSpeciesInterface::initializeMassMatrices( EMFields&     a_emfields,
	                                    const Real&         a_dt,	
		                            const std::string&  a_restart_file_name )
{
   CH_TIME("PicSpeciesInterface::initializeMassMatrices()");

   bool init_mass_matrices = false;
   if(m_use_mass_matrices) init_mass_matrices = true;
#ifdef MASS_MATRIX_TEST
   CH_assert(!m_use_mass_matrices);
   init_mass_matrices = true;
#endif

   if(!init_mass_matrices) return;
   
   if(!procID()) cout << "Initializing mass matrices..." << endl;

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
         m_ncomp_xx[dir] = 3 + tsc;
         m_ncomp_xy[dir] = 4 + tsc;
         m_ncomp_xz[dir] = 4 + tsc;
         //
         m_ncomp_yx[dir] = 4 + tsc;
         m_ncomp_yy[dir] = 3 + tsc;
         m_ncomp_yz[dir] = 3 + tsc;
         //
         m_ncomp_zx[dir] = 4 + tsc;
         m_ncomp_zy[dir] = 3 + tsc;
         m_ncomp_zz[dir] = 3 + tsc;
      }
      if(dir==1) {
         m_ncomp_xx[dir] = 3 + tsc;
         m_ncomp_xy[dir] = 4 + tsc;
         m_ncomp_xz[dir] = 3 + tsc;
         //
         m_ncomp_yx[dir] = 4 + tsc;
         m_ncomp_yy[dir] = 3 + tsc;
         m_ncomp_yz[dir] = 4 + tsc;
         //
         m_ncomp_zx[dir] = 3 + tsc;
         m_ncomp_zy[dir] = 4 + tsc;
         m_ncomp_zz[dir] = 3 + tsc;
      }
      if(dir==2) {
         m_ncomp_xx[dir] = 3 + tsc;
         m_ncomp_xy[dir] = 3 + tsc;
         m_ncomp_xz[dir] = 4 + tsc;
         //
         m_ncomp_yx[dir] = 3 + tsc;
         m_ncomp_yy[dir] = 3 + tsc;
         m_ncomp_yz[dir] = 4 + tsc;
         //
         m_ncomp_zx[dir] = 4 + tsc;
         m_ncomp_zy[dir] = 4 + tsc;
         m_ncomp_zz[dir] = 3 + tsc;
      }
   }
   
   if(interp_type==CIC) {
      CH_assert(ghosts>=2);
   }
   if(interp_type==CC0) {
      m_ncomp_xx[0] = 5; // 2 cell crossings permitted
      CH_assert(ghosts>=2);
      CH_assert(SpaceDim==1);
   }
   if(interp_type==CC1 && SpaceDim==1) {
      CH_assert(ghosts>=2);
      int maxXings = ghosts - 1;
      m_ncomp_xx[0] = 3 + 2*maxXings;
      m_ncomp_xy[0] = 2 + 2*maxXings;
      m_ncomp_xz[0] = 2 + 2*maxXings;
      m_ncomp_yx[0] = 2 + 2*maxXings;
      m_ncomp_zx[0] = 2 + 2*maxXings;
   }
   if(interp_type==CC1 && SpaceDim==2) {
      CH_assert(ghosts>=3);
      int maxXings = ghosts - 2;
      for (int dir=0; dir<SpaceDim; ++dir) {
         if(dir==0) {
            m_ncomp_xx[dir] = 3 + 2*maxXings; 
            m_ncomp_xy[dir] = 4 + 2*maxXings; 
            m_ncomp_xz[dir] = 2 + 2*maxXings;
            //
            m_ncomp_yx[dir] = 4 + 2*maxXings; 
            m_ncomp_yy[dir] = 5 + 2*maxXings; 
            m_ncomp_yz[dir] = 3 + 2*maxXings;
            //
            m_ncomp_zx[dir] = 2 + 2*maxXings; 
            m_ncomp_zy[dir] = 3 + 2*maxXings; 
            m_ncomp_zz[dir] = 3;
         } 
         if(dir==1) {
            m_ncomp_xx[dir] = 5 + 2*maxXings; 
            m_ncomp_xy[dir] = 4 + 2*maxXings; 
            m_ncomp_xz[dir] = 3 + 2*maxXings;
            //
            m_ncomp_yx[dir] = 4 + 2*maxXings; 
            m_ncomp_yy[dir] = 3 + 2*maxXings; 
            m_ncomp_yz[dir] = 2 + 2*maxXings;
            //
            m_ncomp_zx[dir] = 3 + 2*maxXings; 
            m_ncomp_zy[dir] = 2 + 2*maxXings; 
            m_ncomp_zz[dir] = 3;
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
   m_sigma_xx.define(grids,m_ncomp_xx.product(),ghostVect);
   m_sigma_xy.define(grids,m_ncomp_xy.product(),ghostVect);
   m_sigma_xz.define(grids,m_ncomp_xz.product(),ghostVect);   
#if CH_SPACEDIM==1
   m_sigma_yx.define(grids,m_ncomp_yx.product(),ghostVect);
   m_sigma_yy.define(grids,m_ncomp_yy.product(),ghostVect);
   m_sigma_yz.define(grids,m_ncomp_yz.product(),ghostVect);
#endif
   m_sigma_zx.define(grids,m_ncomp_zx.product(),ghostVect);
   m_sigma_zy.define(grids,m_ncomp_zy.product(),ghostVect);
   m_sigma_zz.define(grids,m_ncomp_zz.product(),ghostVect);
   
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

   if(!procID()) {
      cout << " ncomp_xx = " << m_ncomp_xx << endl;
      cout << " ncomp_xy = " << m_ncomp_xy << endl;
      cout << " ncomp_xz = " << m_ncomp_xz << endl;
      cout << " ncomp_yx = " << m_ncomp_yx << endl;
      cout << " ncomp_yy = " << m_ncomp_yy << endl;
      cout << " ncomp_yz = " << m_ncomp_yz << endl;
      cout << " ncomp_zx = " << m_ncomp_zx << endl;
      cout << " ncomp_zy = " << m_ncomp_zy << endl;
      cout << " ncomp_zz = " << m_ncomp_zz << endl;
      cout << "Finished initializing mass matrices" << endl << endl;
   }
	 
   a_emfields.initializeMassMatricesForPC( m_ncomp_xx, m_ncomp_xy, m_ncomp_xz,
	 		                   m_ncomp_yx, m_ncomp_yy, m_ncomp_yz,
         			           m_ncomp_zx, m_ncomp_zy, m_ncomp_zz );
   
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
	 const Real dummy_time = 0.0;
         spA->applyBCs(false,dummy_time);
         spB->applyBCs(false,dummy_time);

         // don't forget to set pointers back to NULL and delete
         partA_ptr = NULL;
         partB_ptr = NULL;
         delete partA_ptr;
         delete partB_ptr;

         if(!procID()) cout << endl;
      }

   }

}

void PicSpeciesInterface::computeJfromMassMatrices( const EMFields&  a_emfields )
{
   CH_TIME("PicSpeciesInterface::computeJfromMassMatrices()");
   
   if (!a_emfields.advanceE()) return;
   
   const LevelData<EdgeDataBox>& electricField = a_emfields.getFilteredElectricField();
   const LevelData<NodeFArrayBox>& electricField_virtual = a_emfields.getVirtualElectricField();

#ifdef MASS_MATRIX_TEST
   LevelData<EdgeDataBox>& currentDensity = m_currentDensity_TEST;
   LevelData<NodeFArrayBox>& currentDensity_virtual = m_currentDensity_virtual_TEST;
#else
   LevelData<EdgeDataBox>& currentDensity = m_currentDensity;
   LevelData<NodeFArrayBox>& currentDensity_virtual = m_currentDensity_virtual;
#endif      
   

   // zero out the containers for J
   SpaceUtils::zero( currentDensity );
   SpaceUtils::zero( currentDensity_virtual );


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
      // compute Jx = Jx0 + sigma*(E-E0)
      //

      const FArrayBox& sigxx = m_sigma_xx[dit][0];      
      const FArrayBox& sigxy = m_sigma_xy[dit][0];      
      const FArrayBox& sigxz = m_sigma_xz[dit][0];

      Box edge_box0 = Jx.box();
      //if(!procID()) cout << "Jx.box() = " << edge_box0 << endl;
      FORT_COMPUTE_JX_FROM_MASS_MATRIX( CHF_BOX(edge_box0),
                                        CHF_CONST_INTVECT(m_ncomp_xx),
                                        CHF_CONST_INTVECT(m_ncomp_xy),
                                        CHF_CONST_INTVECT(m_ncomp_xz),
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
      // compute Jy = Jy0 + sigma*(E-E0)
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

#if CH_SPACEDIM==1
      Box edge_box1 = Jv.box();
#else
      Box edge_box1 = Jy.box();
#endif
      //if(!procID()) cout << "Jy.box() = " << edge_box1 << endl;
      FORT_COMPUTE_JY_FROM_MASS_MATRIX( CHF_BOX(edge_box1),
                                        CHF_CONST_INTVECT(m_ncomp_yx),
                                        CHF_CONST_INTVECT(m_ncomp_yy),
                                        CHF_CONST_INTVECT(m_ncomp_yz),
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
      // compute Jz = Jz0 + sigma*(E-E0)
      //

      const FArrayBox& sigzx = m_sigma_zx[dit].getFab();
      const FArrayBox& sigzy = m_sigma_zy[dit].getFab();
      const FArrayBox& sigzz = m_sigma_zz[dit].getFab();
      
#if CH_SPACEDIM<3
      Box edge_box2 = Jv.box();
#else
      Box edge_box2 = Jz.box();
#endif
      //if(!procID()) cout << "Jz.box() = " << edge_box2 << endl;
      FORT_COMPUTE_JZ_FROM_MASS_MATRIX( CHF_BOX(edge_box2),
                                        CHF_CONST_INTVECT(m_ncomp_zx),
                                        CHF_CONST_INTVECT(m_ncomp_zy),
                                        CHF_CONST_INTVECT(m_ncomp_zz),
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

}

void PicSpeciesInterface::finalizeSettingJ( LevelData<EdgeDataBox>&    a_currentDensity,
		                            LevelData<NodeFArrayBox>&  a_currentDensity_virtual,
				      const EMFields&                  a_emfields )
{
   CH_TIME("PicSpeciesInterface::finalizeSettingJ()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   // apply BCs to J
   const FieldBC& field_bc = a_emfields.getFieldBC();
   field_bc.applyToJ( a_currentDensity, a_currentDensity_virtual );

   // perform addOp exchange on J
   LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   a_currentDensity.exchange( a_currentDensity.interval(), m_mesh.reverseCopier(), addEdgeOp );
#if CH_SPACEDIM<3
   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   a_currentDensity_virtual.exchange( a_currentDensity_virtual.interval(), m_mesh.reverseCopier(), addNodeOp );
#endif
  
   // Here is where we will apply filtering

   // divide by corrected Jacobian
   const LevelData<EdgeDataBox>& Jec = m_mesh.getCorrectedJec();  
   const LevelData<NodeFArrayBox>& Jnc = m_mesh.getCorrectedJnc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         a_currentDensity[dit][dir].divide(Jec[dit][dir],0,0,1);
      }
#if CH_SPACEDIM<3
      for (int comp=0; comp<a_currentDensity_virtual.nComp(); ++comp) {
        a_currentDensity_virtual[dit].getFab().divide(Jnc[dit].getFab(),0,comp,1);
      }
#endif
   }

}

void PicSpeciesInterface::setChargeDensityOnNodes( const bool  a_use_filtering )
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

      species->setChargeDensityOnNodes( a_use_filtering );
      const LevelData<NodeFArrayBox>& species_rho = species->getChargeDensityOnNodes();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& this_species_rho = species_rho[dit].getFab();
         FArrayBox& this_rho = m_chargeDensity_nodes[dit].getFab();
         this_rho.plus(this_species_rho,0,0,m_chargeDensity_nodes.nComp());
      }

   }

}
  
void PicSpeciesInterface::setCurrentDensity( const EMFields&  a_emfields,
                                             const Real       a_dt,
		                             const bool       a_finalizeJ,
		                             const bool       a_from_explicit_solver )
{
   CH_TIME("PicSpeciesInterface::setCurrentDensity()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
 
   // set the current density member to zero
   SpaceUtils::zero( m_currentDensity );
#if CH_SPACEDIM<3
   SpaceUtils::zero( m_currentDensity_virtual );
#endif

   // loop over all pic species and add the current density 
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); ++sp) {
      const PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if(species->charge()==0) continue;
      species->setCurrentDensity( a_dt, a_from_explicit_solver );
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
      }
#endif
   }
   
   if(a_finalizeJ) { // apply BCs, do addOp exchange, then divide by corrected Jacobian
      finalizeSettingJ( m_currentDensity, m_currentDensity_virtual, a_emfields );
   }

}

void PicSpeciesInterface::setSurfaceChargeOnNodes( const bool  a_use_filtering )
{
   CH_TIME("PicSpeciesInterface::setSurfaceChargeOnNodes()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   // set the surface charge to zero
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& this_sigma = m_surfaceCharge_nodes[dit].getFab();
      this_sigma.setVal(0.0);
   }

   // loop over all pic species and add the surface charge
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); ++sp) {

      const PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if(species->charge() == 0) continue;

      const LevelData<NodeFArrayBox>& species_sigma = species->getSurfaceChargeOnNodes();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& this_species_sigma = species_sigma[dit].getFab();
         FArrayBox& this_sigma = m_surfaceCharge_nodes[dit].getFab();
         this_sigma.plus(this_species_sigma,0,0,m_surfaceCharge_nodes.nComp());
      }

   }

#if CH_SPACEDIM>1
   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   m_surfaceCharge_nodes.exchange( m_surfaceCharge_nodes.interval(), 
                                   m_mesh.reverseCopier(), addNodeOp );
   SpaceUtils::exchangeNodeFArrayBox( m_surfaceCharge_nodes );
#endif

}

void PicSpeciesInterface::preRHSOp( const bool       a_from_emjacobian, 
                                    const EMFields&  a_emfields,
                                    const Real       a_dt,
                                    const int        a_nonlinear_iter )
{
   CH_TIME("PicSpeciesInterface::preRHSOp()");
  
   if(a_from_emjacobian && m_freeze_particles_jacobian) return;
  
   // half dt advance of particle positions and velocities
   // and then compute current density at half time step
   //
   // xbar = xn + dt/2*vbar
   // upbar = upn + dt/2*q/m*(E(xbar) + upbar/gammapbar x B(xbar))
   // Jbar = Sum_sSum_p(qp*S(xbar-xg)*upbar/gammap)/dV

   if(a_emfields.useFiltering()) a_emfields.setFilteredFields();

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
         setCurrentDensity( a_emfields, a_dt, false );
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
            species->advancePositionsImplicit( a_dt );
            species->interpolateFieldsToParticles( a_emfields );
            species->addExternalFieldsToParticles( a_emfields ); 
            species->advanceVelocities( a_dt, true );
#ifdef NEW_EXACT_CHARGE_CONSERVATION
            species->advancePositionsImplicit( a_dt );
#endif
         }
         else { species->advanceParticlesIteratively( a_emfields, a_dt ); }
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
         setCurrentDensity( a_emfields, a_dt, false );
#endif
         setMassMatrices( a_emfields, a_dt );
         computeJfromMassMatrices( a_emfields );
      }
      else {
         setCurrentDensity( a_emfields, a_dt, false );
#ifdef MASS_MATRIX_TEST
         setMassMatrices( a_emfields, a_dt );
         computeJfromMassMatrices( a_emfields );
#endif
      }

   }
   
   LevelData<EdgeDataBox>& J = m_currentDensity;
   LevelData<NodeFArrayBox>& Jv = m_currentDensity_virtual;
   addInflowJ( J, Jv, a_emfields, a_dt, a_from_emjacobian );   
   addSubOrbitJ( J, Jv, a_emfields, a_dt, a_from_emjacobian );
   
   // apply BCs, do addOp exchange, then divide by corrected Jacobian
   finalizeSettingJ( J, Jv, a_emfields );

}

void PicSpeciesInterface::filterJ( const EMFields&  a_emfields,
                                   const Real       a_time )
{   
    CH_TIME("PicSpeciesInterface::filterJ()");
    if (a_emfields.filterE_inPlane()) {
        SpaceUtils::exchangeEdgeDataBox(m_currentDensity);
        const FieldBC& field_bc = a_emfields.getFieldBC();
        field_bc.applyEdgeBC( m_currentDensity, a_time );
        SpaceUtils::applyBinomialFilter(m_currentDensity);
    }
    if (a_emfields.filterE_virtual()) {
        SpaceUtils::exchangeNodeFArrayBox(m_currentDensity_virtual);
        const FieldBC& field_bc = a_emfields.getFieldBC();
        field_bc.applyNodeBC( m_currentDensity_virtual, a_time );
        SpaceUtils::applyBinomialFilter(m_currentDensity_virtual);
    }
}
 
void PicSpeciesInterface::postNewtonUpdate( const EMFields&  a_emfields,
                                            const Real       a_time,
                                            const Real       a_dt )
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
      species->advancePositionsImplicit( a_dt );
   }
 
}

void PicSpeciesInterface::setMassMatrices( const EMFields&  a_emfields,
                                           const Real       a_dt )
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
      species->transferFastParticles();
      species->accumulateMassMatrices( m_sigma_xx, m_sigma_xy, m_sigma_xz,
#if CH_SPACEDIM==1
                                       m_sigma_yx, m_sigma_yy, m_sigma_yz,
#endif
                                       m_sigma_zx, m_sigma_zy, m_sigma_zz, 
                                       m_J0, m_J0_virtual, a_emfields, a_dt );
   }
   
   // add ghost cells to valid cells
   //LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   //m_J0.exchange( m_J0.interval(), m_mesh.reverseCopier(), addEdgeOp );

   //LDaddNodeOp<NodeFArrayBox> addNodeOp;
   //m_J0_virtual.exchange( m_J0_virtual.interval(), m_mesh.reverseCopier(), addNodeOp );
   
   // save electric field to compute perturbed E during GMRES
   const LevelData<EdgeDataBox>& electricField = a_emfields.getFilteredElectricField();
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

void PicSpeciesInterface::setMassMatricesForPC( EMFields&  a_emfields )
{
   CH_TIME("PicSpeciesInterface::setMassMatricesForPC()");
   if(!m_use_mass_matrices) return; 

   LevelData<EdgeDataBox>& sigma_xx_pc = a_emfields.getSigmaxxPC();
   LevelData<EdgeDataBox>& sigma_xy_pc = a_emfields.getSigmaxyPC();
   LevelData<EdgeDataBox>& sigma_xz_pc = a_emfields.getSigmaxzPC();
#if CH_SPACEDIM==1
   LevelData<NodeFArrayBox>& sigma_yx_pc = a_emfields.getSigmayxPC();
   LevelData<NodeFArrayBox>& sigma_yy_pc = a_emfields.getSigmayyPC();
   LevelData<NodeFArrayBox>& sigma_yz_pc = a_emfields.getSigmayzPC();
#endif                                     
   LevelData<NodeFArrayBox>& sigma_zx_pc = a_emfields.getSigmazxPC();
   LevelData<NodeFArrayBox>& sigma_zy_pc = a_emfields.getSigmazyPC();
   LevelData<NodeFArrayBox>& sigma_zz_pc = a_emfields.getSigmazzPC();
  
   const int pc_mass_matrix_width = a_emfields.getMassMatrixPCwidth();
   const bool include_ij = a_emfields.includeMassMatrixij();
#if CH_SPACEDIM==1
   const bool use_filtering = a_emfields.useFiltering();      
#endif

   const IntVect& ncomp_xx_pc = a_emfields.getNcompxxPC();
   const IntVect& ncomp_xy_pc = a_emfields.getNcompxyPC();
   const IntVect& ncomp_xz_pc = a_emfields.getNcompxzPC();
   const IntVect& ncomp_yx_pc = a_emfields.getNcompyxPC();
   const IntVect& ncomp_yy_pc = a_emfields.getNcompyyPC();
   const IntVect& ncomp_yz_pc = a_emfields.getNcompyzPC();
   const IntVect& ncomp_zx_pc = a_emfields.getNcompzxPC();
   const IntVect& ncomp_zy_pc = a_emfields.getNcompzyPC();
   const IntVect& ncomp_zz_pc = a_emfields.getNcompzzPC();

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {

#if CH_SPACEDIM==1
      // copy sigma_xx
      const FArrayBox& this_sigma_xx = m_sigma_xx[dit][0];
      FArrayBox& this_sigma_xx_pc = sigma_xx_pc[dit][0];
      int src_xx_offset = (m_sigma_xx.nComp()-1)/2 - pc_mass_matrix_width;
      if(use_filtering) src_xx_offset -= 2;
      if(src_xx_offset<0) { src_xx_offset = 0; }
      for (int n=0; n<ncomp_xx_pc[0]; n++) {
         int src_comp = src_xx_offset + n;
         int dst_comp = n;
#ifdef PRINT_COMPS   
	 if(m_print_comps) {
            cout << "copying MM to PC for sigma_xx: src_comp = " << src_comp << ", dst_comp = " << dst_comp << endl;
	 }
#endif
         this_sigma_xx_pc.copy(this_sigma_xx,src_comp,dst_comp,1);
      }
      
      // copy sigma_yy  
      const FArrayBox& this_sigma_yy = m_sigma_yy[dit].getFab();
      FArrayBox& this_sigma_yy_pc = sigma_yy_pc[dit].getFab();
      int src_yy_offset = (m_sigma_yy.nComp()-1)/2 - pc_mass_matrix_width;
      if(src_yy_offset<0) { src_yy_offset = 0; }
      for (int n=0; n<ncomp_yy_pc[0]; n++) {
         int src_comp = src_yy_offset + n;
         int dst_comp = n;
#ifdef PRINT_COMPS   
	 if(m_print_comps) {
            cout << "copying MM to PC for sigma_yy: src_comp = " << src_comp 
		                              << ", dst_comp = " << dst_comp << endl;
	 }
#endif
         this_sigma_yy_pc.copy(this_sigma_yy,src_comp,dst_comp,1);
      }

      // copy sigma_zz  
      const FArrayBox& this_sigma_zz = m_sigma_zz[dit].getFab();
      FArrayBox& this_sigma_zz_pc = sigma_zz_pc[dit].getFab();
      int src_zz_offset = (m_sigma_zz.nComp()-1)/2 - pc_mass_matrix_width;
      if(src_zz_offset<0) { src_zz_offset = 0; }
      for (int n=0; n<ncomp_zz_pc[0]; n++) {
         int src_comp = src_zz_offset + n;
         int dst_comp = n;
#ifdef PRINT_COMPS   
	 if(m_print_comps) {
            cout << "copying MM to PC for sigma_zz: src_comp = " << src_comp 
		                              << ", dst_comp = " << dst_comp << endl;
	 }
#endif
         this_sigma_zz_pc.copy(this_sigma_zz,src_comp,dst_comp,1);
      }
#elif CH_SPACEDIM==2
      // copy sigma_xx
      const FArrayBox& this_sigma_xx = m_sigma_xx[dit][0];
      FArrayBox& this_sigma_xx_pc = sigma_xx_pc[dit][0];
      int src_xx_offset = (m_sigma_xx.nComp()-1)/2 - (m_ncomp_xx[0]+1)*pc_mass_matrix_width;
      if(src_xx_offset<0) { src_xx_offset = 0; }
      for (int m=0; m<ncomp_xx_pc[1]; m++) {
         for (int n=0; n<ncomp_xx_pc[0]; n++) {
            int src_comp = src_xx_offset + m*m_ncomp_xx[0] + n;
            int dst_comp = m*ncomp_xx_pc[0] + n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_xx: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_xx_pc.copy(this_sigma_xx,src_comp,dst_comp,1);
	 }
      }

      // copy sigma_yy (sigma_yy = sigma_xx[1] in 2D)
      const FArrayBox& this_sigma_yy = m_sigma_xx[dit][1];
      FArrayBox& this_sigma_yy_pc = sigma_xx_pc[dit][1];
      int src_yy_offset = (m_sigma_xx.nComp()-1)/2 - (m_ncomp_yy[0]+1)*pc_mass_matrix_width;
      if(src_yy_offset<0) { src_yy_offset = 0; }
      for (int m=0; m<ncomp_yy_pc[1]; m++) {
	 for (int n=0; n<ncomp_yy_pc[0]; n++) {
	    int src_comp = src_yy_offset + m*m_ncomp_yy[0] + n;
	    int dst_comp = m*ncomp_yy_pc[0] + n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_yy: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_yy_pc.copy(this_sigma_yy,src_comp,dst_comp,1);
	 }
      }

      // copy sigma_zz  
      const FArrayBox& this_sigma_zz = m_sigma_zz[dit].getFab();
      FArrayBox& this_sigma_zz_pc = sigma_zz_pc[dit].getFab();
      int src_zz_offset = (m_sigma_zz.nComp()-1)/2 - (m_ncomp_zz[0]+1)*pc_mass_matrix_width;
      if(src_zz_offset<0) { src_zz_offset = 0; }
      for (int m=0; m<ncomp_zz_pc[1]; m++) {
         for (int n=0; n<ncomp_zz_pc[0]; n++) {
            int src_comp = src_zz_offset + m*m_ncomp_zz[0] + n;
            int dst_comp = m*ncomp_zz_pc[0] + n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_zz: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_zz_pc.copy(this_sigma_zz,src_comp,dst_comp,1);
         }
      }
#endif

      if(include_ij) {
#if CH_SPACEDIM==1
         // copy sigma_xy
         const FArrayBox& this_sigma_xy = m_sigma_xy[dit][0];
         FArrayBox& this_sigma_xy_pc = sigma_xy_pc[dit][0];
         int src_xy_offset = m_sigma_xy.nComp()/2 - pc_mass_matrix_width;
         if(use_filtering) src_xy_offset -= 2;
         if(src_xy_offset<0) { src_xy_offset = 0; }
         for (int n=0; n<ncomp_xy_pc[0]; n++) {
            int src_comp = src_xy_offset + n;
            int dst_comp = n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_xy: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_xy_pc.copy(this_sigma_xy,src_comp,dst_comp,1);
         }

	 // copy sigma_xz
         const FArrayBox& this_sigma_xz = m_sigma_xz[dit][0];
         FArrayBox& this_sigma_xz_pc = sigma_xz_pc[dit][0];
         int src_xz_offset = m_sigma_xz.nComp()/2 - pc_mass_matrix_width;
         if(use_filtering) src_xz_offset -= 2;
         if(src_xz_offset<0) { src_xz_offset = 0; }
	 for (int n=0; n<ncomp_xz_pc[0]; n++) {
	    int src_comp = src_xz_offset + n;
	    int dst_comp = n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_xz: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_xz_pc.copy(this_sigma_xz,src_comp,dst_comp,1);
	 }
         
	 // copy sigma_yx  
         const FArrayBox& this_sigma_yx = m_sigma_yx[dit].getFab();
         FArrayBox& this_sigma_yx_pc = sigma_yx_pc[dit].getFab();
         int src_yx_offset = m_sigma_yx.nComp()/2 - pc_mass_matrix_width;
	 if(src_yx_offset<0) { src_yx_offset = 0; }
         for (int n=0; n<ncomp_yx_pc[0]; n++) {
            int src_comp = src_yx_offset + n;
            int dst_comp = n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_yx: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_yx_pc.copy(this_sigma_yx,src_comp,dst_comp,1);
         }
      
	 // copy sigma_yz  
         const FArrayBox& this_sigma_yz = m_sigma_yz[dit].getFab();
         FArrayBox& this_sigma_yz_pc = sigma_yz_pc[dit].getFab();
         int src_yz_offset = (m_sigma_yz.nComp()-1)/2 - pc_mass_matrix_width;
	 if(src_yz_offset<0) { src_yz_offset = 0; }
         for (int n=0; n<ncomp_yz_pc[0]; n++) {
            int src_comp = src_yz_offset + n;
            int dst_comp = n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_yz: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_yz_pc.copy(this_sigma_yz,src_comp,dst_comp,1);
         }

         // copy sigma_zx  
         int src_zx_offset = m_sigma_zx.nComp()/2 - pc_mass_matrix_width;
         if(src_zx_offset<0) { src_zx_offset = 0; }
         const FArrayBox& this_sigma_zx = m_sigma_zx[dit].getFab();
         FArrayBox& this_sigma_zx_pc = sigma_zx_pc[dit].getFab();
         for (int n=0; n<ncomp_zx_pc[0]; n++) {
            int src_comp = src_zx_offset + n;
            int dst_comp = n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_zx: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_zx_pc.copy(this_sigma_zx,src_comp,dst_comp,1);
         }
      
         // copy sigma_zy  
         int src_zy_offset = (m_sigma_zy.nComp()-1)/2 - pc_mass_matrix_width;
         if(src_zy_offset<0) { src_zy_offset = 0; }
         const FArrayBox& this_sigma_zy = m_sigma_zy[dit].getFab();
         FArrayBox& this_sigma_zy_pc = sigma_zy_pc[dit].getFab();
         for (int n=0; n<ncomp_zy_pc[0]; n++) {
            int src_comp = src_zy_offset + n;
            int dst_comp = n;
#ifdef PRINT_COMPS   
	    if(m_print_comps) {
               cout << "copying MM to PC for sigma_zy: src_comp = " << src_comp 
		                                 << ", dst_comp = " << dst_comp << endl;
	    }
#endif
            this_sigma_zy_pc.copy(this_sigma_zy,src_comp,dst_comp,1);
         }
#elif CH_SPACEDIM==2
         // copy sigma_xy
         const FArrayBox& this_sigma_xy = m_sigma_xy[dit][0];
         FArrayBox& this_sigma_xy_pc = sigma_xy_pc[dit][0];
         int src_xy_offset = (m_sigma_xy.nComp()+m_ncomp_xy[0])/2 - (m_ncomp_xy[0]+1)*pc_mass_matrix_width;
         if(src_xy_offset<0) { src_xy_offset = 0; }
	 for (int m=0; m<ncomp_xy_pc[1]; m++) {
	    for (int n=0; n<ncomp_xy_pc[0]; n++) {
	       int src_comp = src_xy_offset + m*m_ncomp_xy[0] + n;
	       int dst_comp = m*ncomp_xy_pc[0] + n;
#ifdef PRINT_COMPS   
	       if(m_print_comps) {
                  cout << "copying MM to PC for sigma_xy: src_comp = " << src_comp 
			                            << ", dst_comp = " << dst_comp << endl;
	       }
#endif
               this_sigma_xy_pc.copy(this_sigma_xy,src_comp,dst_comp,1);
	    }
	 }

	 // copy sigma_xz
         const FArrayBox& this_sigma_xz = m_sigma_xz[dit][0];
         FArrayBox& this_sigma_xz_pc = sigma_xz_pc[dit][0];
         int src_xz_offset = m_sigma_xz.nComp()/2 - (m_ncomp_xz[0]+1)*pc_mass_matrix_width;
         if(src_xz_offset<0) { src_xz_offset = 0; }
	 for (int m=0; m<ncomp_xz_pc[1]; m++) {
	    for (int n=0; n<ncomp_xz_pc[0]; n++) {
	       int src_comp = src_xz_offset + m*m_ncomp_xz[0] + n;
	       int dst_comp = m*ncomp_xz_pc[0] + n;
#ifdef PRINT_COMPS   
	       if(m_print_comps) {
                  cout << "copying MM to PC for sigma_xz: src_comp = " << src_comp 
			                            << ", dst_comp = " << dst_comp << endl;
	       }
#endif
               this_sigma_xz_pc.copy(this_sigma_xz,src_comp,dst_comp,1);
	    }
	 }

         // copy sigma_yx (sigma_yx = sigma_xy[1] in 2D)
         const FArrayBox& this_sigma_yx = m_sigma_xy[dit][1];
         FArrayBox& this_sigma_yx_pc = sigma_xy_pc[dit][1];
         int src_yx_offset = (m_sigma_xy.nComp()+m_ncomp_yx[0])/2 - (m_ncomp_yx[0]+1)*pc_mass_matrix_width;
         if(src_yx_offset<0) { src_yx_offset = 0; }
	 for (int m=0; m<ncomp_yx_pc[1]; m++) {
	    for (int n=0; n<ncomp_yx_pc[0]; n++) {
	       int src_comp = src_yx_offset + m*m_ncomp_yx[0] + n;
	       int dst_comp = m*ncomp_yx_pc[0] + n;
#ifdef PRINT_COMPS   
	       if(m_print_comps) {
                  cout << "copying MM to PC for sigma_yx: src_comp = " << src_comp 
			                            << ", dst_comp = " << dst_comp << endl;
	       }
#endif
               this_sigma_yx_pc.copy(this_sigma_yx,src_comp,dst_comp,1);
	    }
	 }
	 
	 // copy sigma_yz (sigma_yz = sigam_xz[1] in 2D)
         const FArrayBox& this_sigma_yz = m_sigma_xz[dit][1];
         FArrayBox& this_sigma_yz_pc = sigma_xz_pc[dit][1];
         int src_yz_offset = (m_sigma_xz.nComp()+m_ncomp_yz[0]-1)/2 - (m_ncomp_yz[0]+1)*pc_mass_matrix_width;
         if(src_yz_offset<0) { src_yz_offset = 0; }
	 for (int m=0; m<ncomp_yz_pc[1]; m++) {
	    for (int n=0; n<ncomp_yz_pc[0]; n++) {
	       int src_comp = src_yz_offset + m*m_ncomp_yz[0] + n;
	       int dst_comp = m*ncomp_yz_pc[0] + n;
#ifdef PRINT_COMPS   
	       if(m_print_comps) {
                  cout << "copying MM to PC for sigma_yz: src_comp = " << src_comp 
			                            << ", dst_comp = " << dst_comp << endl;
	       }
#endif
               this_sigma_yz_pc.copy(this_sigma_yz,src_comp,dst_comp,1);
	    }
	 }

         // copy sigma_zx  
         const FArrayBox& this_sigma_zx = m_sigma_zx[dit].getFab();
         FArrayBox& this_sigma_zx_pc = sigma_zx_pc[dit].getFab();
         int src_zx_offset = m_sigma_zx.nComp()/2 - (m_ncomp_zx[0]+1)*pc_mass_matrix_width;
         if(src_zx_offset<0) { src_zx_offset = 0; }
         for (int m=0; m<ncomp_zx_pc[1]; m++) {
            for (int n=0; n<ncomp_zx_pc[0]; n++) {
               int src_comp = src_zx_offset + m*m_ncomp_zx[0] + n;
               int dst_comp = m*ncomp_zx_pc[0] + n;
#ifdef PRINT_COMPS   
	       if(m_print_comps) {
                  cout << "copying MM to PC for sigma_zx: src_comp = " << src_comp 
			                            << ", dst_comp = " << dst_comp << endl;
	       }
#endif
               this_sigma_zx_pc.copy(this_sigma_zx,src_comp,dst_comp,1);
            }
         }
      
         // copy sigma_zy  
         const FArrayBox& this_sigma_zy = m_sigma_zy[dit].getFab();
         FArrayBox& this_sigma_zy_pc = sigma_zy_pc[dit].getFab();
         int src_zy_offset = (m_sigma_zy.nComp()+m_ncomp_zy[0]-1)/2 - (m_ncomp_zy[0]+1)*pc_mass_matrix_width;
         if(src_zy_offset<0) { src_zy_offset = 0; }
         for (int m=0; m<ncomp_zy_pc[1]; m++) {
            for (int n=0; n<ncomp_zy_pc[0]; n++) {
               int src_comp = src_zy_offset + m*m_ncomp_zy[0] + n;
               int dst_comp = m*ncomp_zy_pc[0] + n;
#ifdef PRINT_COMPS   
	       if(m_print_comps) {
                  cout << "copying MM to PC for sigma_zy: src_comp = " << src_comp 
			                            << ", dst_comp = " << dst_comp << endl;
	       }
#endif
               this_sigma_zy_pc.copy(this_sigma_zy,src_comp,dst_comp,1);
            }
         }
#endif
      }

   }
#ifdef PRINT_COMPS   
   m_print_comps = false;
   #undef PRINT_COMPS
#endif

   LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   sigma_xx_pc.exchange( sigma_xx_pc.interval(), m_mesh.reverseCopier(), addEdgeOp );
   if(include_ij) {
      sigma_xy_pc.exchange( sigma_xy_pc.interval(), m_mesh.reverseCopier(), addEdgeOp );
      sigma_xz_pc.exchange( sigma_xz_pc.interval(), m_mesh.reverseCopier(), addEdgeOp );
   }

   LDaddNodeOp<NodeFArrayBox> addNodeOp;
#if CH_SPACEDIM==1
   if(include_ij) {
      sigma_yx_pc.exchange( sigma_yx_pc.interval(), m_mesh.reverseCopier(), addNodeOp );
      sigma_yz_pc.exchange( sigma_yz_pc.interval(), m_mesh.reverseCopier(), addNodeOp );
   }
   sigma_yy_pc.exchange( sigma_yy_pc.interval(), m_mesh.reverseCopier(), addNodeOp );
#endif
   if(include_ij) {
      sigma_zx_pc.exchange( sigma_zx_pc.interval(), m_mesh.reverseCopier(), addNodeOp );
      sigma_zy_pc.exchange( sigma_zy_pc.interval(), m_mesh.reverseCopier(), addNodeOp );
   }
   sigma_zz_pc.exchange( sigma_zz_pc.interval(), m_mesh.reverseCopier(), addNodeOp );

   // extra normal exchange call should only be needed when using filtering
   //SpaceUtils::exchangeEdgeDataBox(sigma_xx_pc);
   //SpaceUtils::exchangeEdgeDataBox(sigma_xy_pc);
   //SpaceUtils::exchangeEdgeDataBox(sigma_xz_pc);

}

void PicSpeciesInterface::addInflowJ( LevelData<EdgeDataBox>&    a_J,
                                      LevelData<NodeFArrayBox>&  a_Jv,
                                const EMFields&                  a_emfields,
                                const Real                       a_dt,
                                const bool                       a_from_emjacobian ) 
{
   CH_TIME("PicSpeciesInterface::addInflowJ()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());

   // if using the sub-orbit model for inflow particles, update the inflow particle
   // quantities, compute their current, and add it to the total J

   // loop over all pic species and add the current density from inflow particles 
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); ++sp) {

      const PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if( species->charge()==0 ) { continue; }
      if( !species->suborbit_inflowJ() ) { continue; }

      // update the inflow particle quantities     
      species->advanceInflowParticlesAndSetJ( a_emfields, a_dt, a_from_emjacobian );
      
      // compute inflow J and add it to the total J     
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

void PicSpeciesInterface::addSubOrbitJ( LevelData<EdgeDataBox>&    a_J,
                                        LevelData<NodeFArrayBox>&  a_Jv,
                                  const EMFields&                  a_emfields,
                                  const Real                       a_dt,
                                  const bool                       a_from_emjacobian ) 
{
   CH_TIME("PicSpeciesInterface::addSubOrbitJ()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   // if using the sub-orbit model for bulk particles, loop over particle sub-orbits, 
   // update the particle quantities, and contribute their current for each sub-orbit

#ifdef MASS_MATRIX_TEST
   LevelData<EdgeDataBox>& J_T = m_currentDensity_TEST;
   LevelData<NodeFArrayBox>& Jv_T = m_currentDensity_virtual_TEST;
#endif

   for (int sp=0; sp<m_pic_species_ptr_vect.size(); ++sp) {

      const PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if( species->charge()==0 ) { continue; }
      if( !species->use_suborbit_model() ) { continue; }

      // advance of suborbit particles and setting of J are done in same routine  
      species->advanceSubOrbitParticlesAndSetJ( a_emfields, a_dt, a_from_emjacobian );

      // add suborbit J to the total J 
      const LevelData<EdgeDataBox>& this_suborbitJ = species->getSubOrbitJ();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
            const FArrayBox& suborbitJ_dir = this_suborbitJ[dit][dir];
            a_J[dit][dir].plus(suborbitJ_dir,0,0,1);
#ifdef MASS_MATRIX_TEST
            J_T[dit][dir].plus(suborbitJ_dir,0,0,1);
#endif
         }
      }
#if CH_SPACEDIM<3
      const LevelData<NodeFArrayBox>& this_suborbitJv = species->getSubOrbitJ_virtual();
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& suborbitJv = this_suborbitJv[dit].getFab();
         a_Jv[dit].getFab().plus(suborbitJv,0,0,a_Jv.nComp());
#ifdef MASS_MATRIX_TEST
         Jv_T[dit].getFab().plus(suborbitJv,0,0,a_Jv.nComp());
#endif
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
      if (a_from_scatterDt) {
         species->setMomentumDensityFromBinFab();
         species->setEnergyDensityFromBinFab();
      }
      else if (a_num_coulomb>0 && species->charge()!=0) {
         species->setMomentumDensityFromBinFab();
         species->setEnergyDensityFromBinFab();
      }
   }
   
   if (a_num_coulomb>0) { setDebyeLength(); }

}

void PicSpeciesInterface::setDebyeLength() 
{
   CH_TIME("PicSpeciesInterface::setDebyeLength()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());

   Real mcSq_eV = Constants::ME*Constants::CVAC*Constants::CVAC*Constants::EV_PER_JOULE;
   SpaceUtils::zero( m_DebyeLength );
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {

      PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if (species->charge()==0) { continue; }

      // job of caller to make sure the momentus are precomputed
      const LevelData<FArrayBox>& numDen = species->getNumberDensityFromBinFab();
      const LevelData<FArrayBox>& momDen = species->getMomentumDensityFromBinFab();
      const LevelData<FArrayBox>& eneDen = species->getEnergyDensityFromBinFab();
      const Real Aconst = Constants::EP0/Constants::QE/(species->charge()*species->charge());

      for (DataIterator dit(grids); dit.ok(); ++dit) {

               FArrayBox& this_LDe = m_DebyeLength[dit];
         const FArrayBox& this_numDen = numDen[dit];
         const FArrayBox& this_momDen = momDen[dit];
         const FArrayBox& this_eneDen = eneDen[dit];

         const Box gridBox = grids.get(dit);
         BoxIterator gbit(gridBox);
         for (gbit.begin(); gbit.ok(); ++gbit) {

            const IntVect ig = gbit();

            Real N = this_numDen.get(ig,0); // number density [1/m^3]
            if (N==0.0) { continue; }
            Real R = 1.0/std::cbrt(4.0/3.0*Constants::PI*N); // atomic spacing [m]
            Real rho = N*species->mass();

            Real rhoUx = this_momDen.get(ig,0); // [m/s*1/m^3]
            Real rhoUy = this_momDen.get(ig,1); // [m/s*1/m^3]
            Real rhoUz = this_momDen.get(ig,2); // [m/s*1/m^3]
            const Real meanE = (rhoUx*rhoUx + rhoUy*rhoUy + rhoUz*rhoUz)/rho/2.0;

            // compute the species temperature
            Real rhoE = 0.0;
            for (int dir=0; dir<3; dir++) { rhoE += this_eneDen.get(ig,dir); }
            Real T_eV = 2.0/3.0*(rhoE - meanE)/N*mcSq_eV;
            T_eV = std::max(T_eV, 0.01);

            //if(!procID()) cout << "JRA: numDen = " << N << endl;
            //if(!procID()) cout << "JRA: meanE[ig="<<ig<<"] = " << meanE << endl;
            //if(!procID()) cout << "JRA: T_eV[ig="<<ig<<"]  = " << T_eV << endl;

            // compute square of screening length for this species
            const Real LDe_sq = std::max(Aconst*T_eV/N,R*R); // [m^2]

            // store sum of 1/LDe^2 in LDe container
            Real sumLDe_sq_inv = this_LDe.get(ig,0);
            sumLDe_sq_inv += 1.0/LDe_sq;
            this_LDe.set(ig,0,sumLDe_sq_inv);

         }

      }

   }

   // invert sum of species inverse Debye lengths squared and take square root
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& this_LDe = m_DebyeLength[dit];
      const Box gridBox = grids.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) {
         const IntVect ig = gbit();
         Real LDe_sq_inv = this_LDe.get(ig,0);
         Real LDe = 1.0/std::sqrt(LDe_sq_inv); // net screening length [m]
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

Real 
PicSpeciesInterface::cyclotronDt( const Real  a_wc0 )
{
   // a_wc0 = max(abs(B_T))*qe/me)*time_scale 
   // i.e., the max cyclotron freq for unit charge and mass in code units

   Real cyclotron_dt = DBL_MAX;
   Real max_qom = 0.0;
   for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
      PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
      if(species->charge()==0) continue;
      if(species->numParticles()==0) continue;
      Real this_qom = std::abs(species->charge()/species->mass());
      max_qom = std::max(this_qom,max_qom);
   }

   if(max_qom>0.0) { cyclotron_dt = 1.0/a_wc0/max_qom; }
   
   return cyclotron_dt;

}

#include "NamespaceFooter.H"

