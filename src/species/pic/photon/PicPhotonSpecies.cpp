#include "PicPhotonSpecies.H"
#include <cmath>
#include "PicnicConstants.H"
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "GridFunctionFactory.H"
#include "BinFabFactory.H"

#include "CH_HDF5.H"
#include "ParticleIO.H"

#include "MathUtils.H"
#include "SpaceUtils.H"

#include "NamespaceHeader.H"

PicPhotonSpecies::PicPhotonSpecies( ParmParse&   a_ppspc,
                              const int          a_species,
                              const string&      a_name,
                              const DomainGrid&  a_mesh )
   : Species(a_ppspc,a_species,a_name,a_mesh),
     m_species_bc(NULL),
     m_meshInterp(nullptr)
{

   // photon species must have zero mass and zero charge
   CH_assert(m_mass==0.0);
   CH_assert(m_charge==0);

   if (!procID()) {
      cout << "  name = " << m_name << endl;
      cout << "  number = " << m_species << endl;
      cout << "  mass/me   = " << m_mass << endl;
      cout << "  charge/qe = " << m_charge << endl;
      cout << "  motion = " << m_motion << endl;
      cout << "  forces = " << m_forces << endl;
      cout << "  scatter = " << m_scatter << endl;
      cout << endl;
   }
   
   createMeshInterp();

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
   const int ghosts(m_mesh.ghosts());

   // initialize the member LevelDatas
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   
   m_Nppc.define(grids,1,ghostVect);
   m_density.define(grids,1,ghostVect);
   m_momentumDensity.define(grids,3,ghostVect);
   m_energyDensity.define(grids,1,ghostVect);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_Nppc[dit].setVal(0.0);
      m_density[dit].setVal(0.0);
      m_momentumDensity[dit].setVal(0.0);
      m_energyDensity[dit].setVal(0.0);
   } 

   // each box has to be square with fixedBoxSize length to use ParticleData()
   int fixedBoxSize;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Box thisbox( grids[dit] ); 
      fixedBoxSize = thisbox.size(0);
      for (int dir=1; dir<SpaceDim; ++dir) { CH_assert(thisbox.size(dir)==fixedBoxSize); }
      break;
   } 
   
   // ParticleData<T> behaves similar to LevelData<FArrayBox>
   m_data.define(grids, domain, fixedBoxSize,
                 meshSpacing, meshOrigin);
   
   // In order to initialize a LevelData<BinFab<T>>, first need to define a
   // BinFabFactory. The factor has a "create" function to define "new" pointers
   // to the BinFab that lives at the box level... Not sure if I need to use this
   // or If I can just do what I'm doing below....
   BinFabFactory<PhotonParticlePtr> bfptrFactory(meshSpacing, meshOrigin);
   m_data_binfab_ptr.define(grids, 1, 0*IntVect::Unit, bfptrFactory);

}

PicPhotonSpecies::~PicPhotonSpecies()
{
   if(m_species_bc!=NULL) {
      delete m_species_bc;
      m_species_bc = NULL;
   }
   if(m_meshInterp!=NULL) {
      delete m_meshInterp;
      m_meshInterp = NULL;
   }
}

void PicPhotonSpecies::initialize( const CodeUnits&    a_units,
                                   const Real          a_time,
                                   const std::string&  a_restart_file_name )
{

   if (!a_restart_file_name.empty()) {
      initializeFromRestartFile( a_time, a_restart_file_name );
   }
   else {

      // initialize boundary probes to zero
      for (int dir=0; dir<SpaceDim; dir++) {
         m_MassOut_lo[dir] = 0.0;
         m_MassOut_hi[dir] = 0.0;
         m_MomXOut_lo[dir] = 0.0;
         m_MomXOut_hi[dir] = 0.0;
         m_MomYOut_lo[dir] = 0.0;
         m_MomYOut_hi[dir] = 0.0;
         m_MomZOut_lo[dir] = 0.0;
         m_MomZOut_hi[dir] = 0.0;
         m_EnergyOut_lo[dir] = 0.0;
         m_EnergyOut_hi[dir] = 0.0;
      }

   }

   m_cvac_norm = a_units.CvacNorm();

   // define BCs object for this species
   int verbosity = 1;
   m_species_bc = new PicPhotonSpeciesBC( m_name, m_mesh, a_units, verbosity );

   int totalParticleCount = m_data.numParticles();
   if(!procID()) {
      cout << "Finished initializing photon species " << m_name  << endl;
      cout << "total particles  = " << totalParticleCount << endl << endl;
   }
   
}

void PicPhotonSpecies::initializeFromRestartFile( const Real          a_time,
                                                  const std::string&  a_restart_file_name )
{
   if(!procID()) cout << "Initializing photon species " << m_name  << " from restart file..." << endl;

   HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );
   
   //  read cummulative boundary diagnostic data
   HDF5HeaderData header;
   
   // set the species group
   std::stringstream s;
   s << "species" << m_species; 
   if(!procID()) cout << s.str() << endl;
   const std::string group_name = std::string( s.str() + "_data");
   handle.setGroup(group_name);
   
   header.readFromFile( handle );

   // initialize boundary probes to zero
   const ProblemDomain& domain(m_mesh.getDomain());
   for (int dir=0; dir<SpaceDim; dir++) {
      if (domain.isPeriodic(dir)) { continue; }
      std::string dirstr = std::to_string(dir);
      m_MassOut_lo[dir] = header.m_real["massOut_lo"+dirstr];
      m_MassOut_hi[dir] = header.m_real["massOut_hi"+dirstr];
      m_MomXOut_lo[dir] = header.m_real["momXOut_lo"+dirstr];
      m_MomXOut_hi[dir] = header.m_real["momXOut_hi"+dirstr];
      m_MomYOut_lo[dir] = header.m_real["momYOut_lo"+dirstr];
      m_MomYOut_hi[dir] = header.m_real["momYOut_hi"+dirstr];
      m_MomZOut_lo[dir] = header.m_real["momZOut_lo"+dirstr];
      m_MomZOut_hi[dir] = header.m_real["momZOut_hi"+dirstr];
      m_EnergyOut_lo[dir] = header.m_real["energyOut_lo"+dirstr];
      m_EnergyOut_hi[dir] = header.m_real["energyOut_hi"+dirstr];
   }

   // read in the particle data
   readParticlesFromHDF( handle, m_data, "particles" );
   handle.close();

   m_data.remapOutcast();
   CH_assert(m_data.isClosed());

}

void PicPhotonSpecies::advancePhotons( const Real  a_full_dt,
                                       const LevelData<FArrayBox>&  a_Ne )
{
   if (!m_motion) { return; }
   CH_TIME("PicPhotonSpecies::advancePhotons()");

   PhotonParticle* this_part_ptr = NULL;

   const string geom_type = m_mesh.geomType();

   // particle arrays are stored as {X,Z,Y} for anticyclic (i.e. cyl_RZ)
   const bool anticyclic = m_mesh.anticyclic();
   int dirpth = 1;
   if (anticyclic) { dirpth = 2; }

   // bin the particles up by grid cell
   binTheParticles();

   const Real cnormDt = m_cvac_norm*a_full_dt;
   const Real mcSq_J = Constants::ME*Constants::CVAC*Constants::CVAC;
   const Real wpe_factor = Constants::QE*std::sqrt(1.0/Constants::ME/Constants::EP0);

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for (DataIterator dit(grids); dit.ok(); ++dit) { // loop over boxes

      BinFab<PhotonParticlePtr>& thisBinFab = m_data_binfab_ptr[dit];

      const Box gridBox = grids.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid cells

         const IntVect ig = gbit(); // grid index
         const List<PhotonParticlePtr>& cell_pList = thisBinFab(ig,0);

         // compute plasma frequency [Hz] in this cell
         const Real cell_ne = a_Ne[dit].get(ig,0);
         const Real wpe = wpe_factor*std::sqrt(cell_ne);

         ListIterator<PhotonParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) { // loop over particles

            PhotonParticlePtr& this_particle = lit();
            this_part_ptr = this_particle.getPointer();
            RealVect& xp = this_part_ptr->position();
            std::array<Real,3>& up = this_part_ptr->velocity();

            // compute photon frequency [Hz]
            const Real upmag = std::sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]); // Ephoton/(me*c^2)
            const Real omega = upmag*mcSq_J/Constants::HBAR;

            // compute group velocity factor sqrt(1 - wpe^2/w^2)
            Real vg_fact = 1.0 - wpe*wpe/omega/omega;
            if (vg_fact<=0.0) { continue; }
            else { vg_fact = std::sqrt(vg_fact); }

            // advance photon positions
            for (int dir=0; dir<SpaceDim; dir++) { xp[dir] += up[dir]/upmag*vg_fact*cnormDt; }

            // account for geometry reference frame
            if (geom_type=="cyl_R" || geom_type=="cyl_RZ") {

               // compute new position in CAR
               const Real yp = up[dirpth]/upmag*vg_fact*cnormDt;
               const Real thetap = std::atan2(yp,xp[0]);
               const Real costhp = std::cos(thetap);
               const Real sinthp = std::sin(thetap);

               // convert CAR velocity vector to CYL
               const Real upr  =  up[0]*costhp + up[dirpth]*sinthp;
               const Real upth = -up[0]*sinthp + up[dirpth]*costhp;

               // rebase to theta = phi = 0 (upx = upr, upy = upth, upz = upph)
               up[0] = upr;
               up[dirpth] = upth;
               xp[0] = std::sqrt(xp[0]*xp[0] + yp*yp);

            }
            else if (geom_type=="sph_R") { // 1D spherical

               // compute new position in CAR
               const Real yp = up[1]/upmag*vg_fact*cnormDt;
               const Real zp = up[2]/upmag*vg_fact*cnormDt;
               const Real rp = std::sqrt(xp[0]*xp[0] + yp*yp + zp*zp);
               const Real thetap = std::atan2(yp,xp[0]);
               const Real costhp = std::cos(thetap);
               const Real sinthp = std::sin(thetap);
               Real sinphp = 0.0;
               if (rp>0.0) { sinphp = zp/rp; }
               const Real phip = std::asin(sinphp);
               const Real cosphp = std::cos(phip);

               // convert CAR velocity vector to SPH
               const Real upr  =  cosphp*(costhp*up[0] + sinthp*up[1]) + sinphp*up[2];
               const Real upth =  -sinthp*up[0] + costhp*up[1];
               const Real upph = -sinphp*(costhp*up[0] + sinthp*up[1]) + cosphp*up[2];

               // rebase to theta = phi = 0 (upx = upr, upy = upth, upz = upph)
               up[0] = upr;
               up[1] = upth;
               up[2] = upph;
               xp[0] = rp;

            }

         }

      }

   }

   this_part_ptr = NULL;
   delete this_part_ptr;

}

void PicPhotonSpecies::applyBCs( const bool  a_intermediate_advance,
                                 const Real  a_time )
{
   CH_TIME("PicPhotonSpecies::applyBCs()");

   if (!m_motion) { return; }

   // gather outcast particles
   m_data.gatherOutcast();

   // apply BCs to outcasts that are also out of bounds
   List<PhotonParticle>& outcast_list = m_data.outcast();
   m_species_bc->apply( outcast_list, a_intermediate_advance, a_time );

   // remap the outcasts
   m_data.remapOutcast();
   if (!m_data.isClosed()) {
      List<PhotonParticle>& outcast_list = m_data.outcast();
      cout << "m_data is not closed in PicPhotonSpecies::applyBCs() " << endl;
      cout << "procID() = " << procID() << endl;
      cout << "species = " << m_species << endl;
      cout << "outcast_list.length() = " << outcast_list.length() << endl;
      ListIterator<PhotonParticle> lit(outcast_list);
      for(lit.begin(); lit.ok(); ++lit) {
         cout << "position = " << lit().position() << endl;
         cout << "velocity = " << lit().velocity()[0] << endl;
      }
      exit(EXIT_FAILURE);
   }

}

void PicPhotonSpecies::removeOutflowParticles()
{
   if (!m_motion) { return; }

   // remove the outcasts that are out of bounds
   m_species_bc->removeOutflowParticles();

   if (!m_data.isClosed()) {
      List<PhotonParticle>& outcast_list = m_data.outcast();
      cout << "procID() = " << procID() << endl;
      cout << "outcast_list.length() = " << outcast_list.length() << endl;
      ListIterator<PhotonParticle> lit(outcast_list);
      for(lit.begin(); lit.ok(); ++lit) {
         cout << "position = " << lit().position() << endl;
         cout << "velocity = " << lit().velocity()[0] << endl;
      }
   }
   CH_assert(m_data.isClosed());

}

void PicPhotonSpecies::binTheParticles()
{
   CH_TIME("PicPhotonSpecies::binTheParticles()");
   
   ///////////////////////////////////////////////////////////
   //
   //   fill BinFab container with pointers to particles
   //
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for (dit.begin(); dit.ok(); ++dit) { // loop over boxes
      
      //m_data_binfab_ptr[dit].reBin(); 
      List<PhotonParticle>& pList = m_data[dit].listItems();
      BinFab<PhotonParticlePtr>& thisBinFab_ptr = m_data_binfab_ptr[dit];
      
      // clear the binfab_ptr container
      const Box gridBox = BL.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices
         const IntVect ig = gbit();
         List<PhotonParticlePtr>& cell_pList_ptr = thisBinFab_ptr(ig,0);
         cell_pList_ptr.clear();
      }

      // refill the binfab_ptr container
      ListIterator<PhotonParticle> li(pList);
      CH_XD::List<PhotonParticlePtr> pListPtr;
      for(li.begin(); li.ok(); ++li) {
         PhotonParticlePtr particlePtr(li());
         pListPtr.append(particlePtr);
      }
      thisBinFab_ptr.addItems(pListPtr);
   }
}

void PicPhotonSpecies::setNppc() const
{

   const DisjointBoxLayout& grids = m_Nppc.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_Nppc = m_Nppc[dit];
      box_Nppc.setVal(0.0);
      const ListBox<PhotonParticle>& box_list = m_data[dit];

      PhotonMomentType thisMoment = photon_number;
      m_meshInterp->photonMoment( box_Nppc,
                                  box_list.listItems(),
                                  thisMoment );
   }
   m_Nppc.exchange();

}

void PicPhotonSpecies::setNumberDensity() const
{
   CH_TIME("PicPhotonSpecies::setNumberDensity()");
    
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids = m_density.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_rho = m_density[dit];
      box_rho.setVal(0.0);
      const ListBox<PhotonParticle>& box_list = m_data[dit];
      
      PhotonMomentType thisMoment = photon_density;
      m_meshInterp->photonMoment( box_rho,
                                  box_list.listItems(),
                                  thisMoment );
 
      const FArrayBox& box_Ja = Jacobian[dit];
      box_rho.divide(box_Ja,0,0,1);
      box_rho.divide(volume_scale);

   }
   m_density.exchange();

}

void PicPhotonSpecies::setNumberDensityFromBinFab() const
{
   CH_TIME("PicPhotonSpecies::setNumberDensityFromBinFab()");
   
   PhotonParticle* this_part_ptr = NULL;  
    
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   const Real dV_mapped = m_mesh.getMappedCellVolume();
   const Real volume_scale = m_mesh.getVolumeScale();
   const Real dV_phys = dV_mapped*volume_scale;
   const Real kernal = 1.0/dV_phys;

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_rho = m_density[dit];
      box_rho.setVal(0.0);
      
      const BinFab<PhotonParticlePtr>& thisBinFab = m_data_binfab_ptr[dit];
      
      const Box gridBox = grids.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices

         const IntVect ig = gbit(); // grid index
         const List<PhotonParticlePtr>& cell_pList = thisBinFab(ig,0);
         
         Real rho0 = 0.0;
         ListIterator<PhotonParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) { // loop over particles in this grid cell
            PhotonParticlePtr& this_particle = lit();
            this_part_ptr = this_particle.getPointer();
            const Real& wp = this_part_ptr->weight();
            rho0 += wp;
         }
         rho0 *= kernal;
         box_rho.set(ig,0,rho0);

      }
 
      const FArrayBox& box_Ja = Jacobian[dit];
      box_rho.divide(box_Ja,0,0,1);

   }

   this_part_ptr = NULL;
   delete this_part_ptr;

}

void PicPhotonSpecies::setEnergyDensity() const
{
   CH_TIME("PicPhotonSpecies::setEnergyDensity()");

   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids = m_energyDensity.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_ene = m_energyDensity[dit];
      box_ene.setVal(0.0);

      const ListBox<PhotonParticle>& box_list = m_data[dit];
      PhotonMomentType thisMoment = photon_energy;
      m_meshInterp->photonMoment( box_ene,
                                  box_list.listItems(),
                                  thisMoment );

      const FArrayBox& box_Ja = Jacobian[dit];
      box_ene.divide(box_Ja,0,0,1);
      box_ene.divide(volume_scale);

   }
   m_energyDensity.exchange();

}

void PicPhotonSpecies::setMomentumDensity() const
{
   CH_TIME("PicPhotonSpecies::setMomentumDensity()");

   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids = m_momentumDensity.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_mom = m_momentumDensity[dit];
      box_mom.setVal(0.0);

      const ListBox<PhotonParticle>& box_list = m_data[dit];
      PhotonMomentType thisMoment = photon_momentum;
      m_meshInterp->photonMoment( box_mom,
                                  box_list.listItems(),
                                  thisMoment );

      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<box_mom.nComp(); ++n) { 
         box_mom.divide(box_Ja,0,n,1);
      }
      box_mom.divide(volume_scale);

   }
   m_momentumDensity.exchange();

}

void PicPhotonSpecies::createMeshInterp()
{

   const ProblemDomain& domain(m_mesh.getDomain()); 
   const int ghosts(m_mesh.ghosts());
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshLower(m_mesh.getXmin());
   const RealVect& meshUpper(m_mesh.getXmax());
  
   if (m_meshInterp!=NULL) { delete m_meshInterp; }
   m_meshInterp = static_cast<MeshInterp*> (new MeshInterp( domain.domainBox(),
                                                            ghosts, meshSpacing,
                                                            meshLower, meshUpper ));

}

void PicPhotonSpecies::globalMoments(std::vector<Real>&  a_global_moments) const
{
   CH_TIME("PicPhotonSpecies::globalMoments()");

   Real mass_local = 0.0;
   Real momX_local = 0.0;
   Real momY_local = 0.0;
   Real momZ_local = 0.0;
   Real energy_local = 0.0;

   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for (dit.begin(); dit.ok(); ++dit) {

      const ListBox<PhotonParticle>& box_list = m_data[dit];
      const List<PhotonParticle>& pList = box_list.listItems();

      ListIterator<PhotonParticle> lit(pList);
      for(lit.begin(); lit.ok(); ++lit) {

         const Real& wp = lit().weight();
         mass_local += wp; // units are weight

         // up is normalized momentum: sqrt(|up|^2) = Ephoton/(me*c^2)
         const std::array<Real,3>& up = lit().velocity();
         momX_local += wp*up[0];
         momY_local += wp*up[1];
         momZ_local += wp*up[2];

         energy_local += wp*std::sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);

      }

   }

   // convert to SI units
   momX_local *= Constants::ME*Constants::CVAC; // [kg-m/s]
   momY_local *= Constants::ME*Constants::CVAC; // [kg-m/s]
   momZ_local *= Constants::ME*Constants::CVAC; // [kg-m/s]
   energy_local *= Constants::ME*Constants::CVAC*Constants::CVAC; // [Joules]

   std::vector<Real> local_moments;
   local_moments.push_back(mass_local);
   local_moments.push_back(momX_local);
   local_moments.push_back(momY_local);
   local_moments.push_back(momZ_local);
   local_moments.push_back(energy_local);

   a_global_moments.resize(local_moments.size());
#ifdef CH_MPI
   MPI_Allreduce( local_moments.data(),
                  a_global_moments.data(),
                  local_moments.size(),
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
   a_global_moments = local_moments;
#endif

}

void PicPhotonSpecies::bdryMoments(std::vector<Real>&  a_bdry_moments)
{
   CH_TIME("PicPhotonSpecies::bdryMoments()");

   const RealVect& delta_MassOut_lo = m_species_bc->getDeltaMassOut_lo();
   const RealVect& delta_MassOut_hi = m_species_bc->getDeltaMassOut_hi();
   const RealVect& delta_MomXOut_lo = m_species_bc->getDeltaMomXOut_lo();
   const RealVect& delta_MomXOut_hi = m_species_bc->getDeltaMomXOut_hi();
   const RealVect& delta_MomYOut_lo = m_species_bc->getDeltaMomYOut_lo();
   const RealVect& delta_MomYOut_hi = m_species_bc->getDeltaMomYOut_hi();
   const RealVect& delta_MomZOut_lo = m_species_bc->getDeltaMomZOut_lo();
   const RealVect& delta_MomZOut_hi = m_species_bc->getDeltaMomZOut_hi();
   const RealVect& delta_EnergyOut_lo = m_species_bc->getDeltaEnergyOut_lo();
   const RealVect& delta_EnergyOut_hi = m_species_bc->getDeltaEnergyOut_hi();

   std::vector<Real> local_delta;
   const Real mass_fact = 1.0;
   const Real momentum_fact = Constants::ME*Constants::CVAC;
   const Real energy_fact = momentum_fact*Constants::CVAC;
   for (int dir=0; dir<SpaceDim; ++dir) {
      local_delta.push_back(mass_fact*delta_MassOut_lo[dir]); // [weight]
      local_delta.push_back(mass_fact*delta_MassOut_hi[dir]); // [weight]
      local_delta.push_back(momentum_fact*delta_MomXOut_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomXOut_hi[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomYOut_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomYOut_hi[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomZOut_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomZOut_hi[dir]); // [kg-m/s]
      local_delta.push_back(energy_fact*delta_EnergyOut_lo[dir]); // [Joules]
      local_delta.push_back(energy_fact*delta_EnergyOut_hi[dir]); // [Joules]
   }

   a_bdry_moments.resize(local_delta.size());
#ifdef CH_MPI
   MPI_Allreduce( local_delta.data(),
                  a_bdry_moments.data(),
                  local_delta.size(),
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
   a_bdry_moments = local_delta;
#endif

   // add previous running sums to bdry moments
   for (int dir=0; dir<SpaceDim; ++dir) {
      a_bdry_moments.at(0+10*dir) += m_MassOut_lo[dir];
      a_bdry_moments.at(1+10*dir) += m_MassOut_hi[dir];
      a_bdry_moments.at(2+10*dir) += m_MomXOut_lo[dir];
      a_bdry_moments.at(3+10*dir) += m_MomXOut_hi[dir];
      a_bdry_moments.at(4+10*dir) += m_MomYOut_lo[dir];
      a_bdry_moments.at(5+10*dir) += m_MomYOut_hi[dir];
      a_bdry_moments.at(6+10*dir) += m_MomZOut_lo[dir];
      a_bdry_moments.at(7+10*dir) += m_MomZOut_hi[dir];
      a_bdry_moments.at(8+10*dir) += m_EnergyOut_lo[dir];
      a_bdry_moments.at(9+10*dir) += m_EnergyOut_hi[dir];
   }

   // update totals and reset deltas to zero
   for (int dir=0; dir<SpaceDim; ++dir) {
      m_MassOut_lo[dir] = a_bdry_moments.at(0+10*dir);
      m_MassOut_hi[dir] = a_bdry_moments.at(1+10*dir);
      m_MomXOut_lo[dir] = a_bdry_moments.at(2+10*dir);
      m_MomXOut_hi[dir] = a_bdry_moments.at(3+10*dir);
      m_MomYOut_lo[dir] = a_bdry_moments.at(4+10*dir);
      m_MomYOut_hi[dir] = a_bdry_moments.at(5+10*dir);
      m_MomZOut_lo[dir] = a_bdry_moments.at(6+10*dir);
      m_MomZOut_hi[dir] = a_bdry_moments.at(7+10*dir);
      m_EnergyOut_lo[dir] = a_bdry_moments.at(8+10*dir);
      m_EnergyOut_hi[dir] = a_bdry_moments.at(9+10*dir);
   }
   m_species_bc->zeroDeltas();

}

#include "NamespaceFooter.H"

