#include "Bremsstrahlung.H"
#include "MathUtils.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "PhotonParticle.H"
#include "PhotonParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"
#include "SpaceUtils.H"

#include <iostream>
#include <fstream>

#include "NamespaceHeader.H"


void Bremsstrahlung::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                         const DomainGrid&           a_mesh )
{
    CH_TIME("Bremsstrahlung::initialize()");
   
    const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
    const PicPhotonSpeciesPtrVect& photon_species_ptr_vect = a_pic_species_intf.getPhotonPtrVect();
    const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
    const int num_species = a_pic_species_intf.numSpecies();   

    // get pointer to species 1 (electron)
    CH_assert(m_sp1<num_species);
    PicChargedSpeciesPtr this_species1(pic_species_ptr_vect[species_map[m_sp1]]);

    // get pointer to species 2 (ion)
    CH_assert(m_sp2<num_species);
    PicChargedSpeciesPtr this_species2(pic_species_ptr_vect[species_map[m_sp2]]);

    // get pointer to species 3 (photon)
    CH_assert(m_sp3<num_species);
    PicPhotonSpeciesPtr this_species3(photon_species_ptr_vect[species_map[m_sp3]]);

    // set the species names
    m_species1_name = this_species1->name();
    m_species2_name = this_species2->name();
    m_species3_name = this_species3->name();

    // set the species masses
    m_mass1 = this_species1->mass();          // species 1 mass / me
    m_mass2 = this_species2->mass();          // species 2 mass / me
    m_mass3 = this_species3->mass();          // species 3 mass / me

    CH_assert(m_mass1==1.0); // species 1 must be electron
    CH_assert(m_mass3==0.0); // species 3 must be photon
    
    m_mcSq = Constants::ME*Constants::CVAC*Constants::CVAC/Constants::QE; // 0.511e6

    // Copy appropriate kdsigdk array for specified atomic number
    m_kdsigdk.resize(m_koT1_grid.size());
    if (m_bremsstrahlung_type==eH) {
       std::copy(m_kdsigdk_eH.begin(), m_kdsigdk_eH.end(), m_kdsigdk.begin());
    }
    else if (m_bremsstrahlung_type==eHe) {
       std::copy(m_kdsigdk_eHe.begin(), m_kdsigdk_eHe.end(), m_kdsigdk.begin());
    }
    else if (m_bremsstrahlung_type==eB) {
       std::copy(m_kdsigdk_eB.begin(), m_kdsigdk_eB.end(), m_kdsigdk.begin());
    }
    else if (m_bremsstrahlung_type==eC) {
       std::copy(m_kdsigdk_eC.begin(), m_kdsigdk_eC.end(), m_kdsigdk.begin());
    }

    // Convert Seltzer and Berger energy-weighted differential cross section to units of [m^2]
    for (int n=0; n<m_koT1_grid.size(); n++) {
       const Real Z = m_atomic_number;
       const Real gamma = 1.0 + m_KEgrid_eV[n]/m_mcSq;
       const Real betaSq = 1.0 - 1.0/gamma/gamma;
       const Real scale_factor = 1.0e-31*Z*Z/betaSq;
       for (auto& elem : m_kdsigdk[n]) {
          elem *= scale_factor;
       }
    }

    m_k_eV.resize(m_koT1_grid.size());   // photon energy grid [eV]
    m_sigmaC.resize(m_koT1_grid.size()); // cumulative cross section [m^2]

    if (m_verbosity) { printParameters(); }

}

void Bremsstrahlung::setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const 
{
   CH_TIME("Bremsstrahlung::setMeanFreeTime()");
   
   const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 

   PicChargedSpeciesPtr this_species1(pic_species_ptr_vect[species_map[m_sp1]]);
   PicChargedSpeciesPtr this_species2(pic_species_ptr_vect[species_map[m_sp2]]);
   
   if (!this_species1->scatter() || !this_species2->scatter()) { return; }
   
   const bool setMoments = false;
   const LevelData<FArrayBox>& numberDensity1 = this_species1->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity1 = this_species1->getEnergyDensity(setMoments);
   
   const LevelData<FArrayBox>& numberDensity2 = this_species2->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity2 = this_species2->getEnergyDensity(setMoments);
 
   setInterMFT(numberDensity1,energyDensity1,numberDensity2,energyDensity2);

}

void Bremsstrahlung::setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                          const LevelData<FArrayBox>&  a_energyDensity1,
                          const LevelData<FArrayBox>&  a_numberDensity2,
                          const LevelData<FArrayBox>&  a_energyDensity2 ) const
{
   CH_TIME("Bremsstrahlung::setInterMFT()");
  
   m_scatter_dt = DBL_MAX; // mean free time [s]

}

void Bremsstrahlung::applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                        const DomainGrid&           a_mesh,
                        const Real                  a_dt_sec ) const
{
    CH_TIME("Bremsstrahlung::applyScattering()");
   
    PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
    PicPhotonSpeciesPtrVect& photon_species_ptr_vect = a_pic_species_intf.getPhotonPtrVect();
    const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
      
    PicChargedSpeciesPtr this_species1(pic_species_ptr_vect[species_map[m_sp1]]);
    PicChargedSpeciesPtr this_species2(pic_species_ptr_vect[species_map[m_sp2]]);
    if (!this_species1->scatter()) { return; }
    if (!this_species2->scatter()) { return; }
   
    PicPhotonSpeciesPtr this_species3(photon_species_ptr_vect[species_map[m_sp3]]);

    interSpeciesBremsstrahlung( *this_species1, *this_species2,
                                *this_species3, a_mesh, a_dt_sec );

}

void Bremsstrahlung::interSpeciesBremsstrahlung( PicChargedSpecies&  a_species1,
                                 PicChargedSpecies&  a_species2, 
                                 PicPhotonSpecies&  a_species3, 
                           const DomainGrid&  a_mesh,
                           const Real         a_dt_sec ) const
{
    CH_TIME("Bremsstrahlung::interSpeciesBremsstrahlung()");

    // predefine some pointers to be used below
    JustinsParticle* part1_ptr = NULL;
    JustinsParticle* part2_ptr = NULL;

    // define references to colliding species and get the density of each species
    LevelData<BinFab<JustinsParticlePtr>>& data1_binfab_ptr = a_species1.partData_binfab();
    LevelData<BinFab<JustinsParticlePtr>>& data2_binfab_ptr = a_species2.partData_binfab();
   
    // get cell-wise number density of electrons
    const LevelData<FArrayBox>& numDen1 = a_species1.getNumberDensity(true);
   
    const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
    const Real volume_scale = a_mesh.getVolumeScale();
    const Real cellV = a_mesh.getMappedCellVolume()*volume_scale;
    // cellV*Jacobian = cell volume in SI units
  
    // predefine some variables
    int numCell1, numCell2;
    Real wp3;
    Real den12, cellV_SI;
    Real g12, arg, q12;
   
    Real KE_eV, sigma;

    // loop over lists in each cell and test shuffle
    const DisjointBoxLayout& grids = data1_binfab_ptr.disjointBoxLayout();
    DataIterator ditg(grids);
    int verbosity=0; // using this as a verbosity flag
    for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes
      
        // initialize lists needed to create new particles
        List<PhotonParticle> species3_pList;

        // get references to the scattering species
        BinFab<JustinsParticlePtr>& thisBinFab1_ptr = data1_binfab_ptr[ditg];
        BinFab<JustinsParticlePtr>& thisBinFab2_ptr = data2_binfab_ptr[ditg];

        std::vector<JustinsParticlePtr> vector_part1_ptrs, vector_part2_ptrs;
        const Box gridBox = grids.get(ditg);
        BoxIterator gbit(gridBox);
        for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
            const IntVect ig = gbit(); // grid index
            cellV_SI = cellV*Jacobian[ditg].get(ig,0);

            // compute cutoff photon energy based on omega = omega_pe
            const Real cell_ne = numDen1[ditg].get(ig,0);
            const Real cell_wpe = Constants::QE*std::sqrt(cell_ne/Constants::ME/Constants::EP0);
            const Real cell_Ec_eV = Constants::HBAR*cell_wpe/Constants::QE;
 
            List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
            List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
            numCell1 = cell_pList1.length();
            numCell2 = cell_pList2.length();
            if (numCell1*numCell2 < 2) { continue; }

            // copy the iterators to a vector in order to shuffle
            vector_part1_ptrs.clear();
            vector_part1_ptrs.reserve(numCell1);
            ListIterator<JustinsParticlePtr> lit1(cell_pList1);
            for (lit1.begin(); lit1.ok(); ++lit1) { vector_part1_ptrs.push_back(lit1()); }
            std::shuffle(vector_part1_ptrs.begin(),vector_part1_ptrs.end(),global_rand_gen);

            // do the same for species 2
            vector_part2_ptrs.clear();
            vector_part2_ptrs.reserve(numCell2);
            ListIterator<JustinsParticlePtr> lit2(cell_pList2);
            for (lit2.begin(); lit2.ok(); ++lit2) { vector_part2_ptrs.push_back(lit2()); }
            std::shuffle(vector_part2_ptrs.begin(),vector_part2_ptrs.end(),global_rand_gen); 

            // loop over electrons            
            unsigned int p1, p2; 
            for (int p=0; p<numCell1; p++) {

                p1 = p;
                p2 = p % numCell2;

                // get particle data for particle from species 1 (electron)
                JustinsParticlePtr& part1 = vector_part1_ptrs[p1];
                part1_ptr = part1.getPointer();
                std::array<Real,3>& up1 = part1_ptr->velocity();
                const Real& wp1 = part1_ptr->weight();

                // get particle data for particle from species 2 (ion)
                JustinsParticlePtr& part2 = vector_part2_ptrs[p2];
                part2_ptr = part2.getPointer();
                std::array<Real,3>& up2 = part2_ptr->velocity();
                const Real& wp2 = part2_ptr->weight();

                // use ion density for collision length
                den12 = wp2*numCell2/cellV_SI;

                // compute gamma for particles 1 (electron) and 2 (ion)
                Real gbp1sq = up1[0]*up1[0] + up1[1]*up1[1] + up1[2]*up1[2];
                Real gbp2sq = up2[0]*up2[0] + up2[1]*up2[1] + up2[2]*up2[2];
                Real gammap1 = std::sqrt(1.0 + gbp1sq);
                Real gammap2 = std::sqrt(1.0 + gbp2sq);

                // save initial ion params as they are used for frame transforming
                const Real gammap2_save = gammap2;
                const std::array<Real,3> up2_save = up2;

                // transform electron to frame where ion is at rest
                std::array<Real,3> up1st;
                Real gammap1st;
                ScatteringUtils::LorentzTransform( gammap1st, up1st, gammap1, up1, gammap2, up2 );

                // compute electron KE in frame where ion is at rest
                Real up1st_sq = up1st[0]*up1st[0] + up1st[1]*up1st[1] + up1st[2]*up1st[2];
                KE_eV = m_mcSq*up1st_sq/(1.0+gammap1st);
                
                // get the cross section
                Real up1st_mag = std::sqrt(up1st_sq);
                g12 = up1st_mag/gammap1st; // magnitude electron beta
                sigma = getSigma( cell_Ec_eV, g12, KE_eV);
                if (sigma==0.0) { continue; }

                // determine if the pair collide and then do the collision
                Real fmulti = m_fmulti; // >= 1.0
                const Real gamma_factor = gammap1st/(gammap1*gammap2);
                arg = fmulti*g12*Constants::CVAC*sigma*den12*a_dt_sec*gamma_factor;
                if (arg > 1.0) {
                    std::cout << "Notice: arg = " << arg << " in interSpeciesBremsstrahlung" << std::endl;
                    std::cout << "        species1 = " << m_species1_name << std::endl;
                    std::cout << "        species2 = " << m_species2_name << std::endl;
                    std::cout << "        g12 = " << g12 << std::endl;
                    std::cout << "        den12 = " << den12 << std::endl;
                    std::cout << "        sigma = " << sigma << std::endl;
                    arg = 1.0;
                    fmulti = arg/(g12*Constants::CVAC*sigma*den12*a_dt_sec*gamma_factor);
                    if (fmulti < 1.0) { fmulti = 1.0; }
                }
                //q12 = 1.0 - exp(-arg);
                q12 = arg;

                const Real rand_num = MathUtils::rand();
                if (rand_num<=q12) { // this pair collides

                    // compute weight for photon
                    wp3 = wp1/fmulti;

                    // compute lab-frame kinetic energy before emission
                    const Real KE_before = wp1*m_mass1*gbp1sq/(1.0 + gammap1)
                                         + wp2*m_mass2*gbp2sq/(1.0 + gammap2);

                    // photon moves in same direction as electron in IRF
                    // compute unit direction vector
                    std::array<Real,3> dir_vec;
                    dir_vec[0] = up1st[0]/up1st_mag;
                    dir_vec[1] = up1st[1]/up1st_mag;
                    dir_vec[2] = up1st[2]/up1st_mag;

                    // get energy of photon (hbar*omega)
                    int index;
                    const Real rand2 = MathUtils::rand();
                    Real Ephoton_norm = MathUtils::linearInterp(index,m_sigmaC,m_k_eV,rand2); // [eV]
                    Ephoton_norm /= m_mcSq; // normalized units

                    // limit Ephoton to KE in weighted center-of-momentum frame
                    // This is needed for physical solution that conserves both momentum and energy
                    const Real mime = wp2*m_mass2/(wp1*m_mass1);
                    const Real KE_cm_norm = mime*(gammap1st-1.0)/(gammap1st - up1st_mag + mime);
                    if (Ephoton_norm > 0.99*KE_cm_norm/fmulti) { Ephoton_norm = KE_cm_norm; }

                    // compute normalized photon momentum vector in IRF
                    std::array<Real,3> up3st;
                    up3st[0] = Ephoton_norm*dir_vec[0];
                    up3st[1] = Ephoton_norm*dir_vec[1];
                    up3st[2] = Ephoton_norm*dir_vec[2];

                    // compute total photon energy in IRF normalized to weighted-electron rest mass energy
                    const Real Eph_norm_weighted = Ephoton_norm/fmulti;

                    // Compute new energy of electron and ion from conservation
                    // of energy and momentum assuming 1D kinematics
                    // gamma_e2 + mi/me*gamma_i2 = gamma_e1 + mi/me - Ephoton_irf_norm
                    // u_e2 + mi/me*u_i2 = u_e1 - Ephoton_irf_norm
                    const Real A0 = up1st_mag - Eph_norm_weighted;
                    const Real B0 = gammap1st - Eph_norm_weighted + mime;
                    const Real D0 = 2.0*(gammap1st - up1st_mag)*Eph_norm_weighted
                                  - 2.0*(gammap1st - Eph_norm_weighted)*mime - 2.0;
                    const Real a0 = 4.0*( 2.0*(gammap1st - Eph_norm_weighted)*mime
                                        - 2.0*(gammap1st - up1st_mag)*Eph_norm_weighted + 1.0 + mime*mime );
                    const Real b0 = 4.0*A0*D0;
                    const Real c0 = 4.0*B0*B0 - D0*D0;
                    const Real root = b0*b0 - 4.0*a0*c0;
                    CH_assert(root>0.0);
                    const Real up1st_mag_new = (-b0 + std::sqrt(root))/(2.0*a0);
                    const Real up2st_mag_new = (A0 - up1st_mag_new)/mime;

                    // reset electron gamma*beta in IRF with post-scatter magnitude
                    up1st[0] = up1st_mag_new*dir_vec[0];
                    up1st[1] = up1st_mag_new*dir_vec[1];
                    up1st[2] = up1st_mag_new*dir_vec[2];
                    Real gammap1st_prime = std::sqrt(1.0 + up1st[0]*up1st[0] + up1st[1]*up1st[1] + up1st[2]*up1st[2]);

                    // set post-scatter ion gamma*beta in IRF
                    std::array<Real,3> up2st;
                    up2st[0] = up2st_mag_new*dir_vec[0];
                    up2st[1] = up2st_mag_new*dir_vec[1];
                    up2st[2] = up2st_mag_new*dir_vec[2];
                    Real gammap2st_prime = std::sqrt(1.0 + up2st[0]*up2st[0] + up2st[1]*up2st[1] + up2st[2]*up2st[2]);

                    // Lorentz transform electron back to lab frame
                    ScatteringUtils::LorentzTransform( gammap1, up1, gammap1st_prime, up1st, gammap2_save, up2_save, true );

                    // Lorentz transform ion back to lab frame
                    ScatteringUtils::LorentzTransform( gammap2, up2, gammap2st_prime, up2st, gammap2_save, up2_save, true );

                    // Lorentz transform photon back to lab frame
                    Real Ephoton_norm_lab;
                    std::array<Real,3> up3;
                    ScatteringUtils::LorentzTransform( Ephoton_norm_lab, up3, Ephoton_norm, up3st, gammap2_save, up2_save, true );

                    // compute kinetic energy after emission
                    gbp1sq = up1[0]*up1[0] + up1[1]*up1[1] + up1[2]*up1[2];
                    gbp2sq = up2[0]*up2[0] + up2[1]*up2[1] + up2[2]*up2[2];
                    const Real KE_after = wp1*m_mass1*gbp1sq/(1.0 + gammap1)
                                        + wp2*m_mass2*gbp2sq/(1.0 + gammap2);

                    // update energy lost to bremsstrahlung
                    const Real deltaE_eV = (KE_after - KE_before)*m_mcSq;
                    //const Real deltaE_eV = -wp3*Ephoton_norm_lab*m_mcSq;
                    m_deltaE_bremsstrahlung -= deltaE_eV*Constants::JOULE_PER_EV; // Joules

                    // create photon
                    if (m_create_photons) {
                        
                        // create the new photon at location of electron
                        const RealVect& Xp1 = part1_ptr->position();
                        PhotonParticle particle3(wp3, Xp1, up3);

                        // need to declare ID to make smoke tests not fail
                        // better fix for future is to add option to remove IDs from outputs
                        uint64_t ID = 0;
                        particle3.setID(ID);

                        // append new particle to the appropriate lists
                        species3_pList.append(particle3);

                    }

                }

            } // end loop over collisions for this cell
            verbosity=0;

        } // end loop over cells
      
        // add the new photons to the appropriate container
        ParticleData<PhotonParticle>& pData3 = a_species3.partData();
        pData3[ditg].addItemsDestructive(species3_pList);
        CH_assert(species3_pList.length()==0);
      
    } // end loop over boxes

    part1_ptr = NULL;
    part2_ptr = NULL;
    delete part1_ptr;
    delete part2_ptr;
   
}

Real Bremsstrahlung::getSigma( const Real  KEcut_eV,
                               const Real  betarel,
                               const Real  KErel_eV ) const
{
   Real sigma = 0.0;
   if (KErel_eV < 1.0e3) { sigma = 0.0; }
   else { sigma = getSigma_SB(KEcut_eV, betarel, KErel_eV); }
   return sigma;
}

Real Bremsstrahlung::getSigma_SB( const Real  KEcut_eV,
                                  const Real  betarel,
                                  const Real  KErel_eV ) const
{
   CH_TIME("Bremsstrahlung::getSigma_SB()");
   
   // cross section for e + ion(Z) -> e + ion(Z) + photon
   // using Seltzer and Berger tables
   
   if (KErel_eV<=KEcut_eV) { return 0.0; }

   // compute the energy-weighted differential cross section
   std::vector<Real> kdsigdk_1D(m_koT1_grid.size());
   if (KErel_eV>=m_KEgrid_eV.back()) { kdsigdk_1D = m_kdsigdk.back(); }
   else {
      // find index for electron energy interpolation
      const int N = m_KEgrid_eV.size();
      int index = N/2;
      while (KErel_eV < m_KEgrid_eV[index]) { index--; }
      while (KErel_eV > m_KEgrid_eV[index+1]) { index++; }

      // compute interpolation weights for k*dsigma/dk
      const Real w0 = (m_KEgrid_eV[index+1] - KErel_eV)/(m_KEgrid_eV[index+1]-m_KEgrid_eV[index]);
      const Real w1 = 1.0 - w0;

      // compute the energy-weighted differential cross section
      for (size_t i = 0; i < kdsigdk_1D.size(); ++i) {
         kdsigdk_1D[i] = w0*m_kdsigdk.at(index)[i] + w1*m_kdsigdk.at(index+1)[i];
      }
   }

   // find lo-index for koT1 cutoff (will typically be 1)
   const Real koT1_cut = std::max(KEcut_eV/KErel_eV,1.0e-4);
   int i0_cut = 0;
   for (int i=1; i<m_koT1_grid.size(); i++) {
      const Real this_val = m_koT1_grid[i];
      if (this_val>koT1_cut) { break; }
      else { i0_cut = i; }
   }
   if (i0_cut == m_koT1_grid.size()-1) { return 0.0; }

   // adjust kdsigdk_1D for cutoff
   const Real w00 = (m_koT1_grid[i0_cut+1]-koT1_cut)/(m_koT1_grid[i0_cut+1]-m_koT1_grid[i0_cut]);
   const Real w01 = 1.0 - w00;
   const Real val_cut = w00*kdsigdk_1D[i0_cut] + w01*kdsigdk_1D[i0_cut+1];
   kdsigdk_1D[i0_cut] = val_cut;
   if (i0_cut>0) {
     for (int i=0; i<i0_cut; i++) { kdsigdk_1D[i] = 0; }
   }

   // define m_k_eV
   for (int i=0; i<m_k_eV.size(); i++) {
      Real this_k_eV = m_koT1_grid[i]*KErel_eV;
      if (i==i0_cut) { this_k_eV = koT1_cut*KErel_eV; }
      m_k_eV[i] = this_k_eV;
   }

   // create cumulative distribution using trapezoidal rule
   for (auto& elem : m_sigmaC) { elem = 0.0; }
   for (int i=i0_cut+1; i<m_koT1_grid.size(); i++) {
      const Real this_dk = m_k_eV[i] - m_k_eV[i-1];
      const Real this_k = (m_k_eV[i] + m_k_eV[i-1])*0.5;
      // using centered k here mitigates divergece effect for k->0
      const Real this_dsigdk = (kdsigdk_1D[i] + kdsigdk_1D[i-1])*0.5/this_k;
      m_sigmaC[i] = m_sigmaC[i-1] + this_dsigdk*this_dk;
   }
   const Real sigmaTot = m_sigmaC.back(); // total cross section [m^2]

   // normalize cumulative cross section vector
   for (auto& elem : m_sigmaC) { elem /= sigmaTot; }

   if (!procID()) {
      //cout << "i0_cut = " << i0_cut << endl;
      //cout << "KErel_eV = " << KErel_eV << endl;
      //cout << "KEcut_eV = " << KEcut_eV << endl;
      //cout << "sigmaTot = " << sigmaTot << endl;
      //cout << "JRA: kdsigdk_1D.size() = " << kdsigdk_1D.size() << endl;
      //cout << "JRA: m_koT1_grid.size() = " << m_koT1_grid.size() << endl;
      //cout << "JRA: m_KEgrid_eV.size() = " << m_KEgrid_eV.size() << endl;
      //cout << "JRA: m_kdsigdk.size() = " << m_kdsigdk.size() << endl;
      //cout << "JRA: m_kdsigdk.at(1).size() = " << m_kdsigdk.at(1).size() << endl;
      //for (auto& elem : kdsigdk_1D) { cout << "JRA: kdsigdk_1D = " << elem << endl; }
      //for (auto& elem : m_k_eV) { cout << "JRA: k_eV = " << elem << endl; }
      //for (auto& elem : m_sigmaC) { cout << "JRA: sigmaC = " << elem << endl; }
   }

   return sigmaTot;   

} 

#include "NamespaceFooter.H"

