#include "IBremsstrahlung.H"
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


void IBremsstrahlung::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                                  const DomainGrid&           a_mesh )
{
    CH_TIME("IBremsstrahlung::initialize()");
   
    const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
    const PicPhotonSpeciesPtrVect& photon_species_ptr_vect = a_pic_species_intf.getPhotonPtrVect();
    const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
    const int num_species = a_pic_species_intf.numSpecies();   
    
    // get pointer to species 1 (photon)
    CH_assert(m_sp1<num_species);
    PicPhotonSpeciesPtr this_species1(photon_species_ptr_vect[species_map[m_sp1]]);

    // get pointer to species 2 (electron)
    CH_assert(m_sp1<num_species);
    PicChargedSpeciesPtr this_species2(pic_species_ptr_vect[species_map[m_sp2]]);

    // set the species names
    m_species1_name = this_species1->name();
    m_species2_name = this_species2->name();

    // set the species masses
    m_mass1 = this_species1->mass();          // species 1 mass / me
    m_mass2 = this_species2->mass();          // species 2 mass / me

    CH_assert(m_mass1==0.0); // species 1 must be photon
    CH_assert(m_mass2==1.0); // species 2 must be electron
    
    m_mcSq_eV = Constants::ME*Constants::CVAC*Constants::CVAC/Constants::QE; // 0.511e6

    if (m_verbosity) { printParameters(); }

}

void IBremsstrahlung::setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const 
{
   CH_TIME("IBremsstrahlung::setMeanFreeTime()");
  
   m_scatter_dt = DBL_MAX; // mean free time [s]

}

void IBremsstrahlung::applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                        const DomainGrid&           a_mesh,
                        const Real                  a_dt_sec ) const
{
    CH_TIME("IBremsstrahlung::applyScattering()");
   
    PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
    PicPhotonSpeciesPtrVect& photon_species_ptr_vect = a_pic_species_intf.getPhotonPtrVect();
    const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
      
    PicPhotonSpeciesPtr this_species1(photon_species_ptr_vect[species_map[m_sp1]]);
    PicChargedSpeciesPtr this_species2(pic_species_ptr_vect[species_map[m_sp2]]);
    if (!this_species1->scatter()) { return; }
    if (!this_species2->scatter()) { return; }
   
    // cell-based electron-ion collision rate [Hz]
    const LevelData<FArrayBox>& nuei = a_pic_species_intf.getnuei();
   
    applyIBremsstrahlung( *this_species1, *this_species2,
                          nuei, a_mesh, a_dt_sec );

}

void IBremsstrahlung::applyIBremsstrahlung( PicPhotonSpecies&      a_species1,
                                            PicChargedSpecies&     a_species2,
                                      const LevelData<FArrayBox>&  a_nuei,
                                      const DomainGrid&            a_mesh,
                                      const Real                   a_dt_sec ) const
{
    CH_TIME("IBremsstrahlung::applyIBremsstrahlung()");

    // predefine some pointers to be used below
    PhotonParticle* part1_ptr = NULL;
    JustinsParticle* part2_ptr = NULL;
    JustinsParticle* part2p1_ptr = NULL;

    // define references to colliding species and get the density of each species
    LevelData<BinFab<PhotonParticlePtr>>& data1_binfab_ptr = a_species1.partData_binfab();
    LevelData<BinFab<JustinsParticlePtr>>& data2_binfab_ptr = a_species2.partData_binfab();

    // get cell-wise number density of electrons
    const LevelData<FArrayBox>& numDenE = a_species2.getNumberDensity(true);

    // predefine some variables
    int numCell1, numCell2;

    // loop over lists in each cell and test shuffle
    const DisjointBoxLayout& grids = data1_binfab_ptr.disjointBoxLayout();
    DataIterator ditg(grids);
    for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

        // initialize lists needed to create new particles
        List<PhotonParticle> species1_pList;

        // get references to the scattering species
        BinFab<PhotonParticlePtr>& thisBinFab1_ptr = data1_binfab_ptr[ditg];
        BinFab<JustinsParticlePtr>& thisBinFab2_ptr = data2_binfab_ptr[ditg];

        std::vector<JustinsParticlePtr> vector_part2_ptrs;

        const Box gridBox = grids.get(ditg);
        BoxIterator gbit(gridBox);
        for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells

            const IntVect ig = gbit(); // grid index

            // compute energy associated with plasma frequency
            const Real cell_ne = numDenE[ditg].get(ig,0);
            const Real cell_wpe = Constants::QE*std::sqrt(cell_ne/Constants::ME/Constants::EP0);
            const Real cell_Ewpe_eV = Constants::HBAR*cell_wpe/Constants::QE;

            List<PhotonParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
            List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
            numCell1 = cell_pList1.length();
            numCell2 = cell_pList2.length();
            if (numCell1*numCell2 < 1) { continue; }

            // copy electron species iterators to a vector in order to randomly select
            vector_part2_ptrs.clear();
            vector_part2_ptrs.reserve(numCell2);
            ListIterator<JustinsParticlePtr> litE(cell_pList2);
            for (litE.begin(); litE.ok(); ++litE) { vector_part2_ptrs.push_back(litE()); }
            std::shuffle(vector_part2_ptrs.begin(),vector_part2_ptrs.end(),global_rand_gen);

            const Real cell_nuei = a_nuei[ditg].get(ig,0);

            // loop over photons, adjust weight based on abosrpotion
            // tabulate total momentum and energy absorbed
            Real sum_deltaE_eV = 0.0;
            Real sum_dPx = 0.0;
            Real sum_dPy = 0.0;
            Real sum_dPz = 0.0;
            ListIterator<PhotonParticlePtr> lit(cell_pList1);
            for (lit.begin(); lit.ok(); ++lit) { // loop over photons in this grid cell

               PhotonParticlePtr& this_particle = lit();
               part1_ptr = this_particle.getPointer();
               Real& wp1 = part1_ptr->weight();
               std::array<Real,3>& up1 = part1_ptr->velocity();

               // compute photon energy in [eV]
               const Real Ephoton_eV = std::sqrt(up1[0]*up1[0] + up1[1]*up1[1] + up1[2]*up1[2])*m_mcSq_eV;
               const Real nuIB = std::pow(cell_Ewpe_eV/Ephoton_eV,2)*cell_nuei*0.5;
               if (!procID()) {
                  //cout << "JRA: cell_ne   = " << cell_ne << endl;
                  //cout << "JRA: cell_nuei = " << cell_nuei << endl;
                  //cout << "JRA: nuIB = " << nuIB << endl;
                  //cout << "JRA: cell_Ewpe_eV = " << cell_Ewpe_eV << endl;
                  //cout << "JRA: Ephoton_eV = " << Ephoton_eV << endl;
                  //cout << "JRA: a_dt_sec = " << a_dt_sec << endl;
               }

               // update photon weight based on absorption
               const Real wp0 = wp1;
               Real dw = wp1*std::min(nuIB*a_dt_sec,1.0); // weight to be removed from photon
               wp1 -= dw;

               // set kill tag for photon if weight is too small
               if (wp1<1.0e-10*wp0) {
                  dw = wp0;
                  wp1 = 0.0;
                  part1_ptr->setKillTag();
                  cell_pList1.remove(part1_ptr);
                  numCell1 = numCell1 - 1;
               }

               // update total energy lost
               const Real this_deltaE_eV = dw*Ephoton_eV;
               sum_deltaE_eV += this_deltaE_eV;

               // update total momentum lost (code units)
               const Real this_dPx = dw*up1[0];
               const Real this_dPy = dw*up1[1];
               const Real this_dPz = dw*up1[2];
               sum_dPx += this_dPx;
               sum_dPy += this_dPy;
               sum_dPz += this_dPz;

            }

            // update probe for energy gain via absorption
            m_deltaE_IBremsstrahlung += sum_deltaE_eV*Constants::JOULE_PER_EV; // Joules

            // loop over electrons and compute total weight
            Real wpE_tot = 0.0;
            for (litE.begin(); litE.ok(); ++litE) {
               JustinsParticlePtr& this_particle = litE();
               part2_ptr = this_particle.getPointer();
               const Real wp2 = part2_ptr->weight();
               wpE_tot += wp2;
            }

            // loop over electrons and correct for momentum conservation
            for (litE.begin(); litE.ok(); ++litE) {

               JustinsParticlePtr& this_particle = litE();
               part2_ptr = this_particle.getPointer();
               std::array<Real,3>& up2 = part2_ptr->velocity();

               // compute energy before momentum correction
               Real gbp2sq = up2[0]*up2[0] +up2[1]*up2[1] + up2[2]*up2[2];
               Real gammap2 = std::sqrt(1.0 + gbp2sq);
               const Real KEbefore_norm = m_mass2*gbp2sq/(1.0 + gammap2);

               // adjust up = gammap*betap to restore momentum conservation
               up2[0] = up2[0] + sum_dPx/wpE_tot;
               up2[1] = up2[1] + sum_dPy/wpE_tot;
               up2[2] = up2[2] + sum_dPz/wpE_tot;

               // compute energy after momentum correction
               gbp2sq = up2[0]*up2[0] +up2[1]*up2[1] + up2[2]*up2[2];
               gammap2 = std::sqrt(1.0 + gbp2sq);
               const Real KEafter_norm = m_mass2*gbp2sq/(1.0 + gammap2);

               // update Energy violation
               sum_deltaE_eV += (KEbefore_norm - KEafter_norm)*m_mcSq_eV;

            }

            // loop over electrons and re-distribute energy
            int success = 0;
            long double deltaE = -sum_deltaE_eV/m_mcSq_eV;
            Real Erel_cumulative = 0.0;
            for (int p=0; p<vector_part2_ptrs.size(); p++) {
               int p1 = p;
               p++;
               if (p==vector_part2_ptrs.size()) {
                  p = 0;
               }
               int p2 = p;

               // get particle data for first particle
               JustinsParticlePtr& part1 = vector_part2_ptrs[p1];
               part2_ptr = part1.getPointer();
               std::array<Real,3>& betap1 = part2_ptr->velocity();
               const Real wp1 = part2_ptr->weight();

               // get particle data for second particle
               JustinsParticlePtr& part2 = vector_part2_ptrs[p2];
               part2p1_ptr = part2.getPointer();
               std::array<Real,3>& betap2 = part2p1_ptr->velocity();
               const Real wp2 = part2p1_ptr->weight();

               // zero angle inelastic scatter
               ScatteringUtils::modEnergyPairwise( betap1, betap2, m_mass2*wp1, m_mass2*wp2,
                                                   m_energy_fraction, Erel_cumulative, deltaE );
               if (deltaE==0.0) {
                  success = 1;
                  break;
               }

            }

            // resort to original method to restore energy if method above
            // does not succeed. Momentum conservation will not be maintained.
            if (!success) {

               cout << "Notice: Pair-wise energy absorbing method failed in IBremsstrahlung::apply() " << endl;
               cout << "        Resorting to original method, which is not momentum-conserving." << endl;

               // set dgamma if all electron are adjusted by same amount
               sum_deltaE_eV = -deltaE*m_mcSq_eV;
               const Real dgamma0 = sum_deltaE_eV/m_mcSq_eV/wpE_tot;

               // loop over electrons and absorb energy lost by photons
               //Real KEbefore_eV = 0.0;
               //Real KEafter_eV = 0.0;
               for (litE.begin(); litE.ok(); ++litE) {

                  JustinsParticlePtr& this_particle = litE();
                  part2_ptr = this_particle.getPointer();
                  const Real wp2 = part2_ptr->weight();
                  std::array<Real,3>& up2 = part2_ptr->velocity();

                  // compute initial gamma for electron
                  Real gbp2sq = up2[0]*up2[0] + up2[1]*up2[1] + up2[2]*up2[2];
                  Real gammap2 = std::sqrt(1.0 + gbp2sq);
                  Real gp2m1 = gbp2sq/(1.0+gammap2); // gamma - 1 = KE/mc2
                  //KEbefore_eV += wp2*gp2m1*m_mcSq_eV;

                  // compute dgamma for this electron and update sum
                  Real dgamma = std::max(m_energy_fraction*gp2m1,dgamma0*1.001);
                  bool break_signal = false;
                  Real this_deltaE_eV = wp2*dgamma*m_mcSq_eV;
                  if (sum_deltaE_eV<this_deltaE_eV) {
                     dgamma = sum_deltaE_eV/m_mcSq_eV/wp2;
                     sum_deltaE_eV = 0.0;
                     break_signal = true;
                  }
                  else { sum_deltaE_eV -= this_deltaE_eV; }

                  // update gamma with new energy
                  gammap2 += dgamma;
                  gp2m1 += dgamma;
                  const Real gbp2sq_new = (1.0+gammap2)*gp2m1;

                  // scale electron proper velocity
                  const Real scale_factor = std::sqrt(gbp2sq_new/gbp2sq);
                  up2[0] = up2[0]*scale_factor;
                  up2[1] = up2[1]*scale_factor;
                  up2[2] = up2[2]*scale_factor;

                  // compute electron energy after to verify it is correct
                  //gbp2sq = up2[0]*up2[0] + up2[1]*up2[1] + up2[2]*up2[2];
                  //gammap2 = std::sqrt(1.0 + gbp2sq);
                  //gp2m1 = gbp2sq/(1.0+gammap2); // gamma - 1 = KE/mc2
                  //KEafter_eV += wp2*gp2m1*m_mcSq_eV;

                  if (break_signal) { break; }

               }

            }

        } // end loop over cells

        // remove photons if tagged to kill
        ParticleData<PhotonParticle>& pData1 = a_species1.partData();
        List<PhotonParticle>& pList1 = pData1[ditg].listItems();
        ListIterator<PhotonParticle> lit1(pList1);
        for (lit1.begin(); lit1.ok();) {
            const int& kill = lit1().killTag();
            if (kill) { pList1.remove(lit1); }
            else { ++lit1; }
        }

    } // end loop over boxes

    part1_ptr = NULL;
    part2_ptr = NULL;
    part2p1_ptr = NULL;
    delete part1_ptr;
    delete part2_ptr;
    delete part2p1_ptr;

}

#include "NamespaceFooter.H"

