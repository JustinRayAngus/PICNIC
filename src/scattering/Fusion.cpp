#include "Fusion.H"
#include "MathUtils.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"
#include "SpaceUtils.H"

#include <iostream>
#include <fstream>

#include "NamespaceHeader.H"


void Fusion::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                         const DomainGrid&           a_mesh )
{
    CH_TIME("Fusion::initialize()");
   
    const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
    const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
    const int num_species = a_pic_species_intf.numSpecies();   

    // get pointer to species 1
    CH_assert(m_sp1<num_species);
    const PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];

    // get pointer to species 2
    CH_assert(m_sp2<num_species);
    const PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];

    // get pointer to species 3
    CH_assert(m_sp3<num_species);
    const PicChargedSpeciesPtr this_species3 = pic_species_ptr_vect[species_map[m_sp3]];

    // get pointer to species 4
    CH_assert(m_sp4<num_species);
    const PicChargedSpeciesPtr this_species4 = pic_species_ptr_vect[species_map[m_sp4]];

    // set the species names
    m_species1_name = this_species1->name();
    m_species2_name = this_species2->name();
    m_species3_name = this_species3->name();
    m_species4_name = this_species4->name();
    m_species3b_name = m_species3_name;
    m_species4b_name = m_species4_name;

    // set the species masses
    m_mass1 = this_species1->mass();          // species 1 mass / me
    m_mass2 = this_species2->mass();          // species 2 mass / me
    m_mass3 = this_species3->mass();          // species 3 mass / me
    m_mass4 = this_species4->mass();          // species 4 mass / me
    m_mass3b = m_mass3;
    m_mass4b = m_mass4;
    m_mu = m_mass1*m_mass2/(m_mass1 + m_mass2); // reduced mass
    
    m_mcSq = Constants::ME*Constants::CVAC*Constants::CVAC/Constants::QE; // 0.511e6 

    // assert species 1 and 2 are what they should be
    if (m_fusion_type==DDa || m_fusion_type==DDb || m_fusion_type==DDab) {
        CH_assert(m_sp1==m_sp2);
        CH_assert(this_species1->charge()==1); // deuterium ion
        CH_assert(m_mass1>3670.0 && m_mass1<3671.0);
    }

    if (m_fusion_type==DDa) {
        CH_assert(this_species3->charge()==1); // T ion
        CH_assert(m_mass3>5496.0 && m_mass3<5497.0);
        CH_assert(this_species4->charge()==1); // proton
        CH_assert(m_mass4>1836.0 && m_mass4<1837.0);
    }
    else if (m_fusion_type==DDb) {
        CH_assert(this_species3->charge()==2); // He3 ion
        CH_assert(m_mass3>5495.0 && m_mass3<5496.0);
        CH_assert(this_species4->charge()==0); // neutron
        CH_assert(m_mass4>1838.0 && m_mass4<1839.0);
    }
    else if (m_fusion_type==DDab) {
        CH_assert(this_species3->charge()==1); // T ion
        CH_assert(m_mass3>5496.0 && m_mass3<5497.0);
        CH_assert(this_species4->charge()==1); // proton
        CH_assert(m_mass4>1836.0 && m_mass4<1837.0);
        //
        const PicChargedSpeciesPtr this_species3b = pic_species_ptr_vect[species_map[m_sp3b]];
        m_species3b_name = this_species3b->name();
        m_mass3b = this_species3b->mass();
        CH_assert(this_species3b->charge()==2); // He3 ion
        CH_assert(m_mass3b>5495.0 && m_mass3b<5496.0);
        //
        const PicChargedSpeciesPtr this_species4b = pic_species_ptr_vect[species_map[m_sp4b]];
        m_species4b_name = this_species4b->name();
        m_mass4b = this_species4b->mass();
        CH_assert(this_species4b->charge()==0); // neutron
        CH_assert(m_mass4b>1838.0 && m_mass4b<1839.0);
    }
    else if (m_fusion_type==DT) {
        // D + T ==> He4(3.5 MeV)  + n(14.1 MeV) 
        CH_assert(this_species1->charge()==1); // deuterium ion
        CH_assert(m_mass1>3670.0 && m_mass1<3671.0);
        CH_assert(this_species2->charge()==1); // tritium ion
        CH_assert(m_mass2>5496.0 && m_mass2<5497.0);
        //
        CH_assert(this_species3->charge()==2); // He4 ion
        CH_assert(m_mass3>7294.0 && m_mass3<7295.0);
        CH_assert(this_species4->charge()==0); // neutron
        CH_assert(m_mass4>1838.0 && m_mass4<1839.0);
    }
    else if (m_fusion_type==DHe3) {
        // D + He3 ==> He4(3.6 MeV)  + p(14.7 MeV) 
        CH_assert(this_species1->charge()==1); // deuterium ion
        CH_assert(m_mass1>3670.0 && m_mass1<3671.0);
        CH_assert(this_species2->charge()==2); // He3 ion
        CH_assert(m_mass2>5495.0 && m_mass2<5496.0);
        //
        CH_assert(this_species3->charge()==2); // He4 ion
        CH_assert(m_mass3>7294.0 && m_mass3<7295.0);
        CH_assert(this_species4->charge()==1); // proton
        CH_assert(m_mass4>1836.0 && m_mass4<1837.0);
    }
    else if (m_fusion_type==TT) {
        // T + T ==> He4(3.79 MeV)  + 2n(3.77 MeV)
        CH_assert(this_species1->charge()==1); // tritium ion
        CH_assert(m_mass1>5496.0 && m_mass1<5497.0);
        CH_assert(this_species2->charge()==1); // tritium ion
        CH_assert(m_mass2>5496.0 && m_mass2<5497.0);
        //
        CH_assert(this_species3->charge()==2); // He4 ion
        CH_assert(m_mass3>7294.0 && m_mass3<7295.0);
        CH_assert(this_species4->charge()==0); // neutron
        CH_assert(m_mass4>1838.0 && m_mass4<1839.0);
        CH_assert(m_sp5<num_species);
        const PicChargedSpeciesPtr this_species5 = pic_species_ptr_vect[species_map[m_sp5]];
        m_species5_name = this_species5->name();
        m_mass5 = this_species5->mass();          // species 5 mass / me
        CH_assert(this_species5->charge()==0); // neutron
        CH_assert(m_mass5>1838.0 && m_mass5<1839.0);
    }

    m_Q = (m_mass1 + m_mass2 - m_mass3 - m_mass4 - m_mass5)*m_mcSq; // reaction energy in eV
    m_Qb = (m_mass1 + m_mass2 - m_mass3b - m_mass4b)*m_mcSq; // reaction energy in eV

    // define the fusion product diagnostic container
    const DisjointBoxLayout& grids(a_mesh.getDBL());
    int n_comps = 1;
    if (m_fusion_type==DDab) { n_comps = 2; };
    m_fusionProducts.define(grids,n_comps,IntVect::Zero);
    SpaceUtils::zero( m_fusionProducts );
 
    if (m_verbosity) { printParameters(); }

}

void Fusion::setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const 
{
   CH_TIME("Fusion::setMeanFreeTime()");

   const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 

   const PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];
   const PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];

   if (!this_species1->scatter() || !this_species2->scatter()) { return; }

   const bool setMoments = false;
   const LevelData<FArrayBox>& numberDensity1 = this_species1->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity1 = this_species1->getEnergyDensity(setMoments);

   const LevelData<FArrayBox>& numberDensity2 = this_species2->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity2 = this_species2->getEnergyDensity(setMoments);

   setInterMFT(numberDensity1,energyDensity1,numberDensity2,energyDensity2);

}

void Fusion::setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                          const LevelData<FArrayBox>&  a_energyDensity1,
                          const LevelData<FArrayBox>&  a_numberDensity2,
                          const LevelData<FArrayBox>&  a_energyDensity2 ) const
{
   CH_TIME("Fusion::setInterMFT()");
  
   Real local_numberDensity1, local_betaSqEff1;
   Real local_numberDensity2, local_betaSqEff2;
   Real g12, sigma;
   Real local_nuMax;
   Real box_nuMax=0.0;
 
   const DisjointBoxLayout& grids = a_numberDensity1.disjointBoxLayout();
   DataIterator ditg(grids);
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity1 = a_numberDensity1[ditg];
      const FArrayBox& this_numberDensity2 = a_numberDensity2[ditg];
      const FArrayBox& this_energyDensity1 = a_energyDensity1[ditg];
      const FArrayBox& this_energyDensity2 = a_energyDensity2[ditg];
     
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         const IntVect ig = gbit(); // grid index
       
         // get local density and temperature and compute local VTeff
         local_numberDensity1 = this_numberDensity1.get(ig,0);
         local_numberDensity2 = this_numberDensity2.get(ig,0);
         if (local_numberDensity1*local_numberDensity2 == 0.0) { continue; }

         local_betaSqEff1 = 0.0;
         local_betaSqEff2 = 0.0;
         for (int dir=0; dir<3; dir++) {
            local_betaSqEff1 += 2.0*this_energyDensity1.get(ig,dir);  
            local_betaSqEff2 += 2.0*this_energyDensity2.get(ig,dir);  
         }
         local_betaSqEff1 /= local_numberDensity1*m_mass1; // |V1/cvac|^2
         local_betaSqEff2 /= local_numberDensity2*m_mass2; // |V2/cvac|^2

         // compute local beta and cross section
         g12 = sqrt(local_betaSqEff1 + local_betaSqEff2);
         Real ratio_b = 0.0;
         sigma = getSigma(ratio_b, g12);

         // compute local local nuMax
         local_nuMax = local_numberDensity2*sigma*g12*Constants::CVAC;

         // update box_nuMax for time step
         box_nuMax = Max(box_nuMax,local_nuMax);
      }
   
   }
   
   Real global_nuMax = box_nuMax;
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]

}

void Fusion::applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                        const DomainGrid&           a_mesh,
                        const Real                  a_dt_sec ) const
{
    CH_TIME("Fusion::applyScattering()");
   
    PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
    const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
      
    PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];
    PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];
    if (!this_species1->scatter()) { return; }
    if (!this_species2->scatter()) { return; }
   
    PicChargedSpeciesPtr this_species3 = pic_species_ptr_vect[species_map[m_sp3]];
    PicChargedSpeciesPtr this_species4 = pic_species_ptr_vect[species_map[m_sp4]];

    if (m_fusion_type==DDa || m_fusion_type==DDb) {
        intraSpeciesFusion( *this_species1,
                            *this_species3, *this_species4, *this_species4,
                            *this_species3, *this_species4, a_mesh, a_dt_sec );
    }
    else if (m_fusion_type==DDab) {
        PicChargedSpeciesPtr this_species3b = pic_species_ptr_vect[species_map[m_sp3b]];
        PicChargedSpeciesPtr this_species4b = pic_species_ptr_vect[species_map[m_sp4b]];
        intraSpeciesFusion( *this_species1,
                            *this_species3, *this_species4, *this_species4,
                            *this_species3b, *this_species4b, a_mesh, a_dt_sec );
    }
    else if (m_fusion_type==TT) {
        PicChargedSpeciesPtr this_species5 = pic_species_ptr_vect[species_map[m_sp5]];
        intraSpeciesFusion( *this_species1,
                            *this_species3, *this_species4, *this_species5,
                            *this_species3, *this_species4, a_mesh, a_dt_sec );
    }
    else {
        interSpeciesFusion( *this_species1, *this_species2, 
                            *this_species3, *this_species4, a_mesh, a_dt_sec );
    }

}

void Fusion::intraSpeciesFusion( PicChargedSpecies&  a_species1,
                                 PicChargedSpecies&  a_species3,
                                 PicChargedSpecies&  a_species4,
                                 PicChargedSpecies&  a_species5,
                                 PicChargedSpecies&  a_species3b,
                                 PicChargedSpecies&  a_species4b,
                           const DomainGrid&  a_mesh,
                           const Real         a_dt_sec ) const
{
    CH_TIME("Fusion::intraSpeciesFusion()");
   
    // predefine some pointers to be used below
    JustinsParticle* part1_ptr = NULL;
    JustinsParticle* part2_ptr = NULL;

    // define references to picSpecies1 (primary species)
    LevelData<BinFab<JustinsParticlePtr>>& data1_binfab_ptr = a_species1.partData_binfab();

    const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
    const Real volume_scale = a_mesh.getVolumeScale();
    const Real cellV = a_mesh.getMappedCellVolume()*volume_scale;
    // cellV*Jacobian = cell volume in SI units
  
    // predefine some variables
    Real g12, g12sq, arg, q12;
   
    Real R, costh, sinth, phi, cosph, sinph, rand_num;

    int numCell, Naa; 
    Real den12, cellV_SI, wpMax, wp34;
    Real KEcm_eV, sigma;
   
    // loop over lists in each cell and test shuffle
    const DisjointBoxLayout& grids = data1_binfab_ptr.disjointBoxLayout();
    DataIterator ditg(grids);

    int verbosity=0; // using this as a verbosity flag
    for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes
      
        // initialize lists needed to create new particles
        List<JustinsParticle> species3_pList, species3b_pList;
        List<JustinsParticle> species4_pList, species4b_pList;
        List<JustinsParticle> species5_pList;

        // get references to the scattering species
        BinFab<JustinsParticlePtr>& thisBinFab1_ptr = data1_binfab_ptr[ditg];

        FArrayBox& this_fusionProducts = m_fusionProducts[ditg];
   
        std::vector<JustinsParticlePtr> vector_part_ptrs;
        const Box gridBox = grids.get(ditg);
        BoxIterator gbit(gridBox);
        for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
            const IntVect ig = gbit(); // grid index
            cellV_SI = cellV*Jacobian[ditg].get(ig,0);
       
            List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
            numCell = cell_pList1.length();
            if (numCell < 2 ) { continue; }
    
            bool odd_NxN = true;
            if (numCell % 2 == 0) { odd_NxN = false; }

            Naa = numCell - 1;

            // copy species 1 iterators to a vector in order to shuffle
            vector_part_ptrs.clear();
            vector_part_ptrs.reserve(numCell);
            ListIterator<JustinsParticlePtr> lit(cell_pList1);
            for (lit.begin(); lit.ok(); ++lit) { vector_part_ptrs.push_back(lit()); }
            std::shuffle(vector_part_ptrs.begin(),vector_part_ptrs.end(),global_rand_gen);
        
            // loop over species 1 vector
            for (int p=0; p<vector_part_ptrs.size(); p++) {

                int p1, p2;
                if (odd_NxN) {
                    p1 = p % 2;
                    if (p==0) { p2 = 1; }
                    else if (p==1) { p2 = 2; }
                    else if (p==2) { p2 = 2; }
                }
                else {
                    p1 = p;
                    p++;
                    p2 = p;
                }

                // get particle data for first (primary) particle
                JustinsParticlePtr& part1 = vector_part_ptrs[p1];
                part1_ptr = part1.getPointer();
                const std::array<Real,3>& up1 = part1_ptr->velocity();
                Real& wp1 = part1_ptr->weight();
            
                // get particle data for second particle
                JustinsParticlePtr& part2 = vector_part_ptrs[p2];
                part2_ptr = part2.getPointer();
                const std::array<Real,3>& up2 = part2_ptr->velocity();
                Real& wp2 = part2_ptr->weight();
                
                wpMax = std::max(wp1,wp2);
                den12 = wpMax*Naa/cellV_SI;
                if (odd_NxN) {
                    den12 = den12/2.0;
                    if (p==2) { odd_NxN = false; }
                    // needed for case where p = 0, 1, or 2 were previously killed 
                    // p = 0: p1 = 0, p2 = 1
                    // p = 1: p1 = 1, p2 = 2
                    // p = 2: p1 = 0, p2 = 2
                    if (part1_ptr->killTag()) { continue; }
                    if (part2_ptr->killTag()) { continue; }
                }

                // compute relative beta and energy
#ifdef RELATIVISTIC_PARTICLES

                // compute gamma for particles 1 and 2
                Real gbp1sq = up1[0]*up1[0] + up1[1]*up1[1] + up1[2]*up1[2];
                Real gbp2sq = up2[0]*up2[0] + up2[1]*up2[1] + up2[2]*up2[2];
                Real gammap1 = std::sqrt(1.0 + gbp1sq);
                Real gammap2 = std::sqrt(1.0 + gbp2sq);
                Real Etot = m_mass1*gammap1 + m_mass2*gammap2;

                // compute center-of-momentum (CoM) energy and proper velocity
                Real px0 = m_mass1*up1[0] + m_mass2*up2[0];
                Real py0 = m_mass1*up1[1] + m_mass2*up2[1];
                Real pz0 = m_mass1*up1[2] + m_mass2*up2[2];
                Real Ecm = std::sqrt(Etot*Etot - px0*px0 - py0*py0 - pz0*pz0);
                Real gammacm = Etot/Ecm;
                std::array<Real,3> ucm; // ucm = gammacm*betacm
                ucm[0] = px0/Ecm;
                ucm[1] = py0/Ecm;
                ucm[2] = pz0/Ecm;

                // transform particles 1 and 2 to CoM frame
                std::array<Real,3> up1st, up2st;
                Real gammap1st, gammap2st;
                ScatteringUtils::LorentzTransform( gammap1st, up1st, gammap1, up1, gammacm, ucm );
                ScatteringUtils::LorentzTransform( gammap2st, up2st, gammap2, up2, gammacm, ucm );

                // compute relative beta in cm frame | beta1st - beta2st|
                std::array<Real,3> betaR;
                g12sq = 0.0;
                Real p1st_sq = 0.0;
                for (int n=0; n<3; n++) {
                    betaR[n] = (up1st[n]/gammap1st - up2st[n]/gammap2st);
                    g12sq = g12sq + betaR[n]*betaR[n];
                    p1st_sq = p1st_sq + up1st[n]*up1st[n];
                }
                g12 = std::sqrt(g12sq);
                p1st_sq = p1st_sq*m_mass1*m_mass1; // square of p1 in cm frame

                // compute KE = Ecm - m1 - m2 in cm frame (using round-off error free method)
                KEcm_eV = m_mcSq*p1st_sq*(1.0/m_mass1/(1.0+gammap1st) + 1.0/m_mass2/(1.0+gammap2st));

#else
                g12sq = 0.0;
                std::array<Real,3> betaR; 
                for (int dir=0; dir<3; dir++) {
                    betaR[dir] = up1[dir] - up2[dir]; 
                    g12sq = g12sq + betaR[dir]*betaR[dir]; 
                }
                g12 = std::sqrt(g12sq);
                KEcm_eV = m_mu*m_mcSq*g12sq/2.0; // relative energy eV
#endif

                // get the cross section
                Real ratio_b = 0.0;
                sigma = getSigma(ratio_b, KEcm_eV);
                if (sigma==0.0) { continue; }

                // determine if the pair collide and then do the collision
                Real fmulti = m_fmulti; // >= 1.0
                arg = fmulti*g12*Constants::CVAC*sigma*den12*a_dt_sec;
#ifdef RELATIVISTIC_PARTICLES
                // See LandL Ch. 2 Sec. 12 on invariant cross section
                // const Real gamma_factor = gammacm; // Wu 2021 Eq 6 (WRONG. Don't use)
                const Real gamma_factor = gammap1st*gammap2st/(gammap1*gammap2); // Perez 2012
                arg *= gamma_factor;
#endif
                if (arg > 1.0) {
                    std::cout << "Notice: arg = " << arg << " in intraSpeciesFusion" << std::endl;
                    std::cout << "        species1 = " << m_species1_name << std::endl;
                    std::cout << "        g12 = " << g12 << std::endl;
                    std::cout << "        den12 = " << den12 << std::endl;
                    std::cout << "        sigma = " << sigma << std::endl;
                    arg = 1.0;
                    fmulti = arg/(g12*Constants::CVAC*sigma*den12*a_dt_sec);
#ifdef RELATIVISTIC_PARTICLES
                    fmulti /= gamma_factor;
#endif
                    if (fmulti < 1.0) { fmulti = 1.0; }
                }
                //q12 = 1.0 - exp(-arg);
                q12 = arg;

                // compute weight for fusion products
                wp34 = std::min(wp1,wp2)/fmulti;

                if (q12 > 0.0) { // update grid-based product density diagnostic
                    const long double new_product = q12*wp34/cellV_SI;
                    if (m_fusion_type==DDab) {
                        const long double new_product_a = (1.0-ratio_b)*new_product;
                        const long double old_product_a = this_fusionProducts.get(ig,0);
                        const long double total_product_a = old_product_a + new_product_a;
                        this_fusionProducts.set(ig,0,total_product_a);
                        //
                        const long double new_product_b = ratio_b*new_product;
                        const long double old_product_b = this_fusionProducts.get(ig,1);
                        const long double total_product_b = old_product_b + new_product_b;
                        this_fusionProducts.set(ig,1,total_product_b);
                    }
                    else {
                        const long double old_product = this_fusionProducts.get(ig,0);
                        const long double total_product = old_product + new_product;
                        this_fusionProducts.set(ig,0,total_product);
                    }
                }

                rand_num = MathUtils::rand();
                if (rand_num<=q12) { // this pair collides
                
                    // determine which channel for type DDab
                    Real mass3 = m_mass3;
                    Real mass4 = m_mass4;
                    Real Q = m_Q;
                    bool channel_b = false;
                    if (m_fusion_type==DDab) {
                        if (MathUtils::rand() <= ratio_b) {
                            channel_b = true;
                            mass3 = m_mass3b;
                            mass4 = m_mass4b;
                            Q = m_Qb;
                        }
                    }

                    // use mass4 = mass4 + mass5 = 2*mass4 for TT kinematics
                    // assuming beta4 = beta5 (both neutrons)
                    if (m_fusion_type==TT) { mass4 += m_mass5; }

                    // compute total lab-frame energy before fusion
                    Real KE_before = 0.0;
#ifdef RELATIVISTIC_PARTICLES
                    KE_before = m_mass1*gbp1sq/(1.0 + gammap1)
                              + m_mass2*gbp2sq/(1.0 + gammap2);
#else
                    for (int n=0; n<3; n++) {
                        KE_before = KE_before + 0.5*m_mass1*up1[n]*up1[n] 
                                              + 0.5*m_mass2*up2[n]*up2[n];
                    }
#endif

                    // isotropic scattering
                    R = MathUtils::rand();
                    costh = 1.0 - 2.0*R;
                    sinth = std::sqrt(1.0 - costh*costh);

                    phi = Constants::TWOPI*MathUtils::rand();
                    cosph = std::cos(phi);
                    sinph = std::sqrt(1.0 - cosph*cosph);

#ifdef RELATIVISTIC_PARTICLES
                    // rotate up1 in CoM frame
                    ScatteringUtils::rotateVelocity(up1st,costh,sinth,cosph,sinph);

                    // compute gammap3 in CoM frame
                    const Real gammap3st = (Ecm*Ecm + mass3*mass3 - mass4*mass4)/(2.0*mass3*Ecm);

                    // compute up3 and up4 in CoM frame 
                    const Real up1st_mag = std::sqrt(up1st[0]*up1st[0] + up1st[1]*up1st[1] + up1st[2]*up1st[2]);
                    const Real up3st_mag = std::sqrt(gammap3st*gammap3st - 1.0);
                    std::array<Real,3> up3st, up4st;
                    up3st[0] = up1st[0]*up3st_mag/up1st_mag;
                    up3st[1] = up1st[1]*up3st_mag/up1st_mag;
                    up3st[2] = up1st[2]*up3st_mag/up1st_mag;
                    //
                    up4st[0] = -mass3/mass4*up3st[0];
                    up4st[1] = -mass3/mass4*up3st[1];
                    up4st[2] = -mass3/mass4*up3st[2];
                    Real gammap4st = std::sqrt(1.0 + up4st[0]*up4st[0] + up4st[1]*up4st[1] + up4st[2]*up4st[2]);

                    // Lorentz transform back to lab frame
                    std::array<Real,3> up3, up4;
                    Real gammap3, gammap4;
                    ScatteringUtils::LorentzTransform( gammap3, up3, gammap3st, up3st, gammacm, ucm, true );
                    ScatteringUtils::LorentzTransform( gammap4, up4, gammap4st, up4st, gammacm, ucm, true );
#else
                    // rotate betaR = beta1 - beta2
                    ScatteringUtils::rotateVelocity(betaR,costh,sinth,cosph,sinph);
                    
                    // compute center-of-mass beta
                    std::array<Real,3> betap_cm;
                    for (int n=0; n<3; n++) {
                        betap_cm[n] = (m_mass1*up1[n]+m_mass2*up2[n])/(m_mass1+m_mass2);
                    }

                    // compute magnitude of up3 in center-of-mass frame
                    const Real Efactor_norm = (KEcm_eV + Q)/m_mcSq;
                    const Real mf1 = 2.0*mass4/mass3/(mass3 + mass4);
                    const Real up3_cm = std::sqrt(mf1*Efactor_norm); 

                    // set velocity of particle 3 in center-of-mass frame
                    Real up3x_cm = up3_cm/g12*betaR[0];
                    Real up3y_cm = up3_cm/g12*betaR[1];
                    Real up3z_cm = up3_cm/g12*betaR[2];
                    
                    // set velocity of particle 4 in center-of-mass frame
                    Real up4x_cm = -mass3/mass4*up3x_cm;
                    Real up4y_cm = -mass3/mass4*up3y_cm;
                    Real up4z_cm = -mass3/mass4*up3z_cm;

                    // convert back to lab frame (with momentum correction needed for non-relativistic)
                    std::array<Real,3> up3, up4;
                    const Real mass_fact = (m_mass1+m_mass2)/(mass3+mass4);
                    up3[0] = up3x_cm + betap_cm[0]*mass_fact;
                    up3[1] = up3y_cm + betap_cm[1]*mass_fact;
                    up3[2] = up3z_cm + betap_cm[2]*mass_fact;

                    up4[0] = up4x_cm + betap_cm[0]*mass_fact;
                    up4[1] = up4y_cm + betap_cm[1]*mass_fact;
                    up4[2] = up4z_cm + betap_cm[2]*mass_fact;
#endif
                    
                    // compute total lab-frame energy after fusion
                    Real KE_after = 0.0;
#ifdef RELATIVISTIC_PARTICLES
                    Real gbp3sq = up3[0]*up3[0] + up3[1]*up3[1] + up3[2]*up3[2];
                    Real gbp4sq = up4[0]*up4[0] + up4[1]*up4[1] + up4[2]*up4[2];
                    KE_after = mass3*gbp3sq/(1.0 + gammap3)
                             + mass4*gbp4sq/(1.0 + gammap4);
#else
                    for (int n=0; n<3; n++) {
                        KE_after = KE_after + 0.5*mass3*up3[n]*up3[n] 
                                            + 0.5*mass4*up4[n]*up4[n];
                    }
#endif
                    const Real deltaE = (KE_after - KE_before)*m_mcSq;
                    //std::cout << "JRA: Q = " << Q << "; deltaE = " << deltaE << std::endl;

                    // update energy gained by fusion 
                    // (not exactly Q for non-relativistic)
                    //m_deltaE_fusion += Q*wp34*Constants::JOULE_PER_EV; // Joules
                    m_deltaE_fusion += deltaE*wp34*Constants::JOULE_PER_EV; // Joules
                    
                    // create the new particles
                    if ((m_fusion_type==DDab && channel_b) || m_fusion_type==DDb) {
                    
                        // create the new helium3 split in two to preserve charge conservation
                        const RealVect& Xp1 = part1_ptr->position();
                        const RealVect& Xp2 = part2_ptr->position();
                        JustinsParticle particle31(0.5*wp34, Xp1, up3);
                        JustinsParticle particle32(0.5*wp34, Xp2, up3);
                    
                        // create the new neutron at centroid location
                        const RealVect Xpn = 0.5*(Xp1 + Xp2);
                        JustinsParticle particle4(wp34, Xpn, up4);
               
                        // need to declare ID to make smoke tests not fail
                        // better fix for future is to add option to remove IDs from outputs
                        uint64_t ID = 0;
                        particle31.setID(ID);
                        particle32.setID(ID);
                        particle4.setID(ID);
 
                        // append new particle to the appropriate lists
                        species3b_pList.append(particle31);
                        species3b_pList.append(particle32);
                        species4b_pList.append(particle4);

                    }
                    else if (m_fusion_type==TT) {

                        // create the new helium4 split in two to preserve charge conservation
                        const RealVect& Xp1 = part1_ptr->position();
                        const RealVect& Xp2 = part2_ptr->position();
                        JustinsParticle particle31(0.5*wp34, Xp1, up3);
                        JustinsParticle particle32(0.5*wp34, Xp2, up3);

                        // create the two new neutron at centroid location
                        const RealVect Xpn = 0.5*(Xp1 + Xp2);
                        JustinsParticle particle4(wp34, Xpn, up4);
                        JustinsParticle particle5(wp34, Xpn, up4);

                        // need to declare ID to make smoke tests not fail
                        // better fix for future is to add option to remove IDs from outputs
                        uint64_t ID = 0;
                        particle31.setID(ID);
                        particle32.setID(ID);
                        particle4.setID(ID);
                        particle5.setID(ID);

                        // append new particle to the appropriate lists
                        species3_pList.append(particle31);
                        species3_pList.append(particle32);
                        species4_pList.append(particle4);
                        species5_pList.append(particle5);

                    }
                    else {

                        // create particle 3 at location of particle 1
                        const RealVect& Xp3 = part1_ptr->position();
                        JustinsParticle particle3(wp34, Xp3, up3);
                    
                        // create particle 4 at location of particle 1
                        const RealVect& Xp4 = part2_ptr->position();
                        JustinsParticle particle4(wp34, Xp4, up4);
               
                        // need to declare ID to make smoke tests not fail
                        // better fix for future is to add option to remove IDs from outputs
                        uint64_t ID = 0;
                        particle3.setID(ID);
                        particle4.setID(ID);
 
                        // append new particle to the appropriate lists
                        species3_pList.append(particle3);
                        species4_pList.append(particle4);

                    }
               
                    // reduce weight of particles 1 and 2, remove particles 
                    // from list and tag to kill if needed
                    wp1 = wp1 - wp34;
                    wp2 = wp2 - wp34;
                    if (wp1<1.0e-8*wp34) {
                        part1_ptr->setKillTag();
                        cell_pList1.remove(part1_ptr);
                        numCell = numCell - 1;
                    }
                    if (wp2<1.0e-8*wp34) {
                        part2_ptr->setKillTag();
                        cell_pList1.remove(part2_ptr);
                        numCell = numCell - 1;
                    }
               
                    // check that there are still particles left in the list      
                    if (numCell==0) { break; }

                }

            } // end loop over collisions for this cell
            verbosity=0;
 
        } // end loop over cells
      
        // add the new particles to the appropriate containers
        ParticleData<JustinsParticle>& pData3 = a_species3.partData();
        ParticleData<JustinsParticle>& pData4 = a_species4.partData();
        pData3[ditg].addItemsDestructive(species3_pList);
        pData4[ditg].addItemsDestructive(species4_pList);
        CH_assert(species3_pList.length()==0);
        CH_assert(species4_pList.length()==0);
        if (m_fusion_type==TT) {
            ParticleData<JustinsParticle>& pData5 = a_species5.partData();
            pData5[ditg].addItemsDestructive(species5_pList);
            CH_assert(species5_pList.length()==0);
        }

        ParticleData<JustinsParticle>& pData3b = a_species3b.partData();
        ParticleData<JustinsParticle>& pData4b = a_species4b.partData();
        pData3b[ditg].addItemsDestructive(species3b_pList);
        pData4b[ditg].addItemsDestructive(species4b_pList);
        CH_assert(species3b_pList.length()==0);
        CH_assert(species4b_pList.length()==0);
      
        // remove particles from species 1 if tagged to kill
        ParticleData<JustinsParticle>& pData1 = a_species1.partData();
        List<JustinsParticle>& pList1 = pData1[ditg].listItems();
        ListIterator<JustinsParticle> lit1(pList1);
        for (lit1.begin(); lit1.ok();) {
            const int& kill = lit1().killTag();
            if (kill) { pList1.remove(lit1); }
            else { ++lit1; }
        }
      
    } // end loop over boxes
   
    part1_ptr = NULL;
    part2_ptr = NULL;
    delete part1_ptr;
    delete part2_ptr;

}
      
void Fusion::interSpeciesFusion( PicChargedSpecies&  a_species1,
                                 PicChargedSpecies&  a_species2, 
                                 PicChargedSpecies&  a_species3, 
                                 PicChargedSpecies&  a_species4, 
                           const DomainGrid&  a_mesh,
                           const Real         a_dt_sec ) const
{
    CH_TIME("Fusion::interSpeciesFusion()");

    // predefine some pointers to be used below
    JustinsParticle* part1_ptr = NULL;
    JustinsParticle* part2_ptr = NULL;

    // define references to colliding species and get the density of each species
    LevelData<BinFab<JustinsParticlePtr>>& data1_binfab_ptr = a_species1.partData_binfab();
    LevelData<BinFab<JustinsParticlePtr>>& data2_binfab_ptr = a_species2.partData_binfab();

    const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
    const Real volume_scale = a_mesh.getVolumeScale();
    const Real cellV = a_mesh.getMappedCellVolume()*volume_scale;
    // cellV*Jacobian = cell volume in SI units
  
    // predefine some variables
    int numCell1, numCell2, Nmin, Nmax;
    Real wpMax, wp34;
    Real den12, cellV_SI;
    Real g12, g12sq, arg, q12;
   
    Real R, costh, sinth, phi, cosph, sinph, rand_num; 
    Real KEcm_eV, sigma;

    // loop over lists in each cell and test shuffle
    const DisjointBoxLayout& grids = data1_binfab_ptr.disjointBoxLayout();
    DataIterator ditg(grids);

    int verbosity=0; // using this as a verbosity flag
    for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes
      
        // initialize lists needed to create new particles
        List<JustinsParticle> species3_pList;
        List<JustinsParticle> species4_pList;

        // get references to the scattering species
        BinFab<JustinsParticlePtr>& thisBinFab1_ptr = data1_binfab_ptr[ditg];
        BinFab<JustinsParticlePtr>& thisBinFab2_ptr = data2_binfab_ptr[ditg];

        FArrayBox& this_fusionProducts = m_fusionProducts[ditg];
     
        std::vector<JustinsParticlePtr> vector_part1_ptrs, vector_part2_ptrs;
        const Box gridBox = grids.get(ditg);
        BoxIterator gbit(gridBox);
        for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
            const IntVect ig = gbit(); // grid index
            cellV_SI = cellV*Jacobian[ditg].get(ig,0);

            List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
            List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
            numCell1 = cell_pList1.length();
            numCell2 = cell_pList2.length();
            if (numCell1*numCell2 < 2) { continue; }
            Nmin = std::min(numCell1,numCell2);
            Nmax = std::max(numCell1,numCell2);
    
            bool Nmin_species1 = false;
            if (Nmin==numCell1) { Nmin_species1 = true; }

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
            
            unsigned int p1, p2; 
            for (int p=0; p<Nmax; p++) {

                if (Nmin_species1) {
                    p1 = p % numCell1;
                    p2 = p;
                }
                else {
                    p1 = p;
                    p2 = p % numCell2;
                }

                // get particle data for particle from species 1
                JustinsParticlePtr& part1 = vector_part1_ptrs[p1];
                part1_ptr = part1.getPointer();
                std::array<Real,3>& up1 = part1_ptr->velocity();
                Real& wp1 = part1_ptr->weight();
            
                // get particle data for particle from species 2
                JustinsParticlePtr& part2 = vector_part2_ptrs[p2];
                part2_ptr = part2.getPointer();
                std::array<Real,3>& up2 = part2_ptr->velocity();
                Real& wp2 = part2_ptr->weight();

                // determine maximum weight of pair and density for collision length
                wpMax = std::max(wp1,wp2);
                den12 = wpMax*Nmin/cellV_SI;

                // compute relative beta and kinetic energy
#ifdef RELATIVISTIC_PARTICLES

                // compute gamma for particles 1 and 2
                Real gbp1sq = up1[0]*up1[0] + up1[1]*up1[1] + up1[2]*up1[2];
                Real gbp2sq = up2[0]*up2[0] + up2[1]*up2[1] + up2[2]*up2[2];
                Real gammap1 = std::sqrt(1.0 + gbp1sq);
                Real gammap2 = std::sqrt(1.0 + gbp2sq);
                Real Etot = m_mass1*gammap1 + m_mass2*gammap2;
                
                // compute center-of-momentum (CoM) energy and proper velocity
                Real px0 = m_mass1*up1[0] + m_mass2*up2[0];
                Real py0 = m_mass1*up1[1] + m_mass2*up2[1];
                Real pz0 = m_mass1*up1[2] + m_mass2*up2[2];
                Real Ecm = std::sqrt(Etot*Etot - px0*px0 - py0*py0 - pz0*pz0);
                Real gammacm = Etot/Ecm;
                std::array<Real,3> ucm; // ucm = gammacm*betacm
                ucm[0] = px0/Ecm;
                ucm[1] = py0/Ecm;
                ucm[2] = pz0/Ecm;
                
                // transform particles 1 and 2 to CoM frame
                std::array<Real,3> up1st, up2st;
                Real gammap1st, gammap2st;
                ScatteringUtils::LorentzTransform( gammap1st, up1st, gammap1, up1, gammacm, ucm );
                ScatteringUtils::LorentzTransform( gammap2st, up2st, gammap2, up2, gammacm, ucm );
                // alternative ways to compute up2st and gammap2st with less flops
                //gammap2st = (Ecm - m_mass1*gammap1st)/m_mass2;
                //for (int n=0; n<3; n++) { up2st[n] = -m_mass1/m_mass2*up1st[n]; }
                //gammap2st = std::sqrt(1.0 + up2st[0]*up2st[0] + up2st[1]*up2st[1] + up2st[2]*up2st[2]);

                // compute relative beta in cm frame | beta1st - beta2st|
                std::array<Real,3> betaR;
                g12sq = 0.0;
                Real p1st_sq = 0.0;
                for (int n=0; n<3; n++) {
                    betaR[n] = (up1st[n]/gammap1st - up2st[n]/gammap2st);
                    g12sq = g12sq + betaR[n]*betaR[n];
                    p1st_sq = p1st_sq + up1st[n]*up1st[n];
                }
                g12 = std::sqrt(g12sq);
                p1st_sq = p1st_sq*m_mass1*m_mass1; // square of p1 in cm frame

                // compute KE = Ecm - m1 - m2 in cm frame (using round-off error free method)
                KEcm_eV = m_mcSq*p1st_sq*(1.0/m_mass1/(1.0+gammap1st) + 1.0/m_mass2/(1.0+gammap2st));

#else
                g12sq = 0.0;
                std::array<Real,3> betaR;
                for (int n=0; n<3; n++) {
                    betaR[n] = up1[n] - up2[n];
                    g12sq = g12sq + betaR[n]*betaR[n]; 
                }
                g12 = std::sqrt(g12sq);
                KEcm_eV = m_mu*m_mcSq*g12sq/2.0; // relative energy eV
#endif
                
                // get the cross section
                Real ratio_b = 0.0;
                sigma = getSigma( ratio_b, KEcm_eV);
                if (sigma==0.0) { continue; }
         
                // determine if the pair collide and then do the collision
                Real fmulti = m_fmulti; // >= 1.0
                arg = fmulti*g12*Constants::CVAC*sigma*den12*a_dt_sec;
#ifdef RELATIVISTIC_PARTICLES
                // See LandL Ch. 2 Sec. 12 on invariant cross section
                // const Real gamma_factor = gammacm; // Wu 2021 Eq 6 (WRONG. Don't use)
                const Real gamma_factor = gammap1st*gammap2st/(gammap1*gammap2); // Perez 2012
                arg *= gamma_factor;
#endif
                if (arg > 1.0) {
                    std::cout << "Notice: arg = " << arg << " in interSpeciesFusion" << std::endl;
                    std::cout << "        species1 = " << m_species1_name << std::endl;
                    std::cout << "        species2 = " << m_species2_name << std::endl;
                    std::cout << "        g12 = " << g12 << std::endl;
                    std::cout << "        den12 = " << den12 << std::endl;
                    std::cout << "        sigma = " << sigma << std::endl;
                    arg = 1.0;
                    fmulti = arg/(g12*Constants::CVAC*sigma*den12*a_dt_sec);
#ifdef RELATIVISTIC_PARTICLES
                    fmulti /= gamma_factor;
#endif
                    if (fmulti < 1.0) { fmulti = 1.0; }
                }
                //q12 = 1.0 - exp(-arg);
                q12 = arg;

                // compute weight for fusion products
                wp34 = std::min(wp1,wp2)/fmulti;

                if (q12 > 0.0) { // update grid-based product density diagnostic
                    const long double new_product = q12*wp34/cellV_SI;
                    const long double old_product = this_fusionProducts.get(ig,0);
                    const long double total_product = old_product + new_product;
                    this_fusionProducts.set(ig,0,total_product);
                }

                rand_num = MathUtils::rand();
                if (rand_num<=q12) { // this pair collides

                    // compute total lab-frame energy before fusion
                    Real KE_before = 0.0;
#ifdef RELATIVISTIC_PARTICLES
                    KE_before = m_mass1*gbp1sq/(1.0 + gammap1) 
                              + m_mass2*gbp2sq/(1.0 + gammap2);
#else
                    for (int n=0; n<3; n++) {
                        KE_before = KE_before + 0.5*m_mass1*up1[n]*up1[n]
                                              + 0.5*m_mass2*up2[n]*up2[n];
                    }
#endif

                    // isotropic scattering
                    R = MathUtils::rand();
                    costh = 1.0 - 2.0*R;
                    sinth = std::sqrt(1.0 - costh*costh);

                    phi = Constants::TWOPI*MathUtils::rand();
                    cosph = std::cos(phi);
                    sinph = std::sqrt(1.0 - cosph*cosph);

#ifdef RELATIVISTIC_PARTICLES
                    
                    // rotate up1 in CoM frame
                    ScatteringUtils::rotateVelocity(up1st,costh,sinth,cosph,sinph);
           
                    // compute gammap3 in CoM frame
                    const Real gammap3st = (Ecm*Ecm + m_mass3*m_mass3 - m_mass4*m_mass4)/(2.0*m_mass3*Ecm);

                    // compute up3 and up4 in CoM frame 
                    const Real up1st_mag = std::sqrt(up1st[0]*up1st[0] + up1st[1]*up1st[1] + up1st[2]*up1st[2]);
                    const Real up3st_mag = std::sqrt(gammap3st*gammap3st - 1.0);
                    std::array<Real,3> up3st, up4st;
                    up3st[0] = up1st[0]*up3st_mag/up1st_mag;
                    up3st[1] = up1st[1]*up3st_mag/up1st_mag;
                    up3st[2] = up1st[2]*up3st_mag/up1st_mag;
                    //
                    up4st[0] = -m_mass3/m_mass4*up3st[0];
                    up4st[1] = -m_mass3/m_mass4*up3st[1];
                    up4st[2] = -m_mass3/m_mass4*up3st[2];
                    Real gammap4st = std::sqrt(1.0 + up4st[0]*up4st[0] + up4st[1]*up4st[1] + up4st[2]*up4st[2]);

                    // Lorentz transform back to lab frame
                    std::array<Real,3> up3, up4;
                    Real gammap3, gammap4;
                    ScatteringUtils::LorentzTransform( gammap3, up3, gammap3st, up3st, gammacm, ucm, true );
                    ScatteringUtils::LorentzTransform( gammap4, up4, gammap4st, up4st, gammacm, ucm, true );
                    // set up4 via momentum conservation
                    //up4[0] = (m_mass1*up1[0] + m_mass2*up2[0] - m_mass3*up3[0])/m_mass4;
                    //up4[1] = (m_mass1*up1[1] + m_mass2*up2[1] - m_mass3*up3[1])/m_mass4;
                    //up4[2] = (m_mass1*up1[2] + m_mass2*up2[2] - m_mass3*up3[2])/m_mass4;
                    
#else
                    // rotate betaR = beta1 - beta2
                    ScatteringUtils::rotateVelocity(betaR,costh,sinth,cosph,sinph);
                    
                    // compute center-of-mass beta
                    std::array<Real,3> betap_cm;
                    for (int n=0; n<3; n++) {
                        betap_cm[n] = (m_mass1*up1[n]+m_mass2*up2[n])/(m_mass1+m_mass2);
                    }

                    // compute magnitude of up3 in center-of-mass frame
                    const Real Efactor_norm = (KEcm_eV + m_Q)/m_mcSq;
                    const Real mf1 = 2.0*m_mass4/m_mass3/(m_mass3 + m_mass4);
                    const Real up3_mag = std::sqrt(mf1*Efactor_norm);

                    // set velocity of particle 3 in center-of-mass frame
                    Real up3x_cm = up3_mag/g12*betaR[0];
                    Real up3y_cm = up3_mag/g12*betaR[1];
                    Real up3z_cm = up3_mag/g12*betaR[2];

                    // set velocity of particle 4 in center-of-mass frame
                    Real up4x_cm = -m_mass3/m_mass4*up3x_cm;
                    Real up4y_cm = -m_mass3/m_mass4*up3y_cm;
                    Real up4z_cm = -m_mass3/m_mass4*up3z_cm;

                    // convert back to lab frame (with momentum correction needed for non-relativistic)
                    std::array<Real,3> up3, up4;
                    const Real mass_fact = (m_mass1+m_mass2)/(m_mass3+m_mass4);
                    up3[0] = up3x_cm + betap_cm[0]*mass_fact;
                    up3[1] = up3y_cm + betap_cm[1]*mass_fact;
                    up3[2] = up3z_cm + betap_cm[2]*mass_fact;

                    up4[0] = up4x_cm + betap_cm[0]*mass_fact;
                    up4[1] = up4y_cm + betap_cm[1]*mass_fact;
                    up4[2] = up4z_cm + betap_cm[2]*mass_fact;
#endif

                    // compute total lab-frame energy after fusion
                    Real KE_after = 0.0;
#ifdef RELATIVISTIC_PARTICLES
                    Real gbp3sq = up3[0]*up3[0] + up3[1]*up3[1] + up3[2]*up3[2];
                    Real gbp4sq = up4[0]*up4[0] + up4[1]*up4[1] + up4[2]*up4[2];
                    KE_after = m_mass3*gbp3sq/(1.0 + gammap3) 
                             + m_mass4*gbp4sq/(1.0 + gammap4);
#else
                    for (int n=0; n<3; n++) {
                        KE_after = KE_after + 0.5*m_mass3*up3[n]*up3[n]
                                            + 0.5*m_mass4*up4[n]*up4[n];
                    }
#endif
                    const Real deltaE = (KE_after - KE_before)*m_mcSq;
                    //std::cout << "JRA: Q = " << m_Q << "; deltaE = " << deltaE << std::endl;
                    //std::cout << "JRA: deltaE - Q = " << deltaE - m_Q << std::endl;

                    // update energy gained by fusion 
                    // (not exactly Q for non-relativistic)
                    //m_deltaE_fusion += Q*wp34*Constants::JOULE_PER_EV; // Joules
                    m_deltaE_fusion += deltaE*wp34*Constants::JOULE_PER_EV; // Joules

                    // create the new particles
                    if (m_fusion_type==DT) {

                        // create the new helium4 split in two to preserve charge conservation
                        const int DT_charge = a_species1.charge() + a_species2.charge();
                        CH_assert(a_species3.charge()==DT_charge);
                        const RealVect& Xp1 = part1_ptr->position();
                        const RealVect& Xp2 = part2_ptr->position();
                        JustinsParticle particle31(0.5*wp34, Xp1, up3);
                        JustinsParticle particle32(0.5*wp34, Xp2, up3);

                        // create the new neutron at centroid location
                        CH_assert(a_species4.charge()==0);
                        const RealVect Xpn = 0.5*(Xp1 + Xp2);
                        JustinsParticle particle4(wp34, Xpn, up4);

                        // need to declare ID to make smoke tests not fail
                        // better fix for future is to add option to remove IDs from outputs
                        uint64_t ID = 0;
                        particle31.setID(ID);
                        particle32.setID(ID);
                        particle4.setID(ID);

                        // append new particle to the appropriate lists
                        species3_pList.append(particle31);
                        species3_pList.append(particle32);
                        species4_pList.append(particle4);

                    }
                    else if (m_fusion_type==DHe3) {

                        // create particle 3 (He4) at location of particle 2 (He2)
                        CH_assert(a_species3.charge()==a_species2.charge());
                        const RealVect& Xp3 = part2_ptr->position();
                        JustinsParticle particle3(wp34, Xp3, up3);

                        // create particle 4 (proton) at location of particle 1 (D)
                        CH_assert(a_species4.charge()==a_species1.charge());
                        const RealVect& Xp4 = part1_ptr->position();
                        JustinsParticle particle4(wp34, Xp4, up4);

                        // need to declare ID to make smoke tests not fail
                        // better fix for future is to add option to remove IDs from outputs
                        uint64_t ID = 0;
                        particle3.setID(ID);
                        particle4.setID(ID);

                        // append new particle to the appropriate lists
                        species3_pList.append(particle3);
                        species4_pList.append(particle4);

                    }

                    // reduce weight of particles 1 and 2, remove particles 
                    // from list and tag to kill if needed
                    wp1 = wp1 - wp34;
                    wp2 = wp2 - wp34;
                    if (wp1<1.0e-8*wp34) {
                        part1_ptr->setKillTag();
                        if (Nmin_species1) {
                            vector_part1_ptrs.erase(vector_part1_ptrs.begin()+p1);
                        }
                        cell_pList1.remove(part1_ptr);
                        numCell1 = numCell1 - 1;
                    }
                    if (wp2<1.0e-8*wp34) {
                        part2_ptr->setKillTag();
                        if (!Nmin_species1) {
                            vector_part2_ptrs.erase(vector_part2_ptrs.begin()+p2);
                        }
                        cell_pList2.remove(part2_ptr);
                        numCell2 = numCell2 - 1;
                    }

                    // check that there are still particles left in the list
                    if (numCell1==0) { break; }
                    if (numCell2==0) { break; }

                }

            } // end loop over collisions for this cell
            verbosity=0;

        } // end loop over cells
      
        // add the new particles to the appropriate containers
        ParticleData<JustinsParticle>& pData3 = a_species3.partData();
        ParticleData<JustinsParticle>& pData4 = a_species4.partData();
        pData3[ditg].addItemsDestructive(species3_pList);
        pData4[ditg].addItemsDestructive(species4_pList);
        CH_assert(species3_pList.length()==0);
        CH_assert(species4_pList.length()==0);
      
        // remove particles from species 1 if tagged to kill
        ParticleData<JustinsParticle>& pData1 = a_species1.partData();
        List<JustinsParticle>& pList1 = pData1[ditg].listItems();
        ListIterator<JustinsParticle> lit1(pList1);
        for (lit1.begin(); lit1.ok();) {
            const int& kill = lit1().killTag();
            if (kill) { pList1.remove(lit1); }
            else { ++lit1; }
        }
      
        // remove particles from species 2 if tagged to kill
        ParticleData<JustinsParticle>& pData2 = a_species2.partData();
        List<JustinsParticle>& pList2 = pData2[ditg].listItems();
        ListIterator<JustinsParticle> lit2(pList2);
        for (lit2.begin(); lit2.ok();) {
            const int& kill = lit2().killTag();
            if (kill) { pList2.remove(lit2); }
            else { ++lit2; }
        }

    } // end loop over boxes
   
    part1_ptr = NULL;
    part2_ptr = NULL;
    delete part1_ptr;
    delete part2_ptr;
   
}

Real Fusion::getSigma( Real& a_ratio_b, const Real  a_Erel_eV ) const
{
   Real sigma = 0.0;
   a_ratio_b = 0.0;
   if      (m_fusion_type==DDa)  { sigma = getDDaSigma(a_Erel_eV); }
   else if (m_fusion_type==DDb)  { sigma = getDDbSigma(a_Erel_eV); }
   else if (m_fusion_type==DDab) { 
       Real sigmaA = 0.0;
       Real sigmaB = 0.0;
       sigmaA = getDDaSigma(a_Erel_eV);
       sigmaB = getDDbSigma(a_Erel_eV);
       sigma = sigmaA + sigmaB;
       if (sigma>0.0) { a_ratio_b = sigmaB/sigma; }
   }
   else if (m_fusion_type==DT)   { sigma = getDTSigma(a_Erel_eV); }
   else if (m_fusion_type==DHe3) { sigma = getDHe3Sigma(a_Erel_eV); }
   else if (m_fusion_type==TT)   { sigma = getTTSigma(a_Erel_eV); }
   return sigma;
}

Real Fusion::getDDaSigma( const Real  a_Erel_eV ) const
{
   CH_TIME("Fusion::getDDaSigma()");
   
   // cross section for D + D -> T(1.01 MeV)   + p(3.02 MeV)
   // See Table IV of Bosch and Hale 1992
   // input energy in the center of mass frame (keV)
   // Note that 1 Barn = 1e-28 m^2
 
   const Real x = a_Erel_eV/1.0e3; // relative energy in keV

   Real sigma = 0.0;
   if (x > 0.5) {

      const Real a1 =  5.5576e4;
      const Real a2 =  2.1054e2;
      const Real a3 = -3.2638e-2;
      const Real a4 =  1.4987e-6;
      const Real a5 =  1.8181e-10;
   
      const Real bg = 31.397;
      const Real sf = a1 + x*(a2 + x*(a3 + x*(a4 + x*a5)));
      sigma = sf / (x * std::exp(bg/std::sqrt(x))); // cross section in mB
      sigma = 0.001*sigma;   // cross section in Barns
      sigma = 1.0e-28*sigma; // cross section m^2

   }
   return sigma;   

} 

Real Fusion::getDDbSigma( const Real  a_Erel_eV ) const
{
   CH_TIME("Fusion::getDDbSigma()");

   // cross section for D + D -> He3(0.82 MeV) + n(2.45 MeV) 
   // See Table IV of Bosch and Hale 1992
   // input energy in the center of mass frame (keV)
   // Note that 1 Barn = 1e-28 m^2
 
   const Real x = a_Erel_eV/1.0e3; // relative energy in keV

   Real sigma = 0.0;
   if (x > 0.5) {

      const Real a1 =  5.3701e4;
      const Real a2 =  3.3027e2;
      const Real a3 = -1.2706e-1;
      const Real a4 =  2.9327e-5;
      const Real a5 = -2.5151e-9;
   
      const Real bg = 31.397;
      const Real sf = a1 + x*(a2 + x*(a3 + x*(a4 + x*a5)));
      sigma = sf / (x * std::exp(bg/std::sqrt(x))); // cross section in mB
      sigma = 0.001*sigma;   // cross section in Barns
      sigma = 1.0e-28*sigma; // cross section m^2

   }
   return sigma;

}


Real Fusion::getDTSigma( const Real  a_Erel_eV ) const
{
   CH_TIME("Fusion::getDTSigma()");
 
   // cross section for D + T -> He4(3.5 MeV) + n(14.1 MeV)
   // See Table IV of Bosch and Hale 1992
   // input energy in the center of mass frame (keV)
   // Note that 1 Barn = 1e-28 m^2
 
   const Real x = a_Erel_eV/1.0e3; // relative energy in keV

   Real sigma = 0.0;
   if (x > 0.5) {
      
      Real a1, a2, a3, a4, b1, b2, b3, b4;

      if ( x < 550.0 ) {
          a1 =  6.927e4;
          a2 =  7.454e8;
          a3 =  2.050e6;
          a4 =  5.2002e4;
          b1 =  6.38e1;
          b2 = -9.95e-1;
          b3 =  6.981e-5;
          b4 =  1.728e-4;
      }
      else { // 550 keV -> 4.7 MeV
          a1 =  -1.4714e6;
          a2 =  0.0;
          a3 =  0.0;
          a4 =  0.0;
          b1 = -8.4127e-3;
          b2 =  4.7983e-6;
          b3 = -1.0748e-9;
          b4 =  8.5184e-14;
      }
   
      const Real bg = 34.3827;
      const Real numer = a1 + x*(a2 + x*(a3 + x*a4));
      const Real denom = 1.0 + x*(b1 + x*(b2 + x*(b3 + x*b4)));
      const Real sf = numer/denom;
      sigma = sf / (x * std::exp(bg/std::sqrt(x))); // cross section in mB
      sigma = 0.001*sigma;   // cross section in Barns
      sigma = 1.0e-28*sigma; // cross section m^2

   }
   return sigma;

}

Real Fusion::getDHe3Sigma( const Real  a_Erel_eV ) const
{
   CH_TIME("Fusion::getDHe3Sigma()");
 
   // cross section for D + He3 -> He4(3.6 MeV)  + p(14.7 MeV) 
   // See Table IV of Bosch and Hale 1992
   // input energy in the center of mass frame (keV)
   // Note that 1 Barn = 1e-28 m^2
 
   const Real x = a_Erel_eV/1.0e3; // relative energy in keV

   Real sigma = 0.0;
   if (x > 0.5) {
 
      Real a1, a2, a3, b1, b2, b3, b4;

      if ( x < 900.0 ) {
          a1 =  5.7501e6;
          a2 =  2.5226e3;
          a3 =  4.5566e1;
          b1 = -3.1995e-3;
          b2 = -8.5530e-6;
          b3 =  5.9014e-8;
          b4 =  0.0;
      }
      else { // 900 keV -> 4.8 MeV
          a1 = -8.3993e5;
          a2 =  0.0;
          a3 =  0.0;
          b1 = -2.6830e-3;
          b2 =  1.1633e-6;
          b3 = -2.1332e-10;
          b4 =  1.4250e-14;
      }  
 
      const Real bg = 68.7508;
      const Real numer = a1 + x*(a2 + x*a3);
      const Real denom = 1.0 + x*(b1 + x*(b2 + x*(b3 + x*b4)));
      const Real sf = numer/denom;
      sigma = sf / (x * std::exp(bg/std::sqrt(x))); // cross section in mB
      sigma = 0.001*sigma;   // cross section in Barns
      sigma = 1.0e-28*sigma; // cross section m^2

   }
   return sigma;
}

Real Fusion::getTTSigma( const Real  a_Erel_eV ) const
{
   CH_TIME("Fusion::getDTTSigma()");

   // cross section for T + T -> He4(3.79 MeV) + 2n(3.77 MeV)
   // Fit not provided in Table IV of Bosch and Hale 1992
   // Fit here obtained from JJVDW
   // input energy in the center of mass frame (keV)
   // Note that 1 Barn = 1e-28 m^2

   const Real x = a_Erel_eV/1.0e3; // relative energy in keV

   Real sigma = 0.0;
   if (x > 0.5) {

      Real a1, a2, a3, a4, a5, b1, b2, b3, b4;

      a1 =  1.55664669e+05;
      a2 = -1.11594261e+01;
      a3 =  1.63727168e-01;
      a4 = -9.37018943e-06;
      a5 =  2.86549948e-08;
      b1 =  2.61124209e-04;
      b2 = -1.89760856e-06;
      b3 =  1.32736774e-09;
      b4 =  1.74931772e-15;

      const Real bg = 38.0031824;
      const Real numer = a1 + x*(a2 + x*(a3 + x*(a4 + x*a5)));
      const Real denom = 1.0 + x*(b1 + x*(b2 + x*(b3 + x*b4)));
      const Real sf = numer/denom;
      sigma = sf / (x * std::exp(bg/std::sqrt(x))); // cross section in mB
      sigma = 0.001*sigma;   // cross section in Barns
      sigma = 1.0e-28*sigma; // cross section m^2

   }
   return sigma;
}

#include "NamespaceFooter.H"

