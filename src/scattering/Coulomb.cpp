
#include "Coulomb.H"
#include "MathUtils.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"
#include <limits>

#include "NamespaceHeader.H"


void Coulomb::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                          const DomainGrid&           a_mesh )
{
   CH_TIME("Coulomb::initialize()");
   
   const PicSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getPtrVect();

   m_count.define(a_mesh.getDBL(),1,IntVect::Zero);

   // get pointer to species 1 and assert collisions allowed
   CH_assert(m_sp1<pic_species_ptr_vect.size());
   PicSpeciesPtr this_species1(pic_species_ptr_vect[m_sp1]);
      
   // get pointer to species 2 and assert collisions allowed
   CH_assert(m_sp2<pic_species_ptr_vect.size());
   PicSpeciesPtr this_species2(pic_species_ptr_vect[m_sp2]);

   // set the species names
   m_species1_name = this_species1->name();
   m_species2_name = this_species2->name();
   
   // set the species charges and assert they are not neutrals
   m_charge1 = this_species1->charge();  // species 1 charge / |qe|
   m_charge2 = this_species2->charge();  // species 2 charge / |qe|
   CH_assert(m_charge1!=0);
   CH_assert(m_charge2!=0);

   // set the species masses
   m_mass1 = this_species1->mass();    // species 1 mass / me
   m_mass2 = this_species2->mass();    // species 2 mass / me
   m_mu = m_mass1*m_mass2/(m_mass1 + m_mass2); // reduced mass
  
   Real qocSq = Constants::QE*Constants::QE/(Constants::CVAC*Constants::CVAC);
   Real b90_codeToPhys = qocSq/(Constants::TWOPI*Constants::EP0*Constants::ME); 
   m_b90_fact = abs(m_charge1*m_charge2)*b90_codeToPhys;
   m_bqm_fact = Constants::HBAR/(2.0*Constants::ME*Constants::CVAC);

   const Real pi = Constants::PI;
   const Real me = Constants::ME;
   const Real hbar = Constants::HBAR;
   const Real cvac = Constants::CVAC;

   if (m_mass1==1.0 || m_mass2==1.0) { // only for electrons

      // set coefficient for Fermi energy calc
      // EF [Joules] = hbar^2/(2*mass)*(3*pi^2*n)^2/3
      m_EF_fact = hbar*hbar/(2.0*me*m_mu)*std::pow(3.0*pi*pi,2.0/3.0)/(me*cvac*cvac);

      // don't include large-angle scattering for electrons
      if (m_exclude_electron_fas &&
         (m_angular_scattering==NANBU_FAS || m_angular_scattering==NANBU_FAS_v2)) {
         m_angular_scattering = NANBU;
      }

   }

   if (m_verbosity) { printParameters(); }

}

void Coulomb::setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const 
{
   CH_TIME("Coulomb::setMeanFreeTime()");
   
   const LevelData<FArrayBox>& LDe_m = a_pic_species_intf.getDebyeLength();
   
   const PicSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getPtrVect();
   PicSpeciesPtr this_species1(pic_species_ptr_vect[m_sp1]);
   PicSpeciesPtr this_species2(pic_species_ptr_vect[m_sp2]);
   
   if (!this_species1->scatter() || !this_species2->scatter()) return;
   
   const bool setMoments = false;
   const LevelData<FArrayBox>& numDen1 = this_species1->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& momDen1 = this_species1->getMomentumDensity(setMoments);
   const LevelData<FArrayBox>& eneDen1 = this_species1->getEnergyDensity(setMoments);
   
   if (m_sp1==m_sp2) { setIntraMFT(LDe_m,numDen1,momDen1,eneDen1); }
   else {
      const LevelData<FArrayBox>& numDen2 = this_species2->getNumberDensity(setMoments);
      const LevelData<FArrayBox>& momDen2 = this_species2->getMomentumDensity(setMoments);
      const LevelData<FArrayBox>& eneDen2 = this_species2->getEnergyDensity(setMoments);
      setInterMFT(LDe_m,numDen1,momDen1,eneDen1,numDen2,momDen2,eneDen2);
   }

}


void Coulomb::setIntraMFT( const LevelData<FArrayBox>&  a_debyeLength,
                           const LevelData<FArrayBox>&  a_numDensity,
                           const LevelData<FArrayBox>&  a_momDensity,
                           const LevelData<FArrayBox>&  a_eneDensity ) const
{
   CH_TIME("Coulomb::setIntraMFT()");
 
   //  tau = 1/nu90 with nu90 = n*vR*sigma90, sigma90 = 8/pi*b90^2*Clog,
   //  b90 = abs(q1*q2)/(2*pi*ep0*muR*vR^2), muR = m1*m2/(m1+m2), and
   //  vR^2 = 3*VT1^2 + 3*VT2^2 + |U1 - U2|^2

   Real T_eV, numDen, rho, meanE, eneDen;
   Real g12sq, g12sq_norm, bmax, bmin;
   Real b90, sigma90, sigma_max, nu90;
   Real box_nuMax=0.0; // for scattering time step calculation

   Real cvacSq  = Constants::CVAC*Constants::CVAC;
   Real mcSq_eV = Constants::ME*Constants::EV_PER_JOULE*cvacSq;
    
   const DisjointBoxLayout& grids = a_numDensity.disjointBoxLayout();
   DataIterator ditg(grids);
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numDensity = a_numDensity[ditg];
      const FArrayBox& this_momDensity = a_momDensity[ditg];
      const FArrayBox& this_eneDensity = a_eneDensity[ditg];
      const FArrayBox& this_LDe = a_debyeLength[ditg];
     
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
         const IntVect ig = gbit(); // grid index
       
         numDen = this_numDensity.get(ig,0); // [1/m^3]
         if (numDen == 0.0) { continue; }
         rho = m_mass1*numDen;

         // compute the atomic spacing
         Real atomic_spacing = 1.0/std::cbrt(4.0/3.0*Constants::PI*numDen); // atomic spacing [m]

         // set maximum impact parameter to LDe (i.e., screening length)
         // which is already limited by atomic spacing for each species
         bmax = this_LDe.get(ig,0);

         // set max sigma based on mfp = atomic spacing
         sigma_max = 1.0/(numDen*atomic_spacing);

         // set the Fermi energy for this cell
         // EF_norm = hbar^2/(2*mass)*(3*pi^2*n)^(2/3)/(me*cvac^2)
         m_EF_norm = m_EF_fact*std::pow(numDen,2.0/3.0);

         // compute Temperature and g12sq
         Real rhoUx = this_momDensity.get(ig,0); // mass1*<betax*numDen>
         Real rhoUy = this_momDensity.get(ig,1); // mass1*<betay*numDen>
         Real rhoUz = this_momDensity.get(ig,2); // mass1*<betaz*numDen>
         meanE = (rhoUx*rhoUx + rhoUy*rhoUy + rhoUz*rhoUz)/rho/2.0;
         
         eneDen = 0.0;
         for (int dir=0; dir<3; dir++) { eneDen += this_eneDensity.get(ig,dir); }
         T_eV = 2.0/3.0*(eneDen - meanE)/numDen*mcSq_eV;
         T_eV = std::max(T_eV,0.01);
         
         // compute g12^2 = 3*T1/m1 + 3*T2/m2 + |U1 - U2|^2 and b90
         g12sq = 6.0*Constants::QE/Constants::ME*T_eV/m_mass1; // [m^2/s^2]
         g12sq_norm = g12sq/cvacSq; // normalized
         b90 = m_b90_fact/(m_mu*g12sq_norm + 2.0*m_EF_norm); // [m]

         // set the Coulomb logarithm
         Real Clog = m_Clog;
         if (Clog==0.0 && g12sq>0.0) { // See Lee and More 1984 Eqs 20-22
            Real bmin_qm = m_bqm_fact/(m_mu*std::sqrt(g12sq_norm));
            bmin = std::max(b90/2.0,bmin_qm);
            Clog = 0.5*std::log(1.0 + bmax*bmax/bmin/bmin);
            Clog = std::max(2.0,Clog);
         }

         // define self-species collision time
         sigma90 = 8.0/Constants::PI*b90*b90*Clog; // [m^2]
         sigma90 = std::min(sigma90,sigma_max);    // [m^2]
         nu90 = sqrt(g12sq)*numDen*sigma90;        // [Hz]
         box_nuMax = Max(box_nuMax,nu90);
            
         bool verbose = false;
         if (!procID() && verbose) {
            cout << "JRA: sp1,sp2 = " << m_sp1 << ", " << m_sp2 << endl;
            cout << "JRA: mass1 = " << m_mass1 << endl;
            cout << "JRA: mass2 = " << m_mass2 << endl;
            cout << "JRA: N[ig="<<ig<<"]  = " << numDen <<" [1/m^3]" << endl;
            cout << "JRA: T[ig="<<ig<<"]  = " << T_eV <<" [eV]" << endl;
            cout << "JRA: bmax[ig="<<ig<<"]  = " << bmax <<" [m]" << endl;
            cout << "JRA: bmin[ig="<<ig<<"]  = " << bmin <<" [m]" << endl;
            cout << "JRA: Clog[ig="<<ig<<"]  = " << Clog << endl;
            cout << "JRA: g12[ig="<<ig<<"]   = " << std::sqrt(g12sq) << " [m/s]" << endl;
            cout << "JRA: b90[ig="<<ig<<"]   = " << b90  <<" [m]" << endl;
            cout << "JRA: tau[ig="<<ig<<"]   = " << 1.0/nu90 <<" [1/s]" << endl;
         }

      }
   
   }
   
   Real global_nuMax = box_nuMax;
//#ifdef CH_MPI
//   MPI_Allreduce( &box_nuMax, &global_nuMax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
//#endif
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]
   
}

void Coulomb::setInterMFT( const LevelData<FArrayBox>&  a_debyeLength,
                           const LevelData<FArrayBox>&  a_numDensity1,
                           const LevelData<FArrayBox>&  a_momDensity1,
                           const LevelData<FArrayBox>&  a_eneDensity1,
                           const LevelData<FArrayBox>&  a_numDensity2,
                           const LevelData<FArrayBox>&  a_momDensity2,
                           const LevelData<FArrayBox>&  a_eneDensity2 ) const
{
   CH_TIME("Coulomb::setInterMFT()");
   
   //  tau = 1/nu90 with nu90 = n*vR*sigma90, sigma90 = 8/pi*b90^2*Clog,
   //  b90 = abs(q1*q2)/(2*pi*ep0*muR*vR^2), muR = m1*m2/(m1+m2), and
   //  vR^2 = 3*VT1^2 + 3*VT2^2 + |U1 - U2|^2
 
   Real T1_eV, numDen1, rho1, eneDen1, meanE1, VT1;
   Real T2_eV, numDen2, rho2, eneDen2, meanE2, VT2;
   Real g12sq, g12sq_norm, bmax, bmin;
   Real b90, sigma90, sigma_max, nu90;
         
   Real box_nuMax=0.0; // for scattering time step calculation
   
   Real cvacSq  = Constants::CVAC*Constants::CVAC;
   Real mcSq_J  = Constants::ME*cvacSq;
   Real mcSq_eV = mcSq_J*Constants::EV_PER_JOULE;

   const DisjointBoxLayout& grids = a_numDensity1.disjointBoxLayout();
   DataIterator ditg(grids);
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numDensity1 = a_numDensity1[ditg];
      const FArrayBox& this_momDensity1 = a_momDensity1[ditg];
      const FArrayBox& this_eneDensity1 = a_eneDensity1[ditg];
      const FArrayBox& this_numDensity2 = a_numDensity2[ditg];
      const FArrayBox& this_momDensity2 = a_momDensity2[ditg];
      const FArrayBox& this_eneDensity2 = a_eneDensity2[ditg];
      const FArrayBox& this_LDe = a_debyeLength[ditg];
     
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         const IntVect ig = gbit(); // grid index
       
         // get local density and temperature and compute local VT
         numDen1 = this_numDensity1.get(ig,0); // [1/m^3]
         numDen2 = this_numDensity2.get(ig,0); // [1/m^3]
         if (numDen1*numDen2 == 0.0) { continue; }
         rho1 = m_mass1*numDen1;
         rho2 = m_mass2*numDen2;

         // compute the atomic spacing associated with minimum density species
         Real minn = std::min(numDen1,numDen2);
         Real atomic_spacing = 1.0/std::cbrt(4.0/3.0*Constants::PI*minn); // atomic spacing [m]

         // set maximum impact parameter to LDe (i.e., screening length)
         // which is already limited by atomic spacing for each species
         bmax = this_LDe.get(ig,0);

         // set max sigma based on mfp = atomic spacing
         sigma_max = 1.0/(minn*atomic_spacing);

         // set the Fermi energy for this cell
         // EF_norm = hbar^2/(2*mass)*(3*pi^2*n)^(2/3)/(me*cvac^2)
         Real maxn = std::max(numDen1,numDen2);
         m_EF_norm = m_EF_fact*std::pow(maxn,2.0/3.0);

         // compute mean energy, temperature, and g12sq
         Real rhoUx1 = this_momDensity1.get(ig,0); // mass1*<betax1*numDen1>
         Real rhoUy1 = this_momDensity1.get(ig,1); // mass1*<betay1*numDen1>
         Real rhoUz1 = this_momDensity1.get(ig,2); // mass1*<betaz1*numDen1>
         meanE1 = (rhoUx1*rhoUx1 + rhoUy1*rhoUy1 + rhoUz1*rhoUz1)/rho1/2.0;
         
         Real rhoUx2 = this_momDensity2.get(ig,0); // mass2*<betax2*numDen2>
         Real rhoUy2 = this_momDensity2.get(ig,1); // mass2*<betay2*numDen2>
         Real rhoUz2 = this_momDensity2.get(ig,2); // mass2*<betaz2*numDen2>
         meanE2 = (rhoUx2*rhoUx2 + rhoUy2*rhoUy2 + rhoUz2*rhoUz2)/rho2/2.0;

         eneDen1 = 0.0;
         eneDen2 = 0.0;
         for (int dir=0; dir<3; dir++) { eneDen1 += this_eneDensity1.get(ig,dir); }
         for (int dir=0; dir<3; dir++) { eneDen2 += this_eneDensity2.get(ig,dir); }
         T1_eV = 2.0/3.0*(eneDen1-meanE1)/numDen1*mcSq_eV;
         T2_eV = 2.0/3.0*(eneDen2-meanE2)/numDen2*mcSq_eV;
         T1_eV = std::max(T1_eV, 0.01);
         T2_eV = std::max(T2_eV, 0.01);
         VT1 = std::sqrt(Constants::QE*T1_eV/(Constants::ME*m_mass1)); // [m/s]
         VT2 = std::sqrt(Constants::QE*T2_eV/(Constants::ME*m_mass2)); // [m/s]
         
         // Nanbu 1998: gab^2 = 3*Ta/ma + 3*Tb/mb + |Ua - Ub|^2 
         g12sq = (3.0*VT1*VT1 + 3.0*VT2*VT2);
         g12sq += std::pow((rhoUx1/rho1 - rhoUx2/rho2),2)*cvacSq;
         g12sq += std::pow((rhoUy1/rho1 - rhoUy2/rho2),2)*cvacSq;
         g12sq += std::pow((rhoUz1/rho1 - rhoUz2/rho2),2)*cvacSq;
         g12sq_norm = g12sq/cvacSq;
         b90 = m_b90_fact/(m_mu*g12sq_norm + 2.0*m_EF_norm); // [m]
         
         // set the Coulomb logarithm
         Real Clog = m_Clog;
         if (Clog==0.0 && g12sq>0.0) { // See Lee and More 1984 Eqs 20-22
            Real bmin_qm = m_bqm_fact/(m_mu*std::sqrt(g12sq_norm));
            bmin = std::max(b90/2.0,bmin_qm);
            Clog = 0.5*std::log(1.0 + bmax*bmax/bmin/bmin);
            Clog = std::max(2.0,Clog);
         }
         
         // define intra-species collision time
         sigma90 = 8.0/Constants::PI*b90*b90*Clog; // [m^2]
         sigma90 = std::min(sigma90,sigma_max);    // [m^2]
         nu90 = sqrt(g12sq)*maxn*sigma90;        // [Hz]
         box_nuMax = Max(box_nuMax,nu90);

         bool verbose = false;
         if (!procID() && verbose) {
            cout << "JRA: sp1,sp2 = " << m_sp1 << ", " << m_sp2 << endl;
            cout << "JRA: mass1 = " << m_mass1 << endl;
            cout << "JRA: mass2 = " << m_mass2 << endl;
            cout << "JRA: N1[ig="<<ig<<"]  = " << numDen1 <<" [1/m^3]" << endl;
            cout << "JRA: N2[ig="<<ig<<"]  = " << numDen2 <<" [1/m^3]" << endl;
            cout << "JRA: T1[ig="<<ig<<"]  = " << T1_eV <<" [eV]" << endl;
            cout << "JRA: T2[ig="<<ig<<"]  = " << T2_eV <<" [eV]" << endl;
            cout << "JRA: bmax[ig="<<ig<<"]  = " << bmax <<" [m]" << endl;
            cout << "JRA: bmin[ig="<<ig<<"]  = " << bmin <<" [m]" << endl;
            cout << "JRA: Clog[ig="<<ig<<"]  = " << Clog << endl;
            cout << "JRA: g12[ig="<<ig<<"]   = " << std::sqrt(g12sq) << " [m/s]" << endl;
            cout << "JRA: b90[ig="<<ig<<"]   = " << b90  <<" [m]" << endl;
            cout << "JRA: tau[ig="<<ig<<"]   = " << 1.0/nu90 <<" [1/s]" << endl;
         }

      }
   
   }
   
   Real global_nuMax = box_nuMax;
//#ifdef CH_MPI
//   MPI_Allreduce( &box_nuMax, &global_nuMax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
//#endif
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]

}

void Coulomb::applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                         const DomainGrid&           a_mesh,
                         const Real                  a_dt_sec ) const
{
   CH_TIME("Coulomb::applyScattering()");
      
   PicSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getPtrVect();
   PicSpeciesPtr this_species1(pic_species_ptr_vect[m_sp1]);
   if (!this_species1->scatter()) return;

   const LevelData<FArrayBox>& LDe_m = a_pic_species_intf.getDebyeLength();
   const Real dt_scatter = a_dt_sec/(double)m_Nsubcycles;

   for (int n=0; n<m_Nsubcycles; n++) {
      if (m_sp1==m_sp2) {
         switch (m_weight_method) {
            case PROBABILISTIC: // scatter larger weight particles with probability wmin/wmax
            applyIntraScattering_PROB( *this_species1, a_mesh, LDe_m, dt_scatter );
            break;
            case CONSERVATIVE: // sentoku and kemp JCP 2008
            applyIntraScattering_SK08( *this_species1, a_mesh, LDe_m, dt_scatter );
            break;
         }
      }
      else {
      PicSpeciesPtr this_species2(pic_species_ptr_vect[m_sp2]);
      if (!this_species2->scatter()) return;
         switch (m_weight_method) {
            case PROBABILISTIC: // scatter larger weight particles with probability wmin/wmax
            applyInterScattering_PROB( *this_species1, *this_species2, a_mesh, LDe_m, dt_scatter );
            break;
            case CONSERVATIVE: // sentoku and kemp JCP 2008
            applyInterScattering_SK08( *this_species1, *this_species2, a_mesh, LDe_m, dt_scatter );
            break;
         }
      }
   }

}

void Coulomb::applyIntraScattering_PROB( PicSpecies&            a_picSpecies,
                                   const DomainGrid&            a_mesh,
                                   const LevelData<FArrayBox>&  a_LDe_m,
                                   const Real                   a_dt_sec ) const
{
   CH_TIME("Coulomb::applyIntraScattering_PROB()");
 
   // predefine some pointers to be used below
   JustinsParticle* part1_ptr = NULL;  
   JustinsParticle* part2_ptr = NULL;  
 
   // define references to a_picSpecies
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab_ptr = a_picSpecies.partData_binfab();
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity = a_picSpecies.getNumberDensity(setMoments);
   
   const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
   const Real volume_scale = a_mesh.getVolumeScale();
   const Real cellV = a_mesh.getMappedCellVolume()*volume_scale;
   // cellV*Jacobian = cell volume in SI units

   // predefine some variables
   int numCell, Naa;
   Real numDen, den12, bmax, sigma_max, cellV_SI, wpMax;

   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data_binfab_ptr.disjointBoxLayout();
   DataIterator ditg(grids);
   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity = numberDensity[ditg];
      const FArrayBox& this_LDe = a_LDe_m[ditg];
      FArrayBox& this_count = m_count[ditg];
     
      BinFab<JustinsParticlePtr>& thisBinFab_ptr = data_binfab_ptr[ditg];
   
      std::vector<JustinsParticlePtr> vector_part_ptrs;
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells

         const IntVect ig = gbit(); // grid index
         cellV_SI = cellV*Jacobian[ditg].get(ig,0);

         // get local density and set bmax
         numDen = this_numberDensity.get(ig,0);
         if (numDen == 0.0) { continue; }

         // compute the atomic spacing
         Real atomic_spacing = 1.0/std::cbrt(4.0/3.0*Constants::PI*numDen); // atomic spacing [m]

         // set maximum impact parameter to LDe (i.e., screening length)
         // which is already limited by atomic spacing for each species
         bmax = this_LDe.get(ig,0);

         // set max sigma based on mfp = atomic spacing
         sigma_max = 1.0/(numDen*atomic_spacing);

         // set the Fermi energy for this cell
         // EF_norm = hbar^2/(2*mass)*(3*pi^2*n)^(2/3)/(me*cvac^2)
         m_EF_norm = m_EF_fact*std::pow(numDen,2.0/3.0);

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         numCell = cell_pList.length();
         if (numCell < 2) { continue; }

         // determine if using NxN pairings 
         bool NxN = m_NxN;
         if (numCell < m_NxN_Nthresh) { NxN = true; }

         // define flag needed to properly treat the first few collisions 
         // for standard order N pairings when numCell is odd
         bool odd_NxN = false;
         if (!NxN && numCell % 2 == 1) { odd_NxN = true; } 
   
         Naa = numCell - 1;
    
         // copy the iterators to a vector in order to shuffle
         vector_part_ptrs.clear();
         vector_part_ptrs.reserve(numCell);
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) { vector_part_ptrs.push_back(lit()); }
         std::shuffle(vector_part_ptrs.begin(),vector_part_ptrs.end(),global_rand_gen); 
         
         Real Wtot0 = 0.0;
         Real Etot0 = 0.0;
         std::array<Real,3> ptot0 = {0.0,0.0,0.0};
         Real Efact = 0.5;
         if (m_enforce_conservations) { // save initial mass, momentum, and Energy
            for (lit.begin(); lit.ok(); ++lit) {
               part1_ptr = lit().getPointer();
               const std::array<Real,3>& betap = part1_ptr->velocity();
               const Real& wp = part1_ptr->weight();
               Wtot0 += std::pow(wp,m_beta_weight_exponent);
               Real gbsq = 0.0;
               for (int n=0; n<3; n++) {
                  ptot0[n] += wp*betap[n];
                  gbsq += betap[n]*betap[n];
               }
#ifdef RELATIVISTIC_PARTICLES
               Real gamma = std::sqrt(1.0 + gbsq);
               Efact = 1.0/(gamma + 1.0);
#endif
               Etot0 += wp*Efact*gbsq;
               if (numCell <= m_conservation_Nmin_save) { // save beta's in case correction fails
                  std::array<Real,3>& betap_old = part1_ptr->velocity_old();
                  for (int n=0; n<3; n++) { betap_old[n] = betap[n]; }
               }
            }
         }

         std::array<Real,3> dBetaAvg = {0.0,0.0,0.0};
         int p1_max = vector_part_ptrs.size()-2;
         for (int p1=0; p1<=p1_max; p1++) { // loop over p1
            
            // get particle data for first particle    
            JustinsParticlePtr& part1 = vector_part_ptrs[p1];
            part1_ptr = part1.getPointer();
            std::array<Real,3>& betap1 = part1_ptr->velocity();
            const Real wp1 = part1_ptr->weight();
            
            // set upper boundar for p2 loop
            int p2_max = p1 + 1;
            if (NxN) { p2_max = vector_part_ptrs.size()-1; }
            else if (odd_NxN) { p2_max = 2; }

            for (int p2=p1+1; p2<=p2_max; p2++) { // loop over p2

               if (!procID() && verbosity && m_sp1==0) { 
                   cout << "JRA: NxN = " << NxN << endl;
                   cout << "JRA: odd_NxN = " << odd_NxN << endl;
                   cout << "JRA: numCell = " << numCell << endl;
                   cout << "JRA: p1, p2  = " << p1 << " " << p2 << endl << endl;
               }

               // get particle data for second particle    
               JustinsParticlePtr& part2 = vector_part_ptrs[p2];
               part2_ptr = part2.getPointer();
               std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real wp2 = part2_ptr->weight();
            
               wpMax = std::max(wp1,wp2);
               if (NxN) { den12 = wpMax/cellV_SI; }
               else if (odd_NxN) { den12 = wpMax*Naa/cellV_SI/2.0; }
               else { den12 = wpMax*Naa/cellV_SI; }
            
#ifdef RELATIVISTIC_PARTICLES   
               if ( (float)wp2 < (float)wp1 ) {
                  bool scatter1 = true;
                  if (MathUtils::rand()>wp2/wp1) { scatter1 = false; }
                  LorentzScatter( betap2, betap1, scatter1, m_mass2, m_mass1, den12, bmax, sigma_max, a_dt_sec );
               }
               else {
                  bool scatter2 = true;
                  if ( (float)wp1 < (float)wp2 ) { if (MathUtils::rand()>wp1/wp2) scatter2 = false;}
                  LorentzScatter( betap1, betap2, scatter2, m_mass1, m_mass2, den12, bmax, sigma_max, a_dt_sec );
               }
#else
               // compute deltaU
               std::array<Real,3> deltaU;
               GalileanScatter( deltaU, betap1, betap2, den12, bmax, sigma_max, a_dt_sec );

               // update particle velocities
               if ( (float)wp1 == (float)wp2 ) {
                  for (int n=0; n<3; n++) { betap1[n] += 0.5*deltaU[n]; }
                  for (int n=0; n<3; n++) { betap2[n] -= 0.5*deltaU[n]; }
               }
               else if ( (float)wp1 < (float)wp2 ) {
                  for (int n=0; n<3; n++) { betap1[n] += 0.5*deltaU[n]; }
                  for (int n=0; n<3; n++) { dBetaAvg[n] += wp1*0.5*deltaU[n]; }
                  if (MathUtils::rand()<wp1/wp2) {
                     for (int n=0; n<3; n++) { betap2[n] -= 0.5*deltaU[n]; }
                     for (int n=0; n<3; n++) { dBetaAvg[n] -= wp2*0.5*deltaU[n]; }
                  }
               }
               else { //if ( (float)wp2 < (float)wp1 ) {
                  if (MathUtils::rand()<wp2/wp1) {
                     for (int n=0; n<3; n++) { betap1[n] += 0.5*deltaU[n]; }
                     for (int n=0; n<3; n++) { dBetaAvg[n] += wp1*0.5*deltaU[n]; }
                  }
                  for (int n=0; n<3; n++) { betap2[n] -= 0.5*deltaU[n]; }
                  for (int n=0; n<3; n++) { dBetaAvg[n] -= wp2*0.5*deltaU[n]; }
               }
#endif

            } // end loop over p2
            
            if (odd_NxN && p1==1) { odd_NxN = false; }
            if (!odd_NxN && !NxN) { ++p1; }

         } // end loop over p1
         verbosity=0;

#ifdef RELATIVISTIC_PARTICLES
         if (m_enforce_conservations) {
         
            std::array<Real,3> ptot1 = {0.0,0.0,0.0};
            for (lit.begin(); lit.ok(); ++lit) {
               part1_ptr = lit().getPointer();
               const std::array<Real,3>& betap = part1_ptr->velocity();
               const Real& wp = part1_ptr->weight();
               for (int n=0; n<3; n++) { ptot1[n] += wp*betap[n]; }
            }

            // compute change in momentum
            for (int n=0; n<3; n++) { dBetaAvg[n] = ptot1[n] - ptot0[n]; }

         }
#endif

         // adjust particle betas such that momentum and energy are identically conserved
         Real dBetaSq = 0.0;
         for (int n=0; n<3; n++) { dBetaSq += dBetaAvg[n]*dBetaAvg[n]; }
         if (dBetaSq>0.0 && m_enforce_conservations) {

            // enforce momentum conservation and compute deltaE
            for (int n=0; n<3; n++) { dBetaAvg[n] /= Wtot0; }
            Real Etot1 = 0.0;
            for (lit.begin(); lit.ok(); ++lit) {
               part1_ptr = lit().getPointer();
               std::array<Real,3>& betap = part1_ptr->velocity();
               const Real& wp = part1_ptr->weight();
               for (int n=0; n<3; n++) { betap[n] -= std::pow(wp,m_beta_weight_exponent-1)*dBetaAvg[n]; }
               Real gbsq = 0.0;
               for (int n=0; n<3; n++) { gbsq += betap[n]*betap[n]; }
#ifdef RELATIVISTIC_PARTICLES
               Real gamma = std::sqrt(1.0 + gbsq);
               Efact = 1.0/(gamma + 1.0);
#endif
               Etot1 += wp*Efact*gbsq;
            }
            long double deltaE = m_mass1*(Etot1 - Etot0);

            // sort particles by weight to bias adjusting particles below to
            // heavier weight ones and to get pairs that have similar weights
            if (m_sort_weighted_particles) {
               std::sort( vector_part_ptrs.begin(), vector_part_ptrs.end(),
                          []( JustinsParticlePtr& lhs, JustinsParticlePtr& rhs )
                               { return lhs.getPointer()->weight() > rhs.getPointer()->weight(); } );
            }

            int count = 0;
            int loop_count = 0;
            Real Erel_cumm = 0.0;
            Real fmult_fact = 1.0;
            for (int p=0; p<vector_part_ptrs.size(); p++) {
               int p1 = p;
               p++;
               if (p==vector_part_ptrs.size()) {
                  loop_count++;
                  p = 0;
               }
               int p2 = p;
            
               // get particle data for first particle    
               JustinsParticlePtr& part1 = vector_part_ptrs[p1];
               part1_ptr = part1.getPointer();
               std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real wp1 = part1_ptr->weight();
            
               // get particle data for second particle    
               JustinsParticlePtr& part2 = vector_part_ptrs[p2];
               part2_ptr = part2.getPointer();
               std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real wp2 = part2_ptr->weight();

               // zero angle inelastic scatter
               modEnergyPairwise( betap1, betap2, m_mass1*wp1, m_mass2*wp2, 
                                  fmult_fact, Erel_cumm, deltaE );
               count += 1;
               if (deltaE==0.0) { break; }
               if (p==vector_part_ptrs.size()-1) {
                  loop_count++;

                  // compute energy frac to correct energy with one more loop
                  const Real energy_frac_eff = std::abs(deltaE)/Erel_cumm;
                  if (energy_frac_eff>m_energy_fraction_max || loop_count>10) {
                     if (numCell>50) {
                        cout << "Notice: energy_frac_eff = " << energy_frac_eff << endl;
                        cout << "        on loop_count = " << loop_count << endl;
                        cout << "        is larger than energy_frac_max = " << m_energy_fraction_max << endl;
                        cout << "        for species1 = " << m_species1_name << endl;
                        cout << "        and species2 = " << m_species2_name << endl;
                        cout << "        num_parts    = " << vector_part_ptrs.size() << endl;
                        cout << "        ig           = " << ig << endl;
                        cout << "        Etot0        = " << Etot0 << endl;
                        cout << "        deltaE       = " << deltaE << endl;
                     }
                     if (numCell <= m_conservation_Nmin_save) {
                        for (lit.begin(); lit.ok(); ++lit) {
                           part1_ptr = lit().getPointer();
                           const std::array<Real,3>& betap_old = part1_ptr->velocity_old();
                           std::array<Real,3>& betap = part1_ptr->velocity();
                           for (int n=0; n<3; n++) { betap[n] = betap_old[n]; }
                        }
                     }
                     break;
                  }
                  else if (energy_frac_eff>m_energy_fraction) {
                     fmult_fact = energy_frac_eff/m_energy_fraction;
                  }

                  //JustinsParticlePtr temp = vector_part_ptrs.front();
                  //std::rotate( vector_part_ptrs.begin(), vector_part_ptrs.begin() + 1, 
                  //             vector_part_ptrs.end() );
                  //vector_part_ptrs.back() = temp;
                  Erel_cumm = 0.0;
                  p = -1;
               }
            }
            this_count.set(ig,0,count);

         }

      } // end loop over cells

   } // end loop over boxes

   if (m_print_correction_count) { printCorrectionCount(); }

   // don't forget to set pointers back to NULL and delete
   part1_ptr = NULL;
   part2_ptr = NULL;
   delete part1_ptr;
   delete part2_ptr;

}

void Coulomb::applyIntraScattering_SK08( PicSpecies&            a_picSpecies, 
                                   const DomainGrid&            a_mesh,
                                   const LevelData<FArrayBox>&  a_LDe_m,
                                   const Real                   a_dt_sec ) const
{
   CH_TIME("Coulomb::applyIntraScattering_SK08()");
 
   // predefine some pointers to be used below
   JustinsParticle* part1_ptr = NULL;  
   JustinsParticle* part2_ptr = NULL;  
 
   // define references to a_picSpecies
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab_ptr = a_picSpecies.partData_binfab();
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity = a_picSpecies.getNumberDensity(setMoments);
   
   const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
   const Real volume_scale = a_mesh.getVolumeScale();
   const Real cellV = a_mesh.getMappedCellVolume()*volume_scale;
   // cellV*Jacobian = cell volume in SI units

   // predefine some variables
   int numCell, Naa;
   Real numDen, den12, bmax, sigma_max, cellV_SI, wpMax, Escatter_sum, Ebefore;
   std::array<Real,3> betap_sum, betap_before;
 
   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data_binfab_ptr.disjointBoxLayout();
   DataIterator ditg(grids);
   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity = numberDensity[ditg];
      const FArrayBox& this_LDe = a_LDe_m[ditg];
     
      BinFab<JustinsParticlePtr>& thisBinFab_ptr = data_binfab_ptr[ditg];
   
      std::vector<JustinsParticlePtr> vector_part_ptrs;
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
         const IntVect ig = gbit(); // grid index
         cellV_SI = cellV*Jacobian[ditg].get(ig,0);
       
         // get local density and temperature and compute local tau
         numDen = this_numberDensity.get(ig,0);
         if (numDen == 0.0) { continue; }

         // compute the atomic spacing
         Real atomic_spacing = 1.0/std::cbrt(4.0/3.0*Constants::PI*numDen); // atomic spacing [m]

         // set maximum impact parameter to LDe (i.e., screening length)
         // which is already limited by atomic spacing for each species
         bmax = this_LDe.get(ig,0);

         // set max sigma based on mfp = atomic spacing
         sigma_max = 1.0/(numDen*atomic_spacing);

         // set the Fermi energy for this cell
         // EF_norm = hbar^2/(2*mass)*(3*pi^2*n)^(2/3)/(me*cvac^2)
         m_EF_norm = m_EF_fact*std::pow(numDen,2.0/3.0);

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         numCell = cell_pList.length();
         if (numCell < 2) { continue; }

         // define flag needed to properly treat the first few collisions when numCell is odd
         bool odd_NxN = true;        
         if (numCell % 2 == 0) { odd_NxN = false; }

         Naa = numCell - 1;

         // copy the iterators to a vector in order to shuffle
         vector_part_ptrs.clear();
         vector_part_ptrs.reserve(numCell);
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) { vector_part_ptrs.push_back(lit()); }
         std::shuffle(vector_part_ptrs.begin(),vector_part_ptrs.end(),global_rand_gen); 

         for (int p=0; p<vector_part_ptrs.size(); p++) { // loop over particle scattering pairs

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

            // get particle data for first particle    
            JustinsParticlePtr& part1 = vector_part_ptrs[p1];
            part1_ptr = part1.getPointer();
            std::array<Real,3>& betap1 = part1_ptr->velocity();
            const Real wp1 = part1_ptr->weight();

            // get particle data for second particle    
            JustinsParticlePtr& part2 = vector_part_ptrs[p2];
            part2_ptr = part2.getPointer();
            std::array<Real,3>& betap2 = part2_ptr->velocity();
            const Real wp2 = part2_ptr->weight();

            wpMax = std::max(wp1,wp2);
            den12 = wpMax*Naa/cellV_SI;
            if (odd_NxN) {
               den12 = den12/2.0;
               if (p==2) { odd_NxN = false; }
            }

            // compute deltaU
            std::array<Real,3> deltaU;
            GalileanScatter( deltaU, betap1, betap2, den12, bmax, sigma_max, a_dt_sec );

            // update particle velocities
            if ( (float)wp1 < (float)wp2 ) {

               const long double ratio = wp1/wp2;
               betap_before = betap2;
               Ebefore = m_mass2*(betap2[0]*betap2[0] + betap2[1]*betap2[1] + betap2[2]*betap2[2])/2.0;
               betap_sum = {0.0,0.0,0.0};
               Escatter_sum = 0.0;

               // scatter particles
               std::array<Real,3> betap2p = betap2;
               for (int n=0; n<3; n++) { betap1[n]  += 0.5*deltaU[n]; }
               for (int n=0; n<3; n++) { betap2p[n] -= 0.5*deltaU[n]; }

               // update betap_sum and Escatter_sum
               for (int n=0; n<3; n++) { betap_sum[n] += betap2p[n]; }
               Escatter_sum += m_mass2*(betap2p[0]*betap2p[0] + betap2p[1]*betap2p[1] + betap2p[2]*betap2p[2])/2.0;
                     
               const Real Eafter = Ebefore + ratio*(Escatter_sum-Ebefore);
               for (int n=0; n<3; n++) { betap2[n] = betap_before[n] + ratio*(betap_sum[n]-betap_before[n]); }

               // correct beta of particle 2 to conserve energy exactly and momentum on average
               enforceEnergyConservation( betap2, m_mass2, Eafter );

            }
            else if ( (float)wp2 < (float)wp1 ) {

               const long double ratio = wp2/wp1;
               betap_before = betap1;
               Ebefore = m_mass1*(betap1[0]*betap1[0] + betap1[1]*betap1[1] + betap1[2]*betap1[2])/2.0;
               betap_sum = {0.0,0.0,0.0};
               Escatter_sum = 0.0;

               // scatter particles
               std::array<Real,3> betap1p = betap1;
               for (int n=0; n<3; n++) { betap1p[n] += 0.5*deltaU[n]; }
               for (int n=0; n<3; n++) { betap2[n]  -= 0.5*deltaU[n]; }
                     
               // update betap_sum and Escatter_sum
               for (int n=0; n<3; n++) { betap_sum[n] += betap1p[n]; }
               Escatter_sum += m_mass1*(betap1p[0]*betap1p[0] + betap1p[1]*betap1p[1] + betap1p[2]*betap1p[2])/2.0;
                     
               const Real Eafter = Ebefore + ratio*(Escatter_sum-Ebefore);
               for (int n=0; n<3; n++) { betap1[n] = betap_before[n] + ratio*(betap_sum[n]-betap_before[n]); }

               // correct beta of particle 1 to conserve energy exactly and momentum on average
               enforceEnergyConservation( betap1, m_mass1, Eafter );

            }
            else {
               for (int n=0; n<3; n++) { betap1[n] += 0.5*deltaU[n]; }
               for (int n=0; n<3; n++) { betap2[n] -= 0.5*deltaU[n]; }
            }

         }
         verbosity=0;

      } // end loop over cells

   } // end loop over boxes
  
 
   // don't forget to set pointers back to NULL and delete
   part1_ptr = NULL;
   part2_ptr = NULL;
   delete part1_ptr;
   delete part2_ptr;
   
}

void Coulomb::applyInterScattering_PROB( PicSpecies&            a_picSpecies1,
                                         PicSpecies&            a_picSpecies2, 
                                   const DomainGrid&            a_mesh,
                                   const LevelData<FArrayBox>&  a_LDe_m,
                                   const Real                   a_dt_sec ) const
{
   CH_TIME("Coulomb::applyInterScattering_PROB()");
   
   CH_assert(m_species1_name==a_picSpecies1.name());
   CH_assert(m_species2_name==a_picSpecies2.name());

   // predefine some pointers to be used below
   JustinsParticle* part1_ptr = NULL;  
   JustinsParticle* part2_ptr = NULL;  
   
   // define references to a_picSpecies1
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab1_ptr = a_picSpecies1.partData_binfab();
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity1 = a_picSpecies1.getNumberDensity(setMoments);
   
   // define references to a_picSpecies2
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab2_ptr = a_picSpecies2.partData_binfab();
   const LevelData<FArrayBox>& numberDensity2 = a_picSpecies2.getNumberDensity(setMoments);
  
   const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
   const Real volume_scale = a_mesh.getVolumeScale();
   const Real cellV = a_mesh.getMappedCellVolume()*volume_scale;
   // cellV*Jacobian = cell volume in SI units
 
   // predefine some variables
   int numCell1, numCell2, Nmax, Nmin;
   Real numDen1, wpMax;
   Real numDen2, den12, bmax, sigma_max, cellV_SI;

   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data_binfab1_ptr.disjointBoxLayout();
   DataIterator ditg(grids);

   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity1 = numberDensity1[ditg];
      const FArrayBox& this_numberDensity2 = numberDensity2[ditg];
      const FArrayBox& this_LDe = a_LDe_m[ditg];
     
      BinFab<JustinsParticlePtr>& thisBinFab1_ptr = data_binfab1_ptr[ditg];
      BinFab<JustinsParticlePtr>& thisBinFab2_ptr = data_binfab2_ptr[ditg];
   
      std::vector<JustinsParticlePtr> vector_part1_ptrs;
      std::vector<JustinsParticlePtr> vector_part2_ptrs;
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
         const IntVect ig = gbit(); // grid index
         cellV_SI = cellV*Jacobian[ditg].get(ig,0);
         
         // get local density and temperature and compute local tau
         numDen1 = this_numberDensity1.get(ig,0);
         numDen2 = this_numberDensity2.get(ig,0);
         if (numDen1*numDen2 == 0.0) { continue; }
         
         // compute the atomic spacing associated with minimum density species
         Real minn = std::min(numDen1,numDen2);
         Real atomic_spacing = 1.0/std::cbrt(4.0/3.0*Constants::PI*minn); // atomic spacing [m]

         // set maximum impact parameter to LDe (i.e., screening length)
         // which is already limited by atomic spacing for each species
         bmax = this_LDe.get(ig,0);

         // set max sigma based on mfp = atomic spacing
         sigma_max = 1.0/(minn*atomic_spacing);

         // set the Fermi energy for this cell
         // EF_norm = hbar^2/(2*mass)*(3*pi^2*n)^(2/3)/(me*cvac^2)
         Real maxn = std::max(numDen1,numDen2);
         m_EF_norm = m_EF_fact*std::pow(maxn,2.0/3.0);

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
         List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
         numCell1 = cell_pList1.length();
         numCell2 = cell_pList2.length();
         if (numCell1*numCell2 < 2) { continue; }
         Nmin = std::min(numCell1,numCell2);
         Nmax = std::max(numCell1,numCell2);

         bool NxN = m_NxN;
         if (Nmin < m_NxN_Nthresh) { NxN = true; }

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
         
         Real Wtot0 = 0.0, wp1_mean = 0.0, wp2_mean = 0.0;
         std::array<Real,3> ptot0 = {0.0,0.0,0.0};
         Real Etot0 = 0.0;
         Real Efact = 0.5;
         if (m_enforce_conservations) { // save initial mass, momentum, and Energy
            Real Wtot01 = 0.0, Etot01 = 0.0;
            std::array<Real,3> ptot01 = {0.0,0.0,0.0};
            wp1_mean = 0.0;
            for (lit1.begin(); lit1.ok(); ++lit1) {
               part1_ptr = lit1().getPointer();
               const std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real& wp1 = part1_ptr->weight();
               wp1_mean += wp1;
               Wtot01 += std::pow(wp1,m_beta_weight_exponent);
               Real gbsq1 = 0.0;
               for (int n=0; n<3; n++) {
                  ptot01[n] += wp1*betap1[n];
                  gbsq1 += betap1[n]*betap1[n];
               }
#ifdef RELATIVISTIC_PARTICLES
               Real gamma1 = std::sqrt(1.0 + gbsq1);
               Efact = 1.0/(gamma1 + 1.0);
#endif
               Etot01 += wp1*Efact*gbsq1;
               if (numCell1 <= m_conservation_Nmin_save || // save beta's in case correction fails
                   numCell2 <= m_conservation_Nmin_save) {
                  std::array<Real,3>& betap1_old = part1_ptr->velocity_old();
                  for (int n=0; n<3; n++) { betap1_old[n] = betap1[n]; }
               }
            }
            wp1_mean /= numCell1; // mean particle weight for species 1

            Real Wtot02 = 0.0, Etot02 = 0.0;
            std::array<Real,3> ptot02 = {0.0,0.0,0.0};
            for (lit2.begin(); lit2.ok(); ++lit2) {
               part2_ptr = lit2().getPointer();
               const std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real& wp2 = part2_ptr->weight();
               wp2_mean += wp2;
               Wtot02 += std::pow(wp2,m_beta_weight_exponent);
               Real gbsq2 = 0.0;
               for (int n=0; n<3; n++) {
                  ptot02[n] += wp2*betap2[n];
                  gbsq2 += betap2[n]*betap2[n];
               }
#ifdef RELATIVISTIC_PARTICLES
               Real gamma2 = std::sqrt(1.0 + gbsq2);
               Efact = 1.0/(gamma2 + 1.0);
#endif
               Etot02 += wp2*Efact*gbsq2;
               if (numCell1 <= m_conservation_Nmin_save ||
                   numCell2 <= m_conservation_Nmin_save) {
                  std::array<Real,3>& betap2_old = part2_ptr->velocity_old();
                  for (int n=0; n<3; n++) { betap2_old[n] = betap2[n]; }
               }
            }
            wp2_mean /= numCell2; // mean particle weight for species 2

            Wtot0 = m_mass1*Wtot01 + m_mass2*Wtot02;
            for (int n=0; n<3; n++) { ptot0[n] = m_mass1*ptot01[n] + m_mass2*ptot02[n]; }
            Etot0 = m_mass1*Etot01 + m_mass2*Etot02;
         }

         unsigned int p1, p2;
         std::array<Real,3> dBetaAvg = {0.0,0.0,0.0};
         for (int p=0; p<Nmax; p++) { // loop over particles with max number
     
            int pmin_start; 
            if (Nmin==numCell1) {
               p1 = p % numCell1;
               p2 = p;
               pmin_start = p1;
            } 
            else {
               p1 = p;
               p2 = p % numCell2;
               pmin_start = p2;
            } 
            int pmin_end = pmin_start; 

            if (NxN) { 
               pmin_start = 0;
               pmin_end = Nmin - 1;
            }
 
            for (int pmin=pmin_start; pmin<=pmin_end; pmin++) { // loop over particles with min number

               if (NxN) {
                  if (Nmin==numCell1) { p1 = pmin; }
                  else { p2 = pmin; }
               }

               if (!procID() && verbosity) { 
                   cout << "JRA: NxN = " << NxN << endl;
                   cout << "JRA: numCell1 = " << numCell1 << endl;
                   cout << "JRA: numCell2 = " << numCell2 << endl;
                   cout << "JRA: p1, p2  = " << p1 << " " << p2 << endl << endl;
               }

               // get particle data for first particle    
               JustinsParticlePtr& part1 = vector_part1_ptrs[p1];
               part1_ptr = part1.getPointer();
               std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real wp1 = part1_ptr->weight();
            
               // get particle data for second particle    
               JustinsParticlePtr& part2 = vector_part2_ptrs[p2];
               part2_ptr = part2.getPointer();
               std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real wp2 = part2_ptr->weight();

               // determine maximum weight of this pair and set scatter flags (rejection method)
               wpMax = std::max(wp1,wp2);
               if (NxN) { den12 = wpMax/cellV_SI; }
               else { den12 = wpMax*Nmin/cellV_SI; }
            
#ifdef RELATIVISTIC_PARTICLES
               if ( (float)wp2 < (float)wp1 ) {
                  bool scatter1 = true;
                  if (MathUtils::rand()>wp2/wp1) { scatter1 = false; }
                  LorentzScatter( betap2, betap1, scatter1, m_mass2, m_mass1, den12, bmax, sigma_max, a_dt_sec );
               }
               else {
                  bool scatter2 = true;
                  if ( (float)wp1 < (float)wp2 ) {if (MathUtils::rand()>wp1/wp2) scatter2 = false;}
                  LorentzScatter( betap1, betap2, scatter2, m_mass1, m_mass2, den12, bmax, sigma_max, a_dt_sec );
               }
#else
               // compute deltaU
               std::array<Real,3> deltaU;
               GalileanScatter( deltaU, betap1, betap2, den12, bmax, sigma_max, a_dt_sec );

               // update particle velocities
               if ( (float)wp1 == (float)wp2 ) {
                  for (int n=0; n<3; n++) { betap1[n] += m_mu/m_mass1*deltaU[n]; }
                  for (int n=0; n<3; n++) { betap2[n] -= m_mu/m_mass2*deltaU[n]; }
               }
               else if ( (float)wp1 < (float)wp2 ) {
                  for (int n=0; n<3; n++) { betap1[n] += m_mu/m_mass1*deltaU[n]; }
                  for (int n=0; n<3; n++) { dBetaAvg[n] += wp1*m_mu*deltaU[n]; }
                  if (MathUtils::rand()<wp1/wp2) {
                     for (int n=0; n<3; n++) { betap2[n] -= m_mu/m_mass2*deltaU[n]; }
                     for (int n=0; n<3; n++) { dBetaAvg[n] -= wp2*m_mu*deltaU[n]; }
                  }
               }
               else { //if ( (float)wp2 < (float)wp1 ) {
                  if (MathUtils::rand()<wp2/wp1) {
                     for (int n=0; n<3; n++) { betap1[n] += m_mu/m_mass1*deltaU[n]; }
                     for (int n=0; n<3; n++) { dBetaAvg[n] += wp1*m_mu*deltaU[n]; }
                  }
                  for (int n=0; n<3; n++) { betap2[n] -= m_mu/m_mass2*deltaU[n]; }
                  for (int n=0; n<3; n++) { dBetaAvg[n] -= wp2*m_mu*deltaU[n]; }
               }
#endif

            } // end loop over pmin

         } // end loop over pmax
         verbosity=0;

#ifdef RELATIVISTIC_PARTICLES
         if (m_enforce_conservations) {

            std::array<Real,3> ptot1 = {0.0,0.0,0.0};
            std::array<Real,3> ptot11 = {0.0,0.0,0.0};
            std::array<Real,3> ptot21 = {0.0,0.0,0.0};
            for (lit1.begin(); lit1.ok(); ++lit1) {
               part1_ptr = lit1().getPointer();
               const std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real& wp1 = part1_ptr->weight();
               for (int n=0; n<3; n++) { ptot11[n] += wp1*betap1[n]; }
            }
            for (int n=0; n<3; n++) { ptot11[n] = m_mass1*ptot11[n]; }

            for (lit2.begin(); lit2.ok(); ++lit2) {
               part2_ptr = lit2().getPointer();
               const std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real& wp2 = part2_ptr->weight();
               for (int n=0; n<3; n++) { ptot21[n] += wp2*betap2[n]; }
            }
            for (int n=0; n<3; n++) { ptot21[n] = m_mass2*ptot21[n]; }

            // compute change in momentum
            for (int n=0; n<3; n++) {
               ptot1[n] = ptot11[n] + ptot21[n];
               dBetaAvg[n] = ptot1[n] - ptot0[n];
            }

         }
#endif

         // adjust particle betas such that momentum and energy are identically conserved
         Real dBetaSq = 0.0;
         for (int n=0; n<3; n++) { dBetaSq += dBetaAvg[n]*dBetaAvg[n]; }
         if (dBetaSq>0.0 && m_enforce_conservations) {
        
            // enforce momentum conservation and compute deltaE
            for (int n=0; n<3; n++) { dBetaAvg[n] /= Wtot0; }

            long double Etot11 = 0.0;
            for (lit1.begin(); lit1.ok(); ++lit1) {
               part1_ptr = lit1().getPointer();
               std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real& wp1 = part1_ptr->weight();
               for (int n=0; n<3; n++) { betap1[n] -= std::pow(wp1,m_beta_weight_exponent-1)*dBetaAvg[n]; }
               Real gbsq1 = 0.0;
               for (int n=0; n<3; n++) { gbsq1 += betap1[n]*betap1[n]; }
#ifdef RELATIVISTIC_PARTICLES
               Real gamma1 = std::sqrt(1.0 + gbsq1);
               Efact = 1.0/(gamma1 + 1.0);
#endif
               Etot11 += wp1*Efact*gbsq1;
            }
            Etot11 *= m_mass1;

            long double Etot12 = 0.0;
            for (lit2.begin(); lit2.ok(); ++lit2) {
               part2_ptr = lit2().getPointer();
               std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real& wp2 = part2_ptr->weight();
               for (int n=0; n<3; n++) { betap2[n] -= std::pow(wp2,m_beta_weight_exponent-1)*dBetaAvg[n]; }
               Real gbsq2 = 0.0;
               for (int n=0; n<3; n++) { gbsq2 += betap2[n]*betap2[n]; }
#ifdef RELATIVISTIC_PARTICLES
               Real gamma2 = std::sqrt(1.0 + gbsq2);
               Efact = 1.0/(gamma2 + 1.0);
#endif
               Etot12 += wp2*Efact*gbsq2;
            }
            Etot12 *= m_mass2;

            const long double Etot1 = Etot11 + Etot12;
            const long double deltaE = Etot1 - Etot0;
            const long double Etotdenom = wp1_mean*Etot11 + wp2_mean*Etot12;
            long double deltaEp1, deltaEp2;
            if (numCell1==1) {
               deltaEp1 = 0.0;
               deltaEp2 = deltaE;
            }
            else if (numCell2==1) {
               deltaEp1 = deltaE;
               deltaEp2 = 0.0;
            }
            else {
               deltaEp1 = wp1_mean*Etot11/Etotdenom*deltaE;
               deltaEp2 = wp2_mean*Etot12/Etotdenom*deltaE;
            }

            // sort particles by weight to bias adjusting particles below to
            // heavier weight ones and to get pairs that have similar weights
            if (m_sort_weighted_particles) {
               std::sort( vector_part1_ptrs.begin(), vector_part1_ptrs.end(),
                          []( JustinsParticlePtr& lhs, JustinsParticlePtr& rhs )
                               { return lhs.getPointer()->weight() > rhs.getPointer()->weight(); } );
               std::sort( vector_part2_ptrs.begin(), vector_part2_ptrs.end(),
                          []( JustinsParticlePtr& lhs, JustinsParticlePtr& rhs )
                               { return lhs.getPointer()->weight() > rhs.getPointer()->weight(); } );
            }
            
            // adjust species 1 particles to absorb deltaEp1
            int loop1_count = 0;
            Real Erel1_cumm = 0.0;
            Real fmult1_fact = 1.0;
            bool correction_failed = false;
            for (int p=0; p<vector_part1_ptrs.size(); p++) {

               if (deltaEp1==0.0) { break; }
               int p1 = p;
               p++;
               if (p==vector_part1_ptrs.size()) {
                  loop1_count++;
                  p = 0;
               }
               int p2 = p;
            
               // get particle data for first particle    
               JustinsParticlePtr& part1 = vector_part1_ptrs[p1];
               part1_ptr = part1.getPointer();
               std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real wp1 = part1_ptr->weight();
            
               // get particle data for second particle    
               JustinsParticlePtr& part2 = vector_part1_ptrs[p2];
               part2_ptr = part2.getPointer();
               std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real wp2 = part2_ptr->weight();

               // zero angle inelastic scatter
               modEnergyPairwise( betap1, betap2, m_mass1*wp1, m_mass1*wp2, 
                                  fmult1_fact, Erel1_cumm, deltaEp1);

               if (p==vector_part1_ptrs.size()-1) {
                  loop1_count++;
                  p = -1;
                  // compute energy frac to correct energy with one more loop
                  const Real energy_frac_eff = std::abs(deltaEp1)/Erel1_cumm;
                  if (energy_frac_eff>m_energy_fraction_max || loop1_count>10) {
                     if (vector_part1_ptrs.size()>50) {
                        cout << "Notice: energy_frac_eff = " << energy_frac_eff << endl;
                        cout << "        on loop_count = " << loop1_count << endl;
                        cout << "        is larger than energy_frac_max = " << m_energy_fraction_max << endl;
                        cout << "        for species1 = " << m_species1_name << endl;
                        cout << "        and species2 = " << m_species2_name << endl;
                        cout << "        num_parts    = " << vector_part1_ptrs.size() << endl;
                        cout << "        ig           = " << ig << endl;
                        cout << "        Etot0        = " << Etot0 << endl;
                        cout << "        Etot1        = " << Etot1 << endl;
                        cout << "        Etot11       = " << Etot11 << endl;
                        cout << "        deltaEp1     = " << deltaEp1 << endl;
                     }
                     correction_failed = true;
                     break;
                  }
                  else if (energy_frac_eff>m_energy_fraction) {
                     fmult1_fact = energy_frac_eff/m_energy_fraction;
                  }
               }

            }
            
            // adjust species 2 particles to absorb deltaEp1
            int loop2_count = 0;
            Real Erel2_cumm = 0.0;
            Real fmult2_fact = 1.0;
            for (int p=0; p<vector_part2_ptrs.size(); p++) {

               if (deltaEp2==0.0) { break; }
               if (correction_failed) { break; }
               int p1 = p;
               p++;
               if (p==vector_part2_ptrs.size()) {
                  loop2_count++;
                  p = 0;
               }
               int p2 = p;
            
               // get particle data for first particle    
               JustinsParticlePtr& part1 = vector_part2_ptrs[p1];
               part1_ptr = part1.getPointer();
               std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real wp1 = part1_ptr->weight();
            
               // get particle data for second particle    
               JustinsParticlePtr& part2 = vector_part2_ptrs[p2];
               part2_ptr = part2.getPointer();
               std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real wp2 = part2_ptr->weight();

               // zero angle inelastic scatter
               modEnergyPairwise( betap1, betap2, m_mass2*wp1, m_mass2*wp2,
                                  fmult2_fact, Erel2_cumm, deltaEp2 );

               if (p==vector_part2_ptrs.size()-1) {
                  loop2_count++;
                  p = -1;
                  // compute energy frac to correct energy with one more loop
                  const Real energy_frac_eff = std::abs(deltaEp2)/Erel2_cumm;
                  if (energy_frac_eff>m_energy_fraction_max || loop2_count>10) {
                     if (vector_part2_ptrs.size()>50) {
                        cout << "Notice: energy_frac_eff = " << energy_frac_eff << endl;
                        cout << "        on loop_count = " << loop2_count << endl;
                        cout << "        is larger than energy_frac_max = " << m_energy_fraction_max << endl;
                        cout << "        for species1 = " << m_species1_name << endl;
                        cout << "        and species2 = " << m_species2_name << endl;
                        cout << "        num_parts    = " << vector_part2_ptrs.size() << endl;
                        cout << "        ig           = " << ig << endl;
                        cout << "        Etot0        = " << Etot0 << endl;
                        cout << "        Etot1        = " << Etot1 << endl;
                        cout << "        Etot12       = " << Etot12 << endl;
                        cout << "        deltaEp2     = " << deltaEp2 << endl;
                     }
                     correction_failed = true;
                     break;
                  }
                  else if (energy_frac_eff>m_energy_fraction) {
                     fmult2_fact = energy_frac_eff/m_energy_fraction;
                  }
               }

            }

            // Reset betas if energy-correction failed
            if (correction_failed) {
               if (numCell1 <= m_conservation_Nmin_save ||
                   numCell2 <= m_conservation_Nmin_save) {
                  for (lit1.begin(); lit1.ok(); ++lit1) {
                     part1_ptr = lit1().getPointer();
                     const std::array<Real,3>& betap_old = part1_ptr->velocity_old();
                     std::array<Real,3>& betap = part1_ptr->velocity();
                     for (int n=0; n<3; n++) { betap[n] = betap_old[n]; }
                  }
                  for (lit2.begin(); lit2.ok(); ++lit2) {
                      part2_ptr = lit2().getPointer();
                      const std::array<Real,3>& betap_old = part2_ptr->velocity_old();
                      std::array<Real,3>& betap = part2_ptr->velocity();
                      for (int n=0; n<3; n++) { betap[n] = betap_old[n]; }
                   }
               }
            }

         }

      } // end loop over cells

   } // end loop over boxes

   // don't forget to set pointers back to NULL and delete
   part1_ptr = NULL;
   part2_ptr = NULL;
   delete part1_ptr;
   delete part2_ptr;
   
}

void Coulomb::applyInterScattering_SK08( PicSpecies&            a_picSpecies1,
                                         PicSpecies&            a_picSpecies2, 
                                   const DomainGrid&            a_mesh,
                                   const LevelData<FArrayBox>&  a_LDe_m,
                                   const Real                   a_dt_sec ) const
{
   CH_TIME("Coulomb::applyInterScattering_SK08()");
   
   CH_assert(m_species1_name==a_picSpecies1.name());
   CH_assert(m_species2_name==a_picSpecies2.name());

   // predefine some pointers to be used below
   JustinsParticle* part1_ptr = NULL;  
   JustinsParticle* part2_ptr = NULL;  
   
   // define references to a_picSpecies1
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab1_ptr = a_picSpecies1.partData_binfab();
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity1 = a_picSpecies1.getNumberDensity(setMoments);
   
   // define references to a_picSpecies2
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab2_ptr = a_picSpecies2.partData_binfab();
   const LevelData<FArrayBox>& numberDensity2 = a_picSpecies2.getNumberDensity(setMoments);
  
   const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();
   const Real volume_scale = a_mesh.getVolumeScale();
   const Real cellV = a_mesh.getMappedCellVolume()*volume_scale;
   // cellV*Jacobian = cell volume in SI units
 
   // predefine some variables
   int numCell1, numCell2, Nmax, Nmin;
   Real numDen1, wpMax;
   Real numDen2, den12, bmax, sigma_max, cellV_SI;
   Real Escatter_sum, Ebefore;
   std::array<Real,3> betap_sum, betap_before;
 
   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data_binfab1_ptr.disjointBoxLayout();
   DataIterator ditg(grids);
   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity1 = numberDensity1[ditg];
      const FArrayBox& this_numberDensity2 = numberDensity2[ditg];
      const FArrayBox& this_LDe = a_LDe_m[ditg];
     
      BinFab<JustinsParticlePtr>& thisBinFab1_ptr = data_binfab1_ptr[ditg];
      BinFab<JustinsParticlePtr>& thisBinFab2_ptr = data_binfab2_ptr[ditg];
   
      std::vector<JustinsParticlePtr> vector_part1_ptrs;
      std::vector<JustinsParticlePtr> vector_part2_ptrs;
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
         const IntVect ig = gbit(); // grid index
         cellV_SI = cellV*Jacobian[ditg].get(ig,0);
         
         // get local density and temperature and compute local tau
         numDen1 = this_numberDensity1.get(ig,0);
         numDen2 = this_numberDensity2.get(ig,0);
         if (numDen1*numDen2 == 0.0) { continue; }
         
         // compute the atomic spacing associated with minimum density species
         Real minn = std::min(numDen1,numDen2);
         Real atomic_spacing = 1.0/std::cbrt(4.0/3.0*Constants::PI*minn); // atomic spacing [m]

         // set maximum impact parameter to LDe (i.e., screening length)
         // which is already limited by atomic spacing for each species
         bmax = this_LDe.get(ig,0);

         // set max sigma based on mfp = atomic spacing
         sigma_max = 1.0/(minn*atomic_spacing);

         // set the Fermi energy for this cell
         // EF_norm = hbar^2/(2*mass)*(3*pi^2*n)^(2/3)/(me*cvac^2)
         Real maxn = std::max(numDen1,numDen2);
         m_EF_norm = m_EF_fact*std::pow(maxn,2.0/3.0);

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
         List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
         numCell1 = cell_pList1.length();
         numCell2 = cell_pList2.length();
         if (numCell1*numCell2 < 2) { continue; }
         Nmin = std::min(numCell1,numCell2);
         Nmax = std::max(numCell1,numCell2);

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
         for (lit2.begin(); lit2.ok(); ++lit2) vector_part2_ptrs.push_back(lit2());
         std::shuffle(vector_part2_ptrs.begin(),vector_part2_ptrs.end(),global_rand_gen); 
         
         unsigned int p1, p2;
         for (auto p=0; p<Nmax; p++) { // loop over particle scattering pairs
     
            if (Nmin==numCell1) {
               p1 = p % numCell1;
               p2 = p;
            } 
            else {
               p1 = p;
               p2 = p % numCell2;
            } 

            // get particle data for first particle    
            JustinsParticlePtr& part1 = vector_part1_ptrs[p1];
            part1_ptr = part1.getPointer();
            std::array<Real,3>& betap1 = part1_ptr->velocity();
            const Real wp1 = part1_ptr->weight();
            
            // get particle data for second particle    
            JustinsParticlePtr& part2 = vector_part2_ptrs[p2];
            part2_ptr = part2.getPointer();
            std::array<Real,3>& betap2 = part2_ptr->velocity();
            const Real wp2 = part2_ptr->weight();

            // determine maximum weight of this pair and set scatter flags (rejection method)
            wpMax = std::max(wp1,wp2);
            den12 = wpMax*Nmin/cellV_SI;
            
            // compute deltaU
            std::array<Real,3> deltaU;
            GalileanScatter( deltaU, betap1, betap2, den12, bmax, sigma_max, a_dt_sec );

            // update particle velocities
            if ( (float)wp1 < (float)wp2 ) {

               const long double ratio = wp1/wp2;
               betap_before = betap2;
               Ebefore = m_mass2*(betap2[0]*betap2[0] + betap2[1]*betap2[1] + betap2[2]*betap2[2])/2.0;
               betap_sum = {0.0,0.0,0.0};
               Escatter_sum = 0.0;

               // scatter particles
               std::array<Real,3> betap2p = betap2;
               for (int n=0; n<3; n++) { betap1[n]  += m_mu/m_mass1*deltaU[n]; }
               for (int n=0; n<3; n++) { betap2p[n] -= m_mu/m_mass2*deltaU[n]; }

               // update betap_sum and Escatter_sum
               for (int n=0; n<3; n++) { betap_sum[n] += betap2p[n]; }
               Escatter_sum += m_mass2*(betap2p[0]*betap2p[0] + betap2p[1]*betap2p[1] + betap2p[2]*betap2p[2])/2.0;
                     
               const Real Eafter = Ebefore + ratio*(Escatter_sum-Ebefore);
               for (int n=0; n<3; n++) { betap2[n] = betap_before[n] + ratio*(betap_sum[n]-betap_before[n]); }

               // correct beta of particle 2 to conserve energy exactly and momentum on average
               enforceEnergyConservation( betap2, m_mass2, Eafter );
                    
            }
            else if ( (float)wp2 < (float)wp1 ) {

               const long double ratio = wp2/wp1;
               betap_before = betap1;
               Ebefore = m_mass1*(betap1[0]*betap1[0] + betap1[1]*betap1[1] + betap1[2]*betap1[2])/2.0;
               betap_sum = {0.0,0.0,0.0};
               Escatter_sum = 0.0;

               // scatter particles
               std::array<Real,3> betap1p = betap1;
               for (int n=0; n<3; n++) { betap1p[n] += m_mu/m_mass1*deltaU[n]; }
               for (int n=0; n<3; n++) { betap2[n]  -= m_mu/m_mass2*deltaU[n]; }
                     
               // update betap_sum and Escatter_sum
               for (int n=0; n<3; n++) { betap_sum[n] += betap1p[n]; }
               Escatter_sum += m_mass1*(betap1p[0]*betap1p[0] + betap1p[1]*betap1p[1] + betap1p[2]*betap1p[2])/2.0;
                     
               const Real Eafter = Ebefore + ratio*(Escatter_sum-Ebefore);
               for (int n=0; n<3; n++) { betap1[n] = betap_before[n] + ratio*(betap_sum[n]-betap_before[n]); }

               // correct beta of particle 1 to conserve energy exactly and momentum on average
               enforceEnergyConservation( betap1, m_mass1, Eafter );

            }
            else {
               for (int n=0; n<3; n++) { betap1[n] += m_mu/m_mass1*deltaU[n]; }
               for (int n=0; n<3; n++) { betap2[n] -= m_mu/m_mass2*deltaU[n]; }
            }

         }
         verbosity=0;

      } // end loop over cells

   } // end loop over boxes

   // don't forget to set pointers back to NULL and delete
   part1_ptr = NULL;
   part2_ptr = NULL;
   delete part1_ptr;
   delete part2_ptr;
   
}

void Coulomb::GalileanScatter( std::array<Real,3>&  a_deltaU,
                         const std::array<Real,3>&  a_vp1,
                         const std::array<Real,3>&  a_vp2,
                         const Real                 a_den12,
                         const Real                 a_bmax,
                         const Real                 a_sigma_max,
                         const Real                 a_dt_sec ) const
{
   CH_TIME("Coulomb::GalileanScatter()");
   
   // compute relative velocity (actually beta)
   Real ux = a_vp1[0]-a_vp2[0];
   Real uy = a_vp1[1]-a_vp2[1];
   Real uz = a_vp1[2]-a_vp2[2];
   Real u = sqrt(ux*ux + uy*uy + uz*uz);

   // don't proceed if vrel = 0 or if vrel is relatively close to machine round-off
   if (u<=std::numeric_limits<Real>::min()) { return; } // avoid divide by zero
   const Real vsum = std::sqrt(a_vp1[0]*a_vp1[0] + a_vp1[1]*a_vp1[1] + a_vp1[2]*a_vp1[2])
                   + std::sqrt(a_vp2[0]*a_vp2[0] + a_vp2[1]*a_vp2[1] + a_vp2[2]*a_vp2[2]);
   if (u<=1.0e-14*vsum) { return; }

   // compute s12 = 2.0*deltasq_var
   Real b0 = m_b90_fact/(m_mu*u*u + 2.0*m_EF_norm); // b0 = q1*q2/(4*pi*ep0*WR) = q1*q2/(2*pi*ep0*mu*u^2)
   const Real bmin_qm = m_bqm_fact/(m_mu*u + std::sqrt(2.0*m_EF_norm*m_mu));
   const Real bmin = std::max(b0/2.0,bmin_qm);
   Real Clog = m_Clog;
   if (Clog==0.0 && u>0.0) {
      Clog = 0.5*std::log(1.0 + a_bmax*a_bmax/bmin/bmin);
      Clog = std::max(2.0,Clog);
   }

   // compute s12 = 2*<delta^2> with sigma_eff limited by sigma_max where mfp = atomic spacing
   b0 = m_b90_fact/(m_mu*u*u); // don't use Fermi energy in b0 here
   Real sigma_eff = Constants::PI*b0*b0*Clog;
   sigma_eff = std::min(sigma_eff,a_sigma_max);
   m_s12 = sigma_eff * a_den12 * u*Constants::CVAC * a_dt_sec;

   switch (m_angular_scattering) {
      case TAKIZUKA:
         if (m_s12<2.0) { // sample from gaussian distribution  
            const Real delta = sqrt(m_s12/2.0)*std::abs(MathUtils::randn());
            const Real deltasq = delta*delta;
            m_sinth = 2.0*delta/(1.0+deltasq); 
            m_costh = 1.0 - 2.0*deltasq/(1.0+deltasq);
         } 
         else { // set random polar angle between zero and pi
            const Real theta = Constants::PI*MathUtils::rand();
            m_costh = cos(theta);
            m_sinth = sin(theta);
         }
         break;
      case NANBU:
         setNANBUcosthsinth(m_s12);
         break;
      case NANBU_FAS:
         setNANBUFAScosthsinth(m_s12, Clog, b0, bmin_qm, a_bmax, sigma_eff);
         break;
      case NANBU_FAS_v2:
         setNANBUFAS_v2_costhsinth(m_s12, Clog, b0, bmin_qm, a_bmax, sigma_eff);
         break;
      case ISOTROPIC:
         const Real theta = Constants::PI*MathUtils::rand();
         m_costh = cos(theta);
         m_sinth = sin(theta);
         break;
   }

   // set random azimuthal angle
   m_phi = Constants::TWOPI*MathUtils::rand();
   m_cosphi = cos(m_phi);
   m_sinphi = sin(m_phi);

   // define deltaU
   ScatteringUtils::computeDeltaU(a_deltaU,ux,uy,uz,m_costh,m_sinth,m_cosphi,m_sinphi);
               
}

void Coulomb::LorentzScatter( std::array<Real,3>&  a_up1,
                              std::array<Real,3>&  a_up2,
                        const bool                 a_scatter2,
                        const long double          a_mass1,
                        const long double          a_mass2,
                        const Real                 a_den12,
                        const Real                 a_bmax,
                        const Real                 a_sigma_max,
                        const Real                 a_dt_sec ) const
{
   CH_TIME("Coulomb::LorentzScatter()");
   
   // Relativistic scattering method for two particles with weight2 > weight1

   long double g1, g2, vcmsq, Etot, gcm, g1st, g2st;
   long double ucmdotup, vrelst, vrelst_invar, muRst, upst_fact, upstsq, denom;
   std::array<Real,3> ptot, vcm, upst;

   // compute the lab frame total energy and momentum
   Real gb1sq = a_up1[0]*a_up1[0] + a_up1[1]*a_up1[1] + a_up1[2]*a_up1[2];
   Real gb2sq = a_up2[0]*a_up2[0] + a_up2[1]*a_up2[1] + a_up2[2]*a_up2[2];
   g1 = sqrt(1.0 + gb1sq);
   g2 = sqrt(1.0 + gb2sq);
   Etot = g1*a_mass1 + g2*a_mass2;
   for (int n=0; n<3; n++) { ptot[n] = a_mass1*a_up1[n] + a_mass2*a_up2[n]; }
 
   // compute center-of-momentum velocity (beta) and gamma
   for (int n=0; n<3; n++) { vcm[n] = ptot[n]/Etot; }
   vcmsq = vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2];
   gcm = 1.0/std::sqrt(1.0 - vcmsq);

   // compute gamma* and up* for particles 1 and 2 ( * = CM frame )
   ucmdotup = gcm*(vcm[0]*a_up2[0] + vcm[1]*a_up2[1] + vcm[2]*a_up2[2]);
   g2st = gcm*g2 - ucmdotup;

   ucmdotup = gcm*(vcm[0]*a_up1[0] + vcm[1]*a_up1[1] + vcm[2]*a_up1[2]);
   g1st = gcm*g1 - ucmdotup;
 
   // compute up* = p1*/m1 = -m2/m2*p2
   upst_fact = gcm*(ucmdotup/(1.0+gcm) - g1);
   for (int n=0; n<3; n++) { upst[n] = a_up1[n] + upst_fact*vcm[n]; }
 
   // compute non-invariant relative velocity in CM frame |v1* - v2*| = |p1*|/muR*
   muRst = g1st*a_mass1*g2st*a_mass2/(g1st*a_mass1 + g2st*a_mass2);
   upstsq = upst[0]*upst[0] + upst[1]*upst[1] + upst[2]*upst[2];
   vrelst = std::sqrt(upstsq)*a_mass1/muRst;

   // don't proceed if vrel = 0 or if vrel is relatively close to machine round-off
   if (vrelst<=std::numeric_limits<Real>::min()) { return; } // avoid divide by zero
   const Real vsum = std::sqrt(gb1sq)/g1 + std::sqrt(gb2sq)/g2;
   if (vrelst<=1.0e-14*vsum) { return; }

   // compute invariant relative velocity in CM frame |v1* - v2*|/(1 - v1*cdotv2*)
   denom = 1.0 + upstsq*a_mass1/a_mass2/g1st/g2st;
   vrelst_invar = vrelst/denom;
   
   Real b0 = m_b90_fact/(muRst*vrelst*vrelst_invar + 2.0*m_EF_norm);
   const Real bmin_qm = m_bqm_fact/(muRst*vrelst + std::sqrt(2.0*m_EF_norm*muRst));
   const Real bmin = std::max(b0/2.0,bmin_qm);
   Real Clog = m_Clog;
   if (Clog==0.0 && upstsq>0.0) {
      Clog = 0.5*std::log(1.0 + a_bmax*a_bmax/bmin/bmin);
      Clog = std::max(2.0,Clog);
   }

   // compute s12 = 2*<delta^2> with sigma_eff limited by sigma_max where mfp = atomic spacing
   b0 = m_b90_fact/(muRst*vrelst*vrelst_invar); // dont use Fermi energy for b0 here
   Real sigma_eff = Constants::PI*b0*b0*Clog;
   sigma_eff = std::min(sigma_eff,a_sigma_max);
   m_s12 = sigma_eff * a_den12 * vrelst*Constants::CVAC * a_dt_sec;
   m_s12 *= g1st*g2st/g1/g2;

   switch (m_angular_scattering) {
      case TAKIZUKA:
         if (m_s12<2.0) { // sample from gaussian distribution  
            const Real delta = sqrt(m_s12/2.0)*std::abs(MathUtils::randn());
            const Real deltasq = delta*delta;
            m_sinth = 2.0*delta/(1.0+deltasq); 
            m_costh = 1.0 - 2.0*deltasq/(1.0+deltasq);
         } 
         else { // set random polar angle between zero and pi
            const Real theta = Constants::PI*MathUtils::rand();
            m_costh = cos(theta);
            m_sinth = sin(theta);
         }
         break;
      case NANBU:
         setNANBUcosthsinth(m_s12);
         break;
      case NANBU_FAS:
         setNANBUFAScosthsinth(m_s12, Clog, b0, bmin_qm, a_bmax, sigma_eff);
         break;
      case NANBU_FAS_v2:
         setNANBUFAS_v2_costhsinth(m_s12, Clog, b0, bmin_qm, a_bmax, sigma_eff);
         break;
      case ISOTROPIC:
         const Real theta = Constants::PI*MathUtils::rand();
         m_costh = cos(theta);
         m_sinth = sin(theta);
         break;
   }

   // set random azimuthal angle
   m_phi = Constants::TWOPI*MathUtils::rand();
   m_cosphi = cos(m_phi);
   m_sinphi = sin(m_phi);

   // rotate upst for particle 1 by scattering angles
   ScatteringUtils::rotateVelocity(upst,m_costh,m_sinth,m_cosphi,m_sinphi);

   // Lorentz transform particle 1 back to lab frame
   ucmdotup = gcm*(vcm[0]*upst[0] + vcm[1]*upst[1] + vcm[2]*upst[2]);
   upst_fact = gcm*(ucmdotup/(1.0+gcm) + g1st);
   for (int n=0; n<3; n++) { a_up1[n] = upst[n] + upst_fact*vcm[n]; }
  
   // Lorentz transform particle 2 back to lab frame (p2st + p1st = 0)
   if (a_scatter2) {
      //for (int n=0; n<3; n++) { upst[n] *= -a_mass1/a_mass2; }
      //ucmdotup *= -a_mass1/a_mass2;
      //upst_fact = gcm*(ucmdotup/(1.0+gcm) + g2st);
      //for (int n=0; n<3; n++) { a_up2[n] = upst[n] + upst_fact*vcm[n]; }
      // set proper velocity of particle 2 using momentum conservation
      for (int n=0; n<3; n++) { a_up2[n] = (ptot[n] - a_mass1*a_up1[n])/a_mass2; }
   }
 
}

#include "NamespaceFooter.H"

