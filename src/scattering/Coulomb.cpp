
#include "Coulomb.H"
#include "MathUtils.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"

#include "NamespaceHeader.H"


void Coulomb::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                          const DomainGrid&           a_mesh )
{
   CH_TIME("Coulomb::initialize()");
   
   const PicSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getPtrVect();
   
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
   
   if (m_verbosity)  printParameters();

}

void Coulomb::setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const 
{
   CH_TIME("Coulomb::setMeanFreeTime()");
   
   const LevelData<FArrayBox>& LDe_m = a_pic_species_intf.getDebyeLength();
   
   const PicSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getPtrVect();
   PicSpeciesPtr this_species1(pic_species_ptr_vect[m_sp1]);
   PicSpeciesPtr this_species2(pic_species_ptr_vect[m_sp2]);
   
   if(!this_species1->scatter() || !this_species2->scatter()) return;
   
   const bool setMoments = false;
   const LevelData<FArrayBox>& numDen1 = this_species1->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& momDen1 = this_species1->getMomentumDensity(setMoments);
   const LevelData<FArrayBox>& eneDen1 = this_species1->getEnergyDensity(setMoments);
   
   if(m_sp1==m_sp2) setIntraMFT(LDe_m,numDen1,momDen1,eneDen1);
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
   Real b90, sigma90, nu90;
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
         if(numDen == 0.0) continue;
	 rho = m_mass1*numDen;
	 
	 // set maximum impact parameter
	 Real LDe = this_LDe.get(ig,0);
	 Real atomic_spacing = std::pow(4.189*numDen,-1.0/3.0); // atomic spacing [m]
	 bmax = std::max(LDe,atomic_spacing);
	    
	 // compute Temperature and g12sq
	 Real rhoUx = this_momDensity.get(ig,0); // mass1*<betax*numDen>
	 Real rhoUy = this_momDensity.get(ig,1); // mass1*<betay*numDen>
	 Real rhoUz = this_momDensity.get(ig,2); // mass1*<betaz*numDen>
	 meanE = (rhoUx*rhoUx + rhoUy*rhoUy + rhoUz*rhoUz)/rho/2.0;
         
         eneDen = 0.0;
         for(int dir=0; dir<3; dir++) eneDen += this_eneDensity.get(ig,dir);  
         T_eV = 2.0/3.0*(eneDen - meanE)/numDen*mcSq_eV;
	 T_eV = std::max(T_eV,0.01);
	 
	 // compute g12^2 = 3*T1/m1 + 3*T2/m2 + |U1 - U2|^2 and b90
         g12sq = 6.0*Constants::QE/Constants::ME*T_eV/m_mass1; // [m^2/s^2]
         g12sq_norm = g12sq/cvacSq; // normalized
         b90 = m_b90_fact/(m_mu*g12sq_norm); // [m]

	 // set the Coulomb logarithm
	 Real Clog = m_Clog;
	 if(Clog==0.0 && g12sq>0.0) { // See Lee and More 1984 Eqs 20-22
	    Real bmin_qm = m_bqm_fact/(m_mu*std::sqrt(g12sq_norm));
            bmin = std::max(b90/2.0,bmin_qm);
            Clog = 0.5*std::log(1.0 + bmax*bmax/bmin/bmin);
            Clog = std::max(2.0,Clog);
	 }

         // define self-species collision time
	 sigma90 = 8.0/Constants::PI*b90*b90*Clog; // [m^2]
         nu90 = sqrt(g12sq)*numDen*sigma90;        // [Hz]
         box_nuMax = Max(box_nuMax,nu90);
            
	 bool verbose = false;
	 if(!procID() && verbose) {
            cout << "JRA: sp1,sp2 = " << m_sp1 << ", " << m_sp2 << endl;
            cout << "JRA: mass1 = " << m_mass1 << endl;
            cout << "JRA: mass2 = " << m_mass2 << endl;
            cout << "JRA: N[ig="<<ig<<"]  = " << numDen <<" [1/m^3]" << endl;
            cout << "JRA: T[ig="<<ig<<"]  = " << T_eV <<" [eV]" << endl;
            cout << "JRA: LDe[ig="<<ig<<"]   = " << LDe  <<" [m]" << endl;
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
   Real b90, sigma90, nu90;
         
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
         if(numDen1*numDen2 == 0.0) continue;
	 rho1 = m_mass1*numDen1;
	 rho2 = m_mass2*numDen2;

	 // set maximum impact parameter
	 Real LDe = this_LDe.get(ig,0);
	 Real Nmax = std::max(numDen1,numDen2);
	 Real atomic_spacing = std::pow(4.189*Nmax,-1.0/3.0); // atomic spacing [m]
	 bmax = std::max(LDe,atomic_spacing);
	 
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
         for( int dir=0; dir<3; dir++) eneDen1 += this_eneDensity1.get(ig,dir);  
         for( int dir=0; dir<3; dir++) eneDen2 += this_eneDensity2.get(ig,dir);  
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
         b90 = m_b90_fact/(m_mu*g12sq_norm); // [m]
         
	 // set the Coulomb logarithm
	 Real Clog = m_Clog;
	 if(Clog==0.0 && g12sq>0.0) { // See Lee and More 1984 Eqs 20-22
	    Real bmin_qm = m_bqm_fact/(m_mu*std::sqrt(g12sq_norm));
            bmin = std::max(b90/2.0,bmin_qm);
            Clog = 0.5*std::log(1.0 + bmax*bmax/bmin/bmin);
            Clog = std::max(2.0,Clog);
	 }
         
	 // define intra-species collision time
	 sigma90 = 8.0/Constants::PI*b90*b90*Clog; // [m^2]
         nu90 = sqrt(g12sq)*Nmax*sigma90;        // [Hz]
         box_nuMax = Max(box_nuMax,nu90);
         
         bool verbose = false;
         if(!procID() && verbose) {
            cout << "JRA: sp1,sp2 = " << m_sp1 << ", " << m_sp2 << endl;
            cout << "JRA: mass1 = " << m_mass1 << endl;
            cout << "JRA: mass2 = " << m_mass2 << endl;
            cout << "JRA: N1[ig="<<ig<<"]  = " << numDen1 <<" [1/m^3]" << endl;
            cout << "JRA: N2[ig="<<ig<<"]  = " << numDen2 <<" [1/m^3]" << endl;
            cout << "JRA: T1[ig="<<ig<<"]  = " << T1_eV <<" [eV]" << endl;
            cout << "JRA: T2[ig="<<ig<<"]  = " << T2_eV <<" [eV]" << endl;
            cout << "JRA: LDe[ig="<<ig<<"]   = " << LDe  <<" [m]" << endl;
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
   if(!this_species1->scatter()) return;
  
   const LevelData<FArrayBox>& LDe_m = a_pic_species_intf.getDebyeLength();

   if(m_sp1==m_sp2) {
     switch (m_weight_method) {
       case PROBABILISTIC: // scatter larger weight particles probabilistically
//#define ALT_METHOD1
#ifdef ALT_METHOD1
       // scatters all particles with own normalized mfp
       applyIntraScattering_PROB_alt( *this_species1, a_mesh, LDe_m, a_dt_sec );
#else
       applyIntraScattering_PROB( *this_species1, a_mesh, LDe_m, a_dt_sec );
#endif
       break;
       case CONSERVATIVE: // sentoku and kemp JCP 2008
       applyIntraScattering_SK08( *this_species1, a_mesh, LDe_m, a_dt_sec );
       break;
     }
   }
   else {
     PicSpeciesPtr this_species2(pic_species_ptr_vect[m_sp2]);
     if(!this_species2->scatter()) return;
     switch (m_weight_method) {
       case PROBABILISTIC: // scatter larger weight particles probabilistically
       applyInterScattering_PROB( *this_species1, *this_species2, a_mesh, LDe_m, a_dt_sec );
       break;
       case CONSERVATIVE: // sentoku and kemp JCP 2008
       applyInterScattering_SK08( *this_species1, *this_species2, a_mesh, LDe_m, a_dt_sec );
       break;
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
   int numCell, Naa, Nsubcycle=1;
   Real numDen, den12, bmax, cellV_SI, wpMax, wpMin;
 
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
       
         // get local density and set bmax
         numDen = this_numberDensity.get(ig,0);
         if(numDen == 0.0) continue;
	 
	 // set maximum impact parameter
	 Real LDe = this_LDe.get(ig,0);
	 Real atomic_spacing = std::pow(4.189*numDen,-1.0/3.0); // atomic spacing [m]
	 bmax = std::max(LDe,atomic_spacing);
	 
         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         numCell = cell_pList.length();
         if(numCell < 2) continue;
          
         // define flag needed to properly treat the first few collisions when numCell is odd
         bool odd_NxN = true;        
         if(numCell % 2 == 0) odd_NxN = false;
   
         Naa = numCell - 1;
    
         // copy the iterators to a vector in order to shuffle
         vector_part_ptrs.clear();
         vector_part_ptrs.reserve(numCell);
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) vector_part_ptrs.push_back(lit());
         std::shuffle(vector_part_ptrs.begin(),vector_part_ptrs.end(),global_rand_gen); 
	 
	 Real Wtot0 = 0.0;
	 Real Etot0 = 0.0;
	 if(m_enforce_conservations) { // save the initial mass and Energy
            for (lit.begin(); lit.ok(); ++lit) {
               part1_ptr = lit().getPointer();
               const std::array<Real,3>& betap = part1_ptr->velocity();
               const Real& wp = part1_ptr->weight();
	       Wtot0 += wp;
	       for (int n=0; n<3; n++) Etot0 += wp/2.0*betap[n]*betap[n];
	    }
	 }

	 std::array<Real,3> dBetaAvg = {0.0,0.0,0.0};
         for (int p=0; p<vector_part_ptrs.size(); p++) { // loop over particle scattering pairs
            
            int p1, p2; 
	    if(odd_NxN) {
               p1 = p % 2;
               if(p==0) p2 = 1;
	       else if(p==1) p2 = 2;
	       else if(p==2) p2 = 2;
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
	    if(odd_NxN) {
	       den12 = den12/2.0;
	       if(p==2) odd_NxN = false;
	    }
	    
	    if(m_Nsubcycle_max>1) {
               wpMin = std::min(wp1,wp2);
      	       Nsubcycle = std::round(wpMax/wpMin);
      	       Nsubcycle = std::min(Nsubcycle,m_Nsubcycle_max);
               den12 = den12/Nsubcycle;
	    }

	    for (int m=0; m<Nsubcycle; m++) {

#ifdef RELATIVISTIC_PARTICLES   
               if( (float)wp2 < (float)wp1 ) {
                  bool scatter1 = true;
                  if(MathUtils::rand()>wp2/wp1) scatter1 = false; 
                  LorentzScatter( betap2, betap1, scatter1, m_mass2, m_mass1, den12, bmax, a_dt_sec );
               }
               else {
                  bool scatter2 = true;
		  if( (float)wp1 < (float)wp2 ) {if(MathUtils::rand()>wp1/wp2) scatter2 = false;}
                  LorentzScatter( betap1, betap2, scatter2, m_mass1, m_mass2, den12, bmax, a_dt_sec );
               }
#else
               // compute deltaU
               std::array<Real,3> deltaU;
	       GalileanScatter( deltaU, betap1, betap2, den12, bmax, a_dt_sec );

               // update particle velocities
               if( (float)wp1 == (float)wp2 ) {
		  for (int n=0; n<3; n++) betap1[n] += 0.5*deltaU[n];
                  for (int n=0; n<3; n++) betap2[n] -= 0.5*deltaU[n];
	       }
	       else if( (float)wp1 < (float)wp2 ) {
		  for (int n=0; n<3; n++) betap1[n] += 0.5*deltaU[n];
	          for (int n=0; n<3; n++) dBetaAvg[n] += wp1*0.5*deltaU[n];
		  if(MathUtils::rand()<wp1/wp2) {
                     for (int n=0; n<3; n++) betap2[n] -= 0.5*deltaU[n];
	             for (int n=0; n<3; n++) dBetaAvg[n] -= wp2*0.5*deltaU[n];
		  }
	       }
               else { //if( (float)wp2 < (float)wp1 ) {
		  if(MathUtils::rand()<wp2/wp1) {
		     for (int n=0; n<3; n++) betap1[n] += 0.5*deltaU[n];
	             for (int n=0; n<3; n++) dBetaAvg[n] += wp1*0.5*deltaU[n];
		  }
                  for (int n=0; n<3; n++) betap2[n] -= 0.5*deltaU[n];
	          for (int n=0; n<3; n++) dBetaAvg[n] -= wp2*0.5*deltaU[n];
	       }

	    }
#endif
         }
         verbosity=0;
	 
         // adjust particle betas such that momentum and energy are identically conserved
	 Real dBetaSq = 0.0;
	 for (int n=0; n<3; n++) dBetaSq += dBetaAvg[n]*dBetaAvg[n];
	 if(dBetaSq>0.0 && m_enforce_conservations) {

	    // enforce momentum conservation and compute deltaE
	    for (int n=0; n<3; n++) dBetaAvg[n] /= Wtot0;
	    Real Etot1 = 0.0;
	    for (lit.begin(); lit.ok(); ++lit) {
               part1_ptr = lit().getPointer();
               std::array<Real,3>& betap = part1_ptr->velocity();
               const Real& wp = part1_ptr->weight();
	       for (int n=0; n<3; n++) betap[n] -= dBetaAvg[n];
	       for (int n=0; n<3; n++) Etot1 += wp/2.0*betap[n]*betap[n];
	    }
	    long double deltaE = m_mass1*(Etot1 - Etot0);

            // sort particles by weight to bias adjusting particles below to
	    // heavier weight ones and to get pairs that have similar weights
	    if(m_sort_weighted_particles) {
	       std::sort( vector_part_ptrs.begin(), vector_part_ptrs.end(),
	                  []( JustinsParticlePtr& lhs, JustinsParticlePtr& rhs )
	     	          { return lhs.getPointer()->weight() > rhs.getPointer()->weight(); } );
	    }

            int count = 0;
	    int loop_count = 0;
            for (int p=0; p<vector_part_ptrs.size(); p++) {
               int p1 = p;
   	       p++;
	       if(p==vector_part_ptrs.size()) {
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
               modEnergyPairwise(betap1, betap2, m_mass1*wp1, m_mass2*wp2, deltaE);
	       count += 1;
	       if(deltaE==0.0) {
		  //cout << "JRA: ig = " << ig << endl;
		  //cout << "JRA: count = " << count << endl;
		  break;
	       }
	       if(p==vector_part_ptrs.size()-1) {
		  loop_count++;
		  p = -1;
	       }
	       if(loop_count>m_loop_count_max) {
		  cout << "Notice: loop_count > loop_count_max = " << m_loop_count_max << endl;
		  cout << "        for species1 = " << m_species1_name << endl;
		  cout << "        and species2 = " << m_species2_name << endl;
		  cout << "        num_parts    = " << vector_part_ptrs.size() << endl;
		  cout << "        ig           = " << ig << endl;
		  cout << "        Etot0        = " << Etot0 << endl;
		  cout << "        deltaE       = " << deltaE << endl;
		  break;
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

void Coulomb::applyIntraScattering_PROB_alt( PicSpecies&        a_picSpecies, 
                                   const DomainGrid&            a_mesh,
		                   const LevelData<FArrayBox>&  a_LDe_m,
                                   const Real                   a_dt_sec ) const
{
   CH_TIME("Coulomb::applyIntraScattering_PROB_alt()");
	       
   // alternative probabilistic approach where both particles always scatter
   // seems to work ok, but have observed a small net heating associated with
   // net momentum gain for test_T1b. This issue goes away if a sufficient number
   // of subcycles are used such that s12 is approx s21
 
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
   int numCell, Naa, Nsubcycle=1;
   Real numDen, den12, bmax, cellV_SI, wpMax, wpMin;
 
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
       
         // get local density  and set bmax
         numDen = this_numberDensity.get(ig,0);
         if(numDen == 0.0) continue;

	 // set maximum impact parameter
	 Real LDe = this_LDe.get(ig,0);
	 Real atomic_spacing = std::pow(4.189*numDen,-1.0/3.0); // atomic spacing [m]
	 bmax = std::max(LDe,atomic_spacing);

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         numCell = cell_pList.length();
         if(numCell < 2) continue;
          
         // define flag needed to properly treat the first few collisions when numCell is odd
         bool odd_NxN = true;        
         if(numCell % 2 == 0) odd_NxN = false;
   
         Naa = numCell - 1;
    
         // copy the iterators to a vector in order to shuffle
         vector_part_ptrs.clear();
         vector_part_ptrs.reserve(numCell);
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) vector_part_ptrs.push_back(lit());
         std::shuffle(vector_part_ptrs.begin(),vector_part_ptrs.end(),global_rand_gen); 
         
	 Real Wtot0 = 0.0;
	 Real Etot0 = 0.0;
	 if(m_enforce_conservations) { // save the initial mass and Energy
            for (lit.begin(); lit.ok(); ++lit) {
               part1_ptr = lit().getPointer();
               const std::array<Real,3>& betap = part1_ptr->velocity();
               const Real& wp = part1_ptr->weight();
	       Wtot0 += wp;
	       for (int n=0; n<3; n++) Etot0 += wp/2.0*betap[n]*betap[n];
	    }
	 }

	 std::array<Real,3> dBetaAvg = {0.0,0.0,0.0};
         for (int p=0; p<vector_part_ptrs.size(); p++) { // loop over particle scattering pairs
            
            int p1, p2; 
	    if(odd_NxN) {
               p1 = p % 2;
               if(p==0) p2 = 1;
	       else if(p==1) p2 = 2;
	       else if(p==2) p2 = 2;
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
	    if(odd_NxN) {
	       den12 = den12/2.0;
	       if(p==2) odd_NxN = false;
	    }

	    if(m_Nsubcycle_max>1) {
               wpMin = std::min(wp1,wp2);
      	       Nsubcycle = std::floor(wpMax/wpMin)-1;
	       if(Nsubcycle==0) Nsubcycle = 1;
	       else Nsubcycle = std::min(Nsubcycle,m_Nsubcycle_max);
	    }
	    
	    // alternative prob approach where both particles always scatter, seems to work ok, but 
	    // have observed a small net heating for test_T1b and test_T1c if Nsubcycle is too small
	    if( (float)wp1 < (float)wp2 ) {
               std::array<Real,3> deltaU;
	       for (int m=0; m<Nsubcycle; m++) {
	          GalileanScatter( deltaU, betap1, betap2, den12*(1.0-wp1/wp2)/Nsubcycle, bmax, a_dt_sec );
	          for (int n=0; n<3; n++) betap1[n] += 0.5*deltaU[n];
	          for (int n=0; n<3; n++) dBetaAvg[n] += wp1*0.5*deltaU[n];
	       }
	       GalileanScatter( deltaU, betap1, betap2, den12*wp1/wp2, bmax, a_dt_sec );
	       for (int n=0; n<3; n++) betap1[n] += 0.5*deltaU[n];
	       for (int n=0; n<3; n++) betap2[n] -= 0.5*deltaU[n];
	       for (int n=0; n<3; n++) {
	          dBetaAvg[n] += 0.5*deltaU[n]*(wp1-wp2);
	       }
	    }
	    else if( (float)wp2 < (float)wp1 ) { 
               std::array<Real,3> deltaU;
	       for (int m=0; m<Nsubcycle; m++) {
	          GalileanScatter( deltaU, betap1, betap2, den12*(1.0-wp2/wp1)/Nsubcycle, bmax, a_dt_sec );
	          for (int n=0; n<3; n++) betap2[n] -= 0.5*deltaU[n];
	          for (int n=0; n<3; n++) dBetaAvg[n] -= wp2*0.5*deltaU[n];
	       }
	       GalileanScatter( deltaU, betap1, betap2, den12*wp2/wp1, bmax, a_dt_sec );
	       for (int n=0; n<3; n++) betap1[n] += 0.5*deltaU[n];
	       for (int n=0; n<3; n++) betap2[n] -= 0.5*deltaU[n];
	       for (int n=0; n<3; n++) {
	          dBetaAvg[n] += 0.5*deltaU[n]*(wp1-wp2);
	       }
	    }
	    else { //( (float)wp1 == (float)wp2 ) {
               std::array<Real,3> deltaU;
	       GalileanScatter( deltaU, betap1, betap2, den12, bmax, a_dt_sec );
               for (int n=0; n<3; n++) betap1[n] += 0.5*deltaU[n];
               for (int n=0; n<3; n++) betap2[n] -= 0.5*deltaU[n];
	    }

         }
         verbosity=0;
	 
         // adjust particle betas such that momentum and energy are identically conserved
	 Real dBetaSq = 0.0;
	 for (int n=0; n<3; n++) dBetaSq += dBetaAvg[n]*dBetaAvg[n];
	 if(dBetaSq>0.0 && m_enforce_conservations) {

	    // enforce momentum conservation and compute deltaE
	    for (int n=0; n<3; n++) dBetaAvg[n] /= Wtot0;
	    Real Etot1 = 0.0;
	    for (lit.begin(); lit.ok(); ++lit) {
               part1_ptr = lit().getPointer();
               std::array<Real,3>& betap = part1_ptr->velocity();
               const Real& wp = part1_ptr->weight();
	       for (int n=0; n<3; n++) betap[n] -= dBetaAvg[n];
	       for (int n=0; n<3; n++) Etot1 += wp/2.0*betap[n]*betap[n];
	    }
	    long double deltaE = m_mass1*(Etot1 - Etot0);

            // sort particles by weight to bias adjusting particles below to
	    // heavier weight ones and to get pairs that have similar weights
	    if(m_sort_weighted_particles) {
	       std::sort( vector_part_ptrs.begin(), vector_part_ptrs.end(),
	                  []( JustinsParticlePtr& lhs, JustinsParticlePtr& rhs )
	    	          { return lhs.getPointer()->weight() > rhs.getPointer()->weight(); } );
	    }

            int count = 0;            
            for (int p=0; p<vector_part_ptrs.size(); p++) {
               int p1 = p;
   	       p++;
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
               modEnergyPairwise(betap1, betap2, m_mass1*wp1, m_mass2*wp2, deltaE);
	       count += 1;
	       if(deltaE==0.0) {
		  //cout << "JRA: ig = " << ig << endl;
		  //cout << "JRA: count = " << count << endl;
		  break;
	       }
	       if(p==vector_part_ptrs.size()-1) p = -1;
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
   int numCell, Naa, Nsubcycle=1;
   Real numDen, den12, bmax, cellV_SI, wpMax, wpMin, Escatter_sum, Ebefore;
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
         if(numDen == 0.0) continue;
	 
	 // set maximum impact parameter
	 Real LDe = this_LDe.get(ig,0);
	 Real atomic_spacing = std::pow(4.189*numDen,-1.0/3.0); // atomic spacing [m]
	 bmax = std::max(LDe,atomic_spacing);

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         numCell = cell_pList.length();
         if(numCell < 2) continue;
          
         // define flag needed to properly treat the first few collisions when numCell is odd
         bool odd_NxN = true;        
         if(numCell % 2 == 0) odd_NxN = false;
   
         Naa = numCell - 1;
    
         // copy the iterators to a vector in order to shuffle
         vector_part_ptrs.clear();
         vector_part_ptrs.reserve(numCell);
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) vector_part_ptrs.push_back(lit());
         std::shuffle(vector_part_ptrs.begin(),vector_part_ptrs.end(),global_rand_gen); 
         
         for (int p=0; p<vector_part_ptrs.size(); p++) { // loop over particle scattering pairs
            
            int p1, p2; 
	    if(odd_NxN) {
               p1 = p % 2;
               if(p==0) p2 = 1;
	       else if(p==1) p2 = 2;
	       else if(p==2) p2 = 2;
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
	    if(odd_NxN) {
	       den12 = den12/2.0;
	       if(p==2) odd_NxN = false;
	    }
	    
	    if(m_Nsubcycle_max>1) {
               wpMin = std::min(wp1,wp2);
      	       Nsubcycle = std::floor(wpMax/wpMin);
      	       Nsubcycle = std::min(Nsubcycle,m_Nsubcycle_max);
               den12 = den12/Nsubcycle;
	    }

	    for (int m=0; m<Nsubcycle; m++) {

               // compute deltaU
               std::array<Real,3> deltaU;
	       GalileanScatter( deltaU, betap1, betap2, den12, bmax, a_dt_sec );

               // update particle velocities
	       if( (float)wp1 < (float)wp2 ) {

	          const long double ratio = wp1/wp2;
		  if(m==0) { // save betap and energy for particle 2 before scatter update
                     betap_before = betap2;               
		     Ebefore = m_mass2*(betap2[0]*betap2[0] + betap2[1]*betap2[1] + betap2[2]*betap2[2])/2.0;
		     betap_sum = {0.0,0.0,0.0};
		     Escatter_sum = 0.0;
		  }

		  // scatter particles
		  std::array<Real,3> betap2p = betap2;
                  for (int n=0; n<3; n++) betap1[n]  += 0.5*deltaU[n];
		  for (int n=0; n<3; n++) betap2p[n] -= 0.5*deltaU[n];

		  // update betap_sum and Escatter_sum
		  for (int n=0; n<3; n++) betap_sum[n] += betap2p[n];
		  Escatter_sum += m_mass2*(betap2p[0]*betap2p[0] + betap2p[1]*betap2p[1] + betap2p[2]*betap2p[2])/2.0;
		     
		  if(m==Nsubcycle-1) {
		     const Real Eafter = Ebefore + ratio*(Escatter_sum-Nsubcycle*Ebefore);
		     for (int n=0; n<3; n++) betap2[n] = betap_before[n] + ratio*(betap_sum[n]-Nsubcycle*betap_before[n]);

		     // correct beta of particle 2 to conserve energy exactly and momentum on average
		     enforceEnergyConservation( betap2, m_mass2, Eafter );
		  }

	       }
	       else if( (float)wp2 < (float)wp1 ) {

	          const long double ratio = wp2/wp1;
	          if(m==0) { // save betap and energy for particle 1 before scatter update
	             betap_before = betap1;
		     Ebefore = m_mass1*(betap1[0]*betap1[0] + betap1[1]*betap1[1] + betap1[2]*betap1[2])/2.0;
		     betap_sum = {0.0,0.0,0.0};
		     Escatter_sum = 0.0;
                  }

		  // scatter particles
		  std::array<Real,3> betap1p = betap1;
                  for (int n=0; n<3; n++) betap1p[n] += 0.5*deltaU[n];
		  for (int n=0; n<3; n++) betap2[n]  -= 0.5*deltaU[n];
		     
		  // update betap_sum and Escatter_sum
		  for (int n=0; n<3; n++) betap_sum[n] += betap1p[n];
		  Escatter_sum += m_mass1*(betap1p[0]*betap1p[0] + betap1p[1]*betap1p[1] + betap1p[2]*betap1p[2])/2.0;
		     
		  if(m==Nsubcycle-1) {
		     const Real Eafter = Ebefore + ratio*(Escatter_sum-Nsubcycle*Ebefore);
		     for (int n=0; n<3; n++) betap1[n] = betap_before[n] + ratio*(betap_sum[n]-Nsubcycle*betap_before[n]);

   		     // correct beta of particle 1 to conserve energy exactly and momentum on average
		     enforceEnergyConservation( betap1, m_mass1, Eafter );
		  }

	       }
	       else {
	          for (int n=0; n<3; n++) betap1[n] += 0.5*deltaU[n];
                  for (int n=0; n<3; n++) betap2[n] -= 0.5*deltaU[n];
	       }

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
   int numCell1, numCell2, pMax, pMin, Nsubcycle=1;
   Real numDen1, wpMax, wpMin;
   Real numDen2, den12, bmax, cellV_SI;
 
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
         if(numDen1*numDen2 == 0.0) continue;
	 
	 // set maximum impact parameter
	 Real LDe = this_LDe.get(ig,0);
	 Real MaxN = std::max(numDen1,numDen2);
	 Real atomic_spacing = std::pow(4.189*MaxN,-1.0/3.0); // atomic spacing [m]
	 bmax = std::max(LDe,atomic_spacing);

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
         List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
         numCell1 = cell_pList1.length();
         numCell2 = cell_pList2.length();
         if(numCell1*numCell2 < 2) continue;
         pMin = std::min(numCell1,numCell2);
         pMax = std::max(numCell1,numCell2);

         // copy the iterators to a vector in order to shuffle
         vector_part1_ptrs.clear();
         vector_part1_ptrs.reserve(numCell1);
         ListIterator<JustinsParticlePtr> lit1(cell_pList1);
         for (lit1.begin(); lit1.ok(); ++lit1) vector_part1_ptrs.push_back(lit1());
         std::shuffle(vector_part1_ptrs.begin(),vector_part1_ptrs.end(),global_rand_gen); 
         
         // do the same for species 2
         vector_part2_ptrs.clear();
         vector_part2_ptrs.reserve(numCell2);
         ListIterator<JustinsParticlePtr> lit2(cell_pList2);
         for (lit2.begin(); lit2.ok(); ++lit2) vector_part2_ptrs.push_back(lit2());
         std::shuffle(vector_part2_ptrs.begin(),vector_part2_ptrs.end(),global_rand_gen); 
	 
	 Real Wtot0 = 0.0;
	 Real Etot0 = 0.0;
	 if(m_enforce_conservations) { // save the initial mass and Energy
	    Real Wtot01 = 0.0, Etot01 = 0.0;
            for (lit1.begin(); lit1.ok(); ++lit1) {
               part1_ptr = lit1().getPointer();
               const std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real& wp1 = part1_ptr->weight();
	       Wtot01 += wp1;
	       for (int n=0; n<3; n++) Etot01 += wp1/2.0*betap1[n]*betap1[n];
	    }

	    Real Wtot02 = 0.0, Etot02 = 0.0;
            for (lit2.begin(); lit2.ok(); ++lit2) {
               part2_ptr = lit2().getPointer();
               const std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real& wp2 = part2_ptr->weight();
	       Wtot02 += wp2;
	       for (int n=0; n<3; n++) Etot02 += wp2/2.0*betap2[n]*betap2[n];
	    }
	    Wtot0 = m_mass1*Wtot01 + m_mass2*Wtot02;
	    Etot0 = m_mass1*Etot01 + m_mass2*Etot02;
         }

         unsigned int p1, p2;
	 std::array<Real,3> dBetaAvg = {0.0,0.0,0.0};
         for (auto p=0; p<pMax; p++) { // loop over particle scattering pairs
     
            if(pMin==numCell1) {
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
	    den12 = wpMax*pMin/cellV_SI;
	    
	    if(m_Nsubcycle_max>1) {
               wpMin = std::min(wp1,wp2);
      	       Nsubcycle = std::round(wpMax/wpMin);
      	       Nsubcycle = std::min(Nsubcycle,m_Nsubcycle_max);
               den12 = den12/Nsubcycle;
	    }

	    for (int m=0; m<Nsubcycle; m++) {

#ifdef RELATIVISTIC_PARTICLES
               if( (float)wp2 < (float)wp1 ) {
	          bool scatter1 = true;
                  if(MathUtils::rand()>wp2/wp1) scatter1 = false; 
                  LorentzScatter( betap2, betap1, scatter1, m_mass2, m_mass1, den12, bmax, a_dt_sec );
               }
               else {
                  bool scatter2 = true;
		  if( (float)wp1 < (float)wp2 ) {if(MathUtils::rand()>wp1/wp2) scatter2 = false;}
                  LorentzScatter( betap1, betap2, scatter2, m_mass1, m_mass2, den12, bmax, a_dt_sec );
               }
#else
               // compute deltaU
               std::array<Real,3> deltaU;
               GalileanScatter( deltaU, betap1, betap2, den12, bmax, a_dt_sec );     

               // update particle velocities
               if( (float)wp1 == (float)wp2 ) {
		  for (int n=0; n<3; n++) betap1[n] += m_mu/m_mass1*deltaU[n];
                  for (int n=0; n<3; n++) betap2[n] -= m_mu/m_mass2*deltaU[n];
	       }
	       else if( (float)wp1 < (float)wp2 ) {
		  for (int n=0; n<3; n++) betap1[n] += m_mu/m_mass1*deltaU[n];
	          for (int n=0; n<3; n++) dBetaAvg[n] += wp1*m_mu*deltaU[n];
		  if(MathUtils::rand()<wp1/wp2) {
                     for (int n=0; n<3; n++) betap2[n] -= m_mu/m_mass2*deltaU[n];
	             for (int n=0; n<3; n++) dBetaAvg[n] -= wp2*m_mu*deltaU[n];
		  }
	       }
               else { //if( (float)wp2 < (float)wp1 ) {
		  if(MathUtils::rand()<wp2/wp1) {
		     for (int n=0; n<3; n++) betap1[n] += m_mu/m_mass1*deltaU[n];
	             for (int n=0; n<3; n++) dBetaAvg[n] += wp1*m_mu*deltaU[n];
		  }
                  for (int n=0; n<3; n++) betap2[n] -= m_mu/m_mass2*deltaU[n];
	          for (int n=0; n<3; n++) dBetaAvg[n] -= wp2*m_mu*deltaU[n];
	       }

	    }
#endif
         }
         verbosity=0;

         // adjust particle betas such that momentum and energy are identically conserved
	 Real dBetaSq = 0.0;
	 for (int n=0; n<3; n++) dBetaSq += dBetaAvg[n]*dBetaAvg[n];
	 if(dBetaSq>0.0 && m_enforce_conservations) {
	
	    // enforce momentum conservation and compute deltaE
	    for (int n=0; n<3; n++) dBetaAvg[n] /= Wtot0;
	    long double Etot11 = 0.0;
	    for (lit1.begin(); lit1.ok(); ++lit1) {
               part1_ptr = lit1().getPointer();
               std::array<Real,3>& betap1 = part1_ptr->velocity();
               const Real& wp1 = part1_ptr->weight();
	       for (int n=0; n<3; n++) betap1[n] -= dBetaAvg[n];
	       for (int n=0; n<3; n++) Etot11 += wp1/2.0*betap1[n]*betap1[n];
	    }
	    Etot11 *= m_mass1;
	    
	    long double Etot12 = 0.0;
	    for (lit2.begin(); lit2.ok(); ++lit2) {
               part2_ptr = lit2().getPointer();
               std::array<Real,3>& betap2 = part2_ptr->velocity();
               const Real& wp2 = part2_ptr->weight();
	       for (int n=0; n<3; n++) betap2[n] -= dBetaAvg[n];
	       for (int n=0; n<3; n++) Etot12 += wp2/2.0*betap2[n]*betap2[n];
	    }
	    Etot12 *= m_mass2;
	    const long double Etot1 = Etot11 + Etot12;
	    const long double deltaE = Etot1 - Etot0;
	    long double deltaEp1, deltaEp2;
            if(numCell1==1) {
               deltaEp1 = 0.0;
	       deltaEp2 = deltaE;
	    }
	    else if(numCell2==1) { 
	       deltaEp1 = deltaE;
               deltaEp2 = 0.0;
	    }
	    else {
	       deltaEp1 = Etot11/Etot1*deltaE;
	       deltaEp2 = Etot12/Etot1*deltaE;
	    }

            // sort particles by weight to bias adjusting particles below to
	    // heavier weight ones and to get pairs that have similar weights
	    if(m_sort_weighted_particles) {
	       std::sort( vector_part1_ptrs.begin(), vector_part1_ptrs.end(),
	                  []( JustinsParticlePtr& lhs, JustinsParticlePtr& rhs )
	     	          { return lhs.getPointer()->weight() > rhs.getPointer()->weight(); } );
	       std::sort( vector_part2_ptrs.begin(), vector_part2_ptrs.end(),
	                  []( JustinsParticlePtr& lhs, JustinsParticlePtr& rhs )
	     	          { return lhs.getPointer()->weight() > rhs.getPointer()->weight(); } );
	    }
	    
	    // adjust species 1 particles to absorb deltaEp1
	    int loop1_count = 0;
            for (int p=0; p<vector_part1_ptrs.size(); p++) {
	       if(deltaEp1==0.0) break;
               int p1 = p;
   	       p++;
	       if(p==vector_part1_ptrs.size()) {
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
               modEnergyPairwise(betap1, betap2, m_mass1*wp1, m_mass1*wp2, deltaEp1);
	       if(p==vector_part1_ptrs.size()-1) {
		  loop1_count++;
		  p = -1;
	       }
	       if(loop1_count>m_loop_count_max) {
		  cout << "Notice: loop1_count > loop_count_max = " << m_loop_count_max << endl;
		  cout << "        for species1 = " << m_species1_name << endl;
		  cout << "        and species2 = " << m_species2_name << endl;
		  cout << "        num_parts1   = " << vector_part1_ptrs.size() << endl;
		  cout << "        ig           = " << ig << endl;
		  cout << "        Etot0        = " << Etot0 << endl;
		  cout << "        Etot1        = " << Etot1 << endl;
		  cout << "        Etot11       = " << Etot11 << endl;
		  cout << "        deltaEp1     = " << deltaEp1 << endl;
		  break;
	       }
	    }
	    
	    // adjust species 2 particles to absorb deltaEp1
	    int loop2_count = 0;
            for (int p=0; p<vector_part2_ptrs.size(); p++) {
	       if(deltaEp2==0.0) break;
               int p1 = p;
   	       p++;
	       if(p==vector_part2_ptrs.size()) {
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
               modEnergyPairwise(betap1, betap2, m_mass2*wp1, m_mass2*wp2, deltaEp2);
	       if(p==vector_part2_ptrs.size()-1) {
		  loop2_count++;
		  p = -1;
	       }
	       if(loop2_count>m_loop_count_max) {
		  cout << "Notice: loop2_count > loop_count_max = " << m_loop_count_max << endl;
		  cout << "        for species1 = " << m_species1_name << endl;
		  cout << "        and species2 = " << m_species2_name << endl;
		  cout << "        num_parts2   = " << vector_part2_ptrs.size() << endl;
		  cout << "        ig           = " << ig << endl;
		  cout << "        Etot0        = " << Etot0 << endl;
		  cout << "        Etot1        = " << Etot1 << endl;
		  cout << "        Etot12       = " << Etot12 << endl;
		  cout << "        deltaEp2     = " << deltaEp2 << endl;
		  break;
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
   int numCell1, numCell2, pMax, pMin, Nsubcycle=1;
   Real numDen1, wpMax, wpMin;
   Real numDen2, den12, bmax, cellV_SI;
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
         if(numDen1*numDen2 == 0.0) continue;
	 
	 // set maximum impact parameter
	 Real LDe = this_LDe.get(ig,0);
	 Real MaxN = std::max(numDen1,numDen2);
	 Real atomic_spacing = std::pow(4.189*MaxN,-1.0/3.0); // atomic spacing [m]
	 bmax = std::max(LDe,atomic_spacing);

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
         List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
         numCell1 = cell_pList1.length();
         numCell2 = cell_pList2.length();
         if(numCell1*numCell2 < 2) continue;
         pMin = std::min(numCell1,numCell2);
         pMax = std::max(numCell1,numCell2);

         // copy the iterators to a vector in order to shuffle
         vector_part1_ptrs.clear();
         vector_part1_ptrs.reserve(numCell1);
         ListIterator<JustinsParticlePtr> lit1(cell_pList1);
         for (lit1.begin(); lit1.ok(); ++lit1) vector_part1_ptrs.push_back(lit1());
         std::shuffle(vector_part1_ptrs.begin(),vector_part1_ptrs.end(),global_rand_gen); 
         
         // do the same for species 2
         vector_part2_ptrs.clear();
         vector_part2_ptrs.reserve(numCell2);
         ListIterator<JustinsParticlePtr> lit2(cell_pList2);
         for (lit2.begin(); lit2.ok(); ++lit2) vector_part2_ptrs.push_back(lit2());
         std::shuffle(vector_part2_ptrs.begin(),vector_part2_ptrs.end(),global_rand_gen); 
         
         unsigned int p1, p2;
         for (auto p=0; p<pMax; p++) { // loop over particle scattering pairs
     
            if(pMin==numCell1) {
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
	    den12 = wpMax*pMin/cellV_SI;
	    
	    if(m_Nsubcycle_max>1) {
               wpMin = std::min(wp1,wp2);
      	       Nsubcycle = std::floor(wpMax/wpMin);
      	       Nsubcycle = std::min(Nsubcycle,m_Nsubcycle_max);
               den12 = den12/Nsubcycle;
	    }

	    for (int m=0; m<Nsubcycle; m++) {

               // compute deltaU
               std::array<Real,3> deltaU;
               GalileanScatter( deltaU, betap1, betap2, den12, bmax, a_dt_sec );     

               // update particle velocities
	       if( (float)wp1 < (float)wp2 ) {

	          const long double ratio = wp1/wp2;
		  if(m==0) { // save betap and energy for particle 2 before scatter update
                     betap_before = betap2;               
		     Ebefore = m_mass2*(betap2[0]*betap2[0] + betap2[1]*betap2[1] + betap2[2]*betap2[2])/2.0;
		     betap_sum = {0.0,0.0,0.0};
		     Escatter_sum = 0.0;
		  }

		  // scatter particles
		  std::array<Real,3> betap2p = betap2;
                  for (int n=0; n<3; n++) betap1[n]  += m_mu/m_mass1*deltaU[n];
		  for (int n=0; n<3; n++) betap2p[n] -= m_mu/m_mass2*deltaU[n];

		  // update betap_sum and Escatter_sum
		  for (int n=0; n<3; n++) betap_sum[n] += betap2p[n];
		  Escatter_sum += m_mass2*(betap2p[0]*betap2p[0] + betap2p[1]*betap2p[1] + betap2p[2]*betap2p[2])/2.0;
		     
		  if(m==Nsubcycle-1) {
		     const Real Eafter = Ebefore + ratio*(Escatter_sum-Nsubcycle*Ebefore);
		     for (int n=0; n<3; n++) betap2[n] = betap_before[n] + ratio*(betap_sum[n]-Nsubcycle*betap_before[n]);

		     // correct beta of particle 2 to conserve energy exactly and momentum on average
		     enforceEnergyConservation( betap2, m_mass2, Eafter );
		  }
                    
	       }
	       else if( (float)wp2 < (float)wp1 ) {

	          const long double ratio = wp2/wp1;
	          if(m==0) { // save betap and energy for particle 1 before scatter update
	             betap_before = betap1;
		     Ebefore = m_mass1*(betap1[0]*betap1[0] + betap1[1]*betap1[1] + betap1[2]*betap1[2])/2.0;
		     betap_sum = {0.0,0.0,0.0};
		     Escatter_sum = 0.0;
                  }

		  // scatter particles
		  std::array<Real,3> betap1p = betap1;
                  for (int n=0; n<3; n++) betap1p[n] += m_mu/m_mass1*deltaU[n];
	          for (int n=0; n<3; n++) betap2[n]  -= m_mu/m_mass2*deltaU[n];
		     
		  // update betap_sum and Escatter_sum
		  for (int n=0; n<3; n++) betap_sum[n] += betap1p[n];
		  Escatter_sum += m_mass1*(betap1p[0]*betap1p[0] + betap1p[1]*betap1p[1] + betap1p[2]*betap1p[2])/2.0;
		     
		  if(m==Nsubcycle-1) {
		     const Real Eafter = Ebefore + ratio*(Escatter_sum-Nsubcycle*Ebefore);
		     for (int n=0; n<3; n++) betap1[n] = betap_before[n] + ratio*(betap_sum[n]-Nsubcycle*betap_before[n]);

   		     // correct beta of particle 1 to conserve energy exactly and momentum on average
		     enforceEnergyConservation( betap1, m_mass1, Eafter );
		  }

	       }
	       else {
	          for (int n=0; n<3; n++) betap1[n] += m_mu/m_mass1*deltaU[n];
                  for (int n=0; n<3; n++) betap2[n] -= m_mu/m_mass2*deltaU[n];
	       }

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
                         const Real                 a_dt_sec ) const
{
   CH_TIME("Coulomb::GalileanScatter()");
   
   Real ux = a_vp1[0]-a_vp2[0];
   Real uy = a_vp1[1]-a_vp2[1];
   Real uz = a_vp1[2]-a_vp2[2];
   Real u = sqrt(ux*ux + uy*uy + uz*uz);

   // compute s12 = 2.0*deltasq_var
   Real b90 = m_b90_fact/(m_mu*u*u); // b90 = q1*q2/(4*pi*ep0*WR) = q1*q2/(2*pi*ep0*mu*u^2)

   Real Clog = m_Clog;
   if(Clog==0.0 && u>0.0) {
      Real bmin_qm = m_bqm_fact/(m_mu*u);
      Real bmin = std::max(b90/2.0,bmin_qm);
      Clog = 0.5*std::log(1.0 + a_bmax*a_bmax/bmin/bmin);
      Clog = std::max(2.0,Clog);
   }
   m_s12 = Constants::PI*b90*b90*a_den12*Clog*u*Constants::CVAC*a_dt_sec;

   switch (m_angular_scattering) {
      case TAKIZUKA:

      if(m_s12<2.0) { // sample from gaussian distribution  
         const Real delta = sqrt(m_s12/2.0)*MathUtils::randn();
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
                        const Real                 a_dt_sec ) const
{
   CH_TIME("Coulomb::LorentzScatter()");
   
   // Relativistic scattering method for two particles with weight2 < weight1

   long double gamma1, gamma2, Etot, gammacm, gamma1st, gamma2st;
   long double vcmdotup, vrelst, muRst, upst_fact, upstsq, denom;
   std::array<Real,3> vcm, upst;

   // compute the lab frame total energy
   gamma1 = sqrt(1.0 + a_up1[0]*a_up1[0] + a_up1[1]*a_up1[1] + a_up1[2]*a_up1[2]);
   gamma2 = sqrt(1.0 + a_up2[0]*a_up2[0] + a_up2[1]*a_up2[1] + a_up2[2]*a_up2[2]);
   Etot = gamma1*a_mass1 + gamma2*a_mass2;
 
   // compute center-of-momentum beta and gamma
   for(int n=0; n<3; n++) vcm[n] = (a_mass1*a_up1[n] + a_mass2*a_up2[n])/Etot;
   gammacm = 1.0/sqrt(1.0-vcm[0]*vcm[0]-vcm[1]*vcm[1]-vcm[2]*vcm[2]);

   // compute gamma* and up* for particles 1 and 2 ( * = CM frame )
   vcmdotup = vcm[0]*a_up2[0] + vcm[1]*a_up2[1] + vcm[2]*a_up2[2];
   gamma2st = gammacm*(gamma2 - vcmdotup);

   vcmdotup = vcm[0]*a_up1[0] + vcm[1]*a_up1[1] + vcm[2]*a_up1[2];
   gamma1st = gammacm*(gamma1 - vcmdotup);
 
   upst_fact = (gammacm/(1.0+gammacm)*vcmdotup - gamma1)*gammacm;
   for(int n=0; n<3; n++) upst[n] = a_up1[n] + upst_fact*vcm[n]; 
 
   // compute relative beta in CM frame (vrelst = |v1st - v2st|/(1 - v1st/cdotv2st))
   muRst = gamma1st*a_mass1*gamma2st*a_mass2/(a_mass2*gamma2st + a_mass1*gamma1st); 
   upstsq = upst[0]*upst[0] + upst[1]*upst[1] + upst[2]*upst[2];
   denom = 1.0 + upstsq*a_mass1/a_mass2/gamma1st/gamma2st;
   vrelst = sqrt(upstsq)*a_mass1/muRst/denom;
   
   Real Clog = m_Clog;
   if(Clog==0.0 && upstsq>0.0) {
      Real b90 = m_b90_fact/(a_mass1*gamma1st + a_mass2*gamma2st)
               *(gamma1st*gamma2st*a_mass2/a_mass1/upstsq + 1.0);
      Real bmin_qm = m_bqm_fact/(a_mass1*std::sqrt(upstsq));
      Real bmin = std::max(b90/2.0,bmin_qm);
      Clog = 0.5*std::log(1.0 + a_bmax*a_bmax/bmin/bmin);
      Clog = std::max(2.0,Clog);
   }

   // compute s12 = 2*<delta^2>
   m_s12 = Constants::PI*m_b90_fact*m_b90_fact*a_den12*Clog*vrelst*Constants::CVAC*a_dt_sec;
   m_s12 *= gamma1st*gamma2st/gamma1/gamma2;
   m_s12 /= pow(muRst*vrelst*vrelst,2);
   
   switch (m_angular_scattering) {
      case TAKIZUKA:

      if(m_s12<2.0) { // sample from gaussian distribution  
         const Real delta = sqrt(m_s12/2.0)*MathUtils::randn();
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
   }
 
   // set random azimuthal angle
   m_phi = Constants::TWOPI*MathUtils::rand();
   m_cosphi = cos(m_phi);
   m_sinphi = sin(m_phi);

   // rotate upst for particle 1 by scattering angles
   ScatteringUtils::rotateVelocity(upst,m_costh,m_sinth,m_cosphi,m_sinphi);

   // Lorentz transform particle 1 back to lab frame
   vcmdotup = vcm[0]*upst[0] + vcm[1]*upst[1] + vcm[2]*upst[2];
   upst_fact = (gammacm/(1.0+gammacm)*vcmdotup + gamma1st)*gammacm;
   for(int n=0; n<3; n++) a_up1[n] = upst[n] + upst_fact*vcm[n]; 
  
   if(a_scatter2) { 
     // Lorentz transform particle 2 back to lab frame (pp2st + pp1st = 0)
      for(int n=0; n<3; n++) upst[n] *= -a_mass1/a_mass2;
      vcmdotup *= -a_mass1/a_mass2;
      upst_fact = (gammacm/(1.0+gammacm)*vcmdotup + gamma2st)*gammacm;
      for(int n=0; n<3; n++) a_up2[n] = upst[n] + upst_fact*vcm[n];
      // set velocity of particle 2 using momentum conservation
      //for(int n=0; n<3; n++) a_up2[n] = (vcm[n]*Etot - a_mass1*a_up1[n])/a_mass2;
   }
 
}

#include "NamespaceFooter.H"

