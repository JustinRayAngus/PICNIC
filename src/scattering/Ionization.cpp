
#include "Ionization.H"
#include "MathUtils.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"

#include <iostream>
#include <fstream>

#include "NamespaceHeader.H"


void Ionization::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                             const DomainGrid&           a_mesh )
{
   CH_TIME("Ionization::initialize()");
   
   const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
   const int num_species = a_pic_species_intf.numSpecies();   

   // get pointer to species 1
   CH_assert(m_sp1<num_species);
   const PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];

   // get pointer to species 2
   CH_assert(m_sp2<num_species);
   const PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];
   
   // get pointer to electron species
   CH_assert(m_spE<num_species);
   const PicChargedSpeciesPtr this_speciesE = pic_species_ptr_vect[species_map[m_spE]];
   
   // get pointer to ion species
   CH_assert(m_spI<num_species);
   const PicChargedSpeciesPtr this_speciesI = pic_species_ptr_vect[species_map[m_spI]];

   // set the species names
   m_species1_name = this_species1->name();
   m_species2_name = this_species2->name();
   m_speciesE_name = this_speciesE->name();
   m_speciesI_name = this_speciesI->name();

   // set the species masses
   m_mass1 = this_species1->mass();          // species 1 mass / me
   m_mass2 = this_species2->mass();          // species 2 mass / me
   m_massE = this_speciesE->mass();          // species 2 mass / me
   m_massI = this_speciesI->mass();          // species 2 mass / me
   m_mu = m_mass1*m_mass2/(m_mass1 + m_mass2); // reduced mass

   // pre-define some constants
   m_mcSq = Constants::ME*Constants::CVAC*Constants::CVAC/Constants::QE; // 0.511e6 
         
   // load cross section from file
   if(m_xsec_type == TEXT_FILE) {
      readCrossSectionFile();
   }
   
   if (m_verbosity)  printParameters();

}

void Ionization::readCrossSectionFile()
{

   std::ifstream xsecFile;
   std::string this_line;
   bool data_start_line_found=false;
   int this_line_number=0;
   Real this_Uizn;
   Real this_energy;
   Real this_xsec;
      
   xsecFile.open(m_xsec_fname,std::ifstream::in);
   if(xsecFile.is_open()) {
      while(getline(xsecFile, this_line)) {
         stringstream ss;
         ss << this_line;
         if(this_line_number==0) CH_assert(this_line=="ionization" || this_line=="IONIZATION")
         if(this_line_number==1) m_reaction_label = this_line;
         if(this_line_number==2) {
            stringstream(this_line) >> this_Uizn;
            CH_assert(m_Uizn==this_Uizn);
         }
         if(data_start_line_found && this_line.compare(1,5,"-----")==0) {
            break;
         }
         if(data_start_line_found) {
            string temp;
            ss >> temp; // first element is energy [eV]
            stringstream(temp) >> this_energy;
            m_Evec.push_back(this_energy);
            ss >> temp; // second element is xsec [m^2]
            stringstream(temp) >> this_xsec;
            m_Qvec.push_back(this_xsec);
         }
         if(!data_start_line_found && this_line.compare(1,5,"-----")==0) {
            data_start_line_found = true;
         }
         ++this_line_number;
      }
      xsecFile.close();
   }
   else {
      cout << "Problem opening cross section file: " << m_xsec_fname << endl;
   }
   
}

void Ionization::setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const 
{
   CH_TIME("Ionization::setMeanFreeTime()");
   
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

void Ionization::setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                              const LevelData<FArrayBox>&  a_energyDensity1,
                              const LevelData<FArrayBox>&  a_numberDensity2,
                              const LevelData<FArrayBox>&  a_energyDensity2 ) const
{
   CH_TIME("Ionization::setInterMFT()");
  
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
         if(local_numberDensity1*local_numberDensity2 == 0.0) continue;

         local_betaSqEff1 = 0.0;
         local_betaSqEff2 = 0.0;
         for( int dir=0; dir<3; dir++) {
            local_betaSqEff1 += 2.0*this_energyDensity1.get(ig,dir);  
            local_betaSqEff2 += 2.0*this_energyDensity2.get(ig,dir);  
         }
         local_betaSqEff1 /= local_numberDensity1*m_mass1; // |V1/cvac|^2
         local_betaSqEff2 /= local_numberDensity2*m_mass2; // |V2/cvac|^2

         // compute local beta and cross section
         g12 = sqrt(local_betaSqEff1 + local_betaSqEff2);
         sigma = getSigma(g12);

         // compute local local nuMax
         local_nuMax = local_numberDensity2*sigma*g12*Constants::CVAC;

         // update box_nuMax for time step
         box_nuMax = Max(box_nuMax,local_nuMax);
      }
   
   }
   
   Real global_nuMax = box_nuMax;
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]

}

void Ionization::applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                            const DomainGrid&           a_mesh,
                            const Real                  a_dt_sec ) const
{
   CH_TIME("Ionization::applyScattering()");
   
   PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 

   PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];
   PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];
   if (!this_species1->scatter()) { return; }
   if (!this_species2->scatter()) { return; }

   PicChargedSpeciesPtr this_speciesE = pic_species_ptr_vect[species_map[m_spE]];
   PicChargedSpeciesPtr this_speciesI = pic_species_ptr_vect[species_map[m_spI]];

   // electron impact ionizaton e + A => e + e + A+
   electronImpact( *this_species1, *this_species2, 
                   *this_speciesE, *this_speciesI, a_mesh, a_dt_sec );

}
      
void Ionization::electronImpact( PicChargedSpecies&  a_picSpecies1,
                                 PicChargedSpecies&  a_picSpecies2, 
                                 PicChargedSpecies&  a_picSpeciesE, 
                                 PicChargedSpecies&  a_picSpeciesI, 
                           const DomainGrid&  a_mesh,
                           const Real         a_dt_sec ) const
{
   CH_TIME("Ionization::electronImpact()");
   
   // predefine some pointers to be used below
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  

   // define references to picSpecies1 (primary species)
   LevelData<BinFab<JustinsParticlePtr>>& data1_binfab_ptr = a_picSpecies1.partData_binfab();
   
   // define references to picSpecies2
   LevelData<BinFab<JustinsParticlePtr>>& data2_binfab_ptr = a_picSpecies2.partData_binfab();
   // get/set density here from binfab_ptr, which  doesn't contain particles created
   // during inelastic events, but does know about those that are killed ... 
   const LevelData<FArrayBox>& numberDensity2 = a_picSpecies2.getNumberDensityFromBinFab(true);
  
   // predefine some variables
   Real local_numberDensity2;
   int local_numCell1, local_numCell2;
   Real wpEI;
   Real g12, arg, q12;
   
   std::array<Real,3> deltaU;
   Real R, costh, sinth, phi, rand_num; 
 
   Real sigma;
   
   const Real cvacSq = Constants::CVAC*Constants::CVAC;
   const Real Uizn_norm = m_Uizn*Constants::QE/cvacSq/Constants::ME;

   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data1_binfab_ptr.disjointBoxLayout();
   DataIterator ditg(grids);

   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes
      
      // initialize lists needed to create new particles
      List<JustinsParticle> new_pList1;
      List<JustinsParticle> new_ele_pList;
      List<JustinsParticle> new_ion_pList;

      // get references to the scattering species
      BinFab<JustinsParticlePtr>& thisBinFab1_ptr = data1_binfab_ptr[ditg];
      BinFab<JustinsParticlePtr>& thisBinFab2_ptr = data2_binfab_ptr[ditg];
      const FArrayBox& this_numberDensity2 = numberDensity2[ditg];
   
      std::vector<JustinsParticlePtr> vector_part2_ptrs;
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
         const IntVect ig = gbit(); // grid index
         local_numberDensity2 = this_numberDensity2.get(ig,0);
       
         List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
         List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
         local_numCell1 = cell_pList1.length();
         local_numCell2 = cell_pList2.length();
         if(local_numCell1 < 1 || local_numCell2 < 1) continue;
    
         // copy species B iterators to a vector in order to randomly select
         vector_part2_ptrs.clear();
         vector_part2_ptrs.reserve(local_numCell2);
         ListIterator<JustinsParticlePtr> lit2(cell_pList2);
         for (lit2.begin(); lit2.ok(); ++lit2) vector_part2_ptrs.push_back(lit2());
        
         // loop over species A particles in this cell
         int random_index2;
         int numCollisions = 0;
         ListIterator<JustinsParticlePtr> lit1(cell_pList1);
         for (lit1.begin(); lit1.ok(); ++lit1) {

            // get particle data for first (primary) particle
            JustinsParticlePtr& this_part1 = lit1();
            this_part1_ptr = this_part1.getPointer();
            std::array<Real,3>& betap1 = this_part1_ptr->velocity();
            Real& wp1 = this_part1_ptr->weight();
            
            // get particle data for second particle
            random_index2 = MathUtils::randInt(0, local_numCell2-1);  
            JustinsParticlePtr& this_part2 = vector_part2_ptrs[random_index2];
            this_part2_ptr = this_part2.getPointer();
            const std::array<Real,3>& betap2 = this_part2_ptr->velocity();
            Real& wp2 = this_part2_ptr->weight();

            // compute relative beta magnitude (non-relativistic)
            g12 = 0.0;
            for (int dir=0; dir<3; dir++) g12 = g12 + pow(betap1[dir]-betap2[dir],2);
            g12 = sqrt(g12);

            // get the cross section
            sigma = getSigma(g12);
            if(sigma==0.0) continue;
         
            // determine if the pair collide and then do the collision
            arg = g12*Constants::CVAC*sigma*local_numberDensity2*a_dt_sec;
            q12 = 1.0 - exp(-arg);
            rand_num = MathUtils::rand();

            if(rand_num<=q12) { // this pair collides

               numCollisions = numCollisions + 1;

               // adjust weights to account for ionization
               if(wp1<=wp2) {
                  wpEI = wp1;
                  wp2 = wp2 - wp1;
               }
               else { // JRA: split primary ele and kill neutral in this case ...
                  wpEI = wp2;
                  Real wpne = wp1 - wp2; // new electron to create
                  wp1 = wp2;
                  wp2 = 0.0; // will be tagged to kill below
                  if(wpne>=1.0e-8*wpEI) {
                     const RealVect& Xpe = this_part1_ptr->position();
                     JustinsParticle particleNE(wpne, Xpe, betap1);
                     new_pList1.append(particleNE);
                  }
               }
               
               // compute new kinetic energy and decide how to partition between
               // the new and old electron

               //compute deltaU
               phi = Constants::TWOPI*MathUtils::rand();
               R = MathUtils::rand();
               costh = 1.0 - 2.0*R;
               sinth = sqrt(1.0 - costh*costh);
               ScatteringUtils::computeDeltaU(deltaU,betap1,betap2,Uizn_norm,m_mu,costh,sinth,phi);

               // update primary particle velocities
               for (int dir=0; dir<3; dir++) betap1[dir] = betap1[dir] + m_mu/m_mass1*deltaU[dir];
               
               // scatter the secondary electron and ion as in excitation
               std::array<Real,3> betapEI = betap2;
               for (int dir=0; dir<3; dir++) betapEI[dir] = betapEI[dir] - m_mu/m_mass2*deltaU[dir];
               
               // create the new ele/ion pair
               const RealVect& Xp = this_part2_ptr->position();
               JustinsParticle particleE(wpEI, Xp, betapEI);
               JustinsParticle particleI(wpEI, Xp, betapEI);
               
               // need to declare ID to make smoke tests not fail
               // better fix for future is to add option to remove IDs from outputs
               uint64_t ID = 0;
               particleE.setID(ID);
               particleI.setID(ID);
 
               // append new particle to the appropriate lists
               new_ele_pList.append(particleE);
               new_ion_pList.append(particleI);
               
               // update energy spent ionizing
               m_deltaE_izn += m_Uizn*wpEI*Constants::JOULE_PER_EV; // Joules
               
               // remove particles from list that are to be killed
               if(wp1<1.0e-8*wpEI) {
                  this_part1_ptr->setKillTag();
                  cell_pList1.remove(this_part1_ptr); // JRA, not sure this works !
                  local_numCell1 = local_numCell1 - 1;
               }
               if(wp2<1.0e-8*wpEI) {
                  vector_part2_ptrs.erase(vector_part2_ptrs.begin()+random_index2);
                  this_part2_ptr->setKillTag();
                  cell_pList2.remove(this_part2_ptr); // JRA, not sure this works !
                  local_numCell2 = local_numCell2 - 1;
               }
               
               // check that there are still particles left in the list      
               if(local_numCell1==0) break;
               if(local_numCell2==0) break;

               // lower the density by wpEI for next collision probability
               // local_numberDensity2 -= wpEI/dV_phys/Jacobian;
               //const LevelData<FArrayBox>& Jacobian = a_mesh.getJcc();  
               //const Real volume_scale = a_mesh.getVolumeScale();
               //const Real dV_mapped = a_mesh.getMappedCellVolume();
               //const Real dV_phys = dV_mapped*volume_scale;

            }

         } // end loop over collisions for this cell
         if(procID()==0 && verbosity) cout << "numCollisions = " << numCollisions << endl;
         verbosity=0;

      } // end loop over cells
      
      // add the new particles to the appropriate containers
      ParticleData<JustinsParticle>& pData1 = a_picSpecies1.partData();
      ParticleData<JustinsParticle>& pDataE = a_picSpeciesE.partData();
      ParticleData<JustinsParticle>& pDataI = a_picSpeciesI.partData();

      pData1[ditg].addItemsDestructive(new_pList1);
      pDataE[ditg].addItemsDestructive(new_ele_pList);
      pDataI[ditg].addItemsDestructive(new_ion_pList);
      CH_assert(new_pList1.length()==0);
      CH_assert(new_ele_pList.length()==0);
      CH_assert(new_ion_pList.length()==0);
      
      // remove particles from species 1 if tagged to kill
      List<JustinsParticle>& pList1 = pData1[ditg].listItems();
      ListIterator<JustinsParticle> lit1(pList1);
      for(lit1.begin(); lit1.ok();) {
         const int& kill = lit1().killTag();
         if(kill) pList1.remove(lit1);
         else ++lit1;
      }
      
      // remove particles from species 2 if tagged to kill
      ParticleData<JustinsParticle>& pData2 = a_picSpecies2.partData();
      List<JustinsParticle>& pList2 = pData2[ditg].listItems();
      ListIterator<JustinsParticle> lit(pList2);
      for(lit.begin(); lit.ok();) {
         const int& kill = lit().killTag();
         if(kill) pList2.remove(lit);
         else ++lit;
      }

   } // end loop over boxes
   
   this_part1_ptr = NULL;
   this_part2_ptr = NULL;
   delete this_part1_ptr;
   delete this_part2_ptr;
   
}

Real Ionization::getSigma( const Real  a_beta ) const
{
   Real sigma = 0.0;
   if(m_xsec_type==LOTZ) sigma = getLotzSigma(a_beta);
   if(m_xsec_type==TEXT_FILE) sigma = getTextSigma(a_beta);
   return sigma;
}

Real Ionization::getLotzSigma( const Real  a_beta ) const
{
   CH_TIME("Ionization::getLotzSigma()");
 
   const Real betaSq = a_beta*a_beta;
   Real gamma;
   if(a_beta<1.0e-3) gamma = 1 + 0.5*betaSq;
   else gamma = 1.0/sqrt(1.0-betaSq);

   const Real KE = m_mu*m_mcSq*(gamma - 1.0);
   const Real Eclassic = m_mu*m_mcSq*betaSq/2.0;

   Real sigma = 0.0;
   if(KE>m_Uizn) {
     Real logTerm = (log(Eclassic/m_Uizn*gamma*gamma)-betaSq)/(Eclassic*m_Uizn);
     //Real logTerm = log(KE/m_Uizn)/(KE*m_Uizn);
     sigma = m_lotz_a*m_lotz_xi*logTerm*(1.0-m_lotz_b*exp(-m_lotz_c*(KE/m_Uizn-1.0))); // [m^2]
   }

   return sigma;   

}  

Real Ionization::getTextSigma( const Real  a_beta ) const
{
   CH_TIME("Ionization::getTextSigma()");
 
   const Real betaSq = a_beta*a_beta;
   const Real KE_cm = m_mu*m_mcSq*betaSq/2.0; // relative energy in cm frame [eV]

   Real sigma;
   const int N = m_Evec.size();
   int index = 0;

   if(KE_cm >= m_Evec.back()) { // ~ log(E)/E
      //sigma = m_Qvec.back();
      sigma = m_Qvec.back()*log(KE_cm)/log(m_Evec.back())*m_Evec.back()/KE_cm;
   }
   else if(KE_cm > m_Uizn) {
      index = N/2;
      while (KE_cm < m_Evec[index]) index--;
      while (KE_cm > m_Evec[index+1]) index++;
      if(m_use_loglog_interp && m_Qvec[index]>0.0) {
         sigma = ScatteringUtils::loglogInterp(m_Evec,m_Qvec,KE_cm,index);
      }
      else {
         sigma = MathUtils::linearInterp(m_Evec,m_Qvec,KE_cm,index);
      }
   }
   else {
      sigma = 0.0;
   }

   return sigma;   

}  

#include "NamespaceFooter.H"

