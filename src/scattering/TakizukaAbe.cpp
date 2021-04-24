
#include "TakizukaAbe.H"
#include "Constants.H"
#include "MathUtils.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ParticleUtils.H"

#include "NamespaceHeader.H"


void TakizukaAbe::initialize( const DomainGrid&         a_mesh,
                              const PicSpeciesPtrVect&  a_picSpeciesPtrVect )
{
   CH_TIME("TakizukaAbe::initialize()");
   
   // get pointer to species 1 and assert collisions allowed
   CH_assert(m_sp1<a_picSpeciesPtrVect.size());
   PicSpeciesPtr this_picSpecies1(a_picSpeciesPtrVect[m_sp1]);
   CH_assert(this_picSpecies1->scatter());
      
   // get pointer to species 2 and assert collisions allowed
   CH_assert(m_sp2<a_picSpeciesPtrVect.size());
   PicSpeciesPtr this_picSpecies2(a_picSpeciesPtrVect[m_sp2]);
   CH_assert(this_picSpecies2->scatter());

   // set the species names
   m_species1_name = this_picSpecies1->name();
   m_species2_name = this_picSpecies2->name();
   
   // set the species charges and assert they are not neutrals
   m_charge1 = this_picSpecies1->charge();  // species 1 charge / |qe|
   m_charge2 = this_picSpecies2->charge();  // species 2 charge / |qe|
   CH_assert(m_charge1!=0);
   CH_assert(m_charge2!=0);

   // set the species masses
   m_mass1 = this_picSpecies1->mass();    // species 1 mass / me
   m_mass2 = this_picSpecies2->mass();    // species 2 mass / me
   m_mu = m_mass1*m_mass2/(m_mass1 + m_mass2); // reduced mass
         
   Real cvacSq = Constants::CVAC*Constants::CVAC;
   m_b90_cofactor = abs(m_charge1*m_charge2)/(m_mu*cvacSq)*m_b90_codeToPhys; 
   
   //
   // set the mean free time
   //
   
   // define references to picSpecies1
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity1 = this_picSpecies1->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity1 = this_picSpecies1->getEnergyDensity(setMoments);
   
   if(m_sp1==m_sp2) {
      setMeanFreeTime(numberDensity1,energyDensity1);
   }
   else {
      // define references to picSpecies2
      const LevelData<FArrayBox>& numberDensity2 = this_picSpecies2->getNumberDensity(setMoments);
      const LevelData<FArrayBox>& energyDensity2 = this_picSpecies2->getEnergyDensity(setMoments);
      setMeanFreeTime(numberDensity1,energyDensity1,numberDensity2,energyDensity2);
   }

   if (m_verbosity)  printParameters();

}

void TakizukaAbe::setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity,
                                   const LevelData<FArrayBox>&  a_energyDensity ) const
{
   CH_TIME("TakizukaAbe::setMeanFreeTime()");
 
   Real Teff_eV, numberDensity, energyDensity, Clog;
   Real tau;
   Real box_nuMax=0.0; // for scattering time step calculation
   
   Real cvacSq = Constants::CVAC*Constants::CVAC;
    
   const DisjointBoxLayout& grids = a_numberDensity.disjointBoxLayout();
   DataIterator ditg(grids);
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity = a_numberDensity[ditg];
      const FArrayBox& this_energyDensity = a_energyDensity[ditg];
     
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         const IntVect ig = gbit(); // grid index
       
         // get local density and temperature and compute local VTeff
         numberDensity = this_numberDensity.get(ig,0); // [1/m^3]
         if(numberDensity == 0.0) continue;
         
         energyDensity = 0.0;
         for( int dir=0; dir<3; dir++) {
            energyDensity = energyDensity + this_energyDensity.get(ig,dir);  
         }
         Teff_eV = Constants::ME*2.0/3.0*energyDensity/numberDensity*cvacSq; // [Joules]
         Teff_eV = Constants::EV_PER_JOULE*Teff_eV; // local temperature [eV]

         Clog = 10.0;
         
         // define self-species collision time
         tau = 3.44e5*pow(Teff_eV,1.5)/(numberDensity*Constants::M3_PER_CM3)/Clog; // [s]
         tau = tau*sqrt(m_mass1/2.0)/pow(m_charge1*m_charge2,2);

         // compute nuMax [Hz]
         box_nuMax = Max(box_nuMax,1.0/tau);
      }
   
   }
   
   Real global_nuMax = box_nuMax;
#ifdef CH_MPI
   MPI_Allreduce( &box_nuMax, &global_nuMax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
#endif
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]
   
}

void TakizukaAbe::setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity1,
                                   const LevelData<FArrayBox>&  a_energyDensity1,
                                   const LevelData<FArrayBox>&  a_numberDensity2,
                                   const LevelData<FArrayBox>&  a_energyDensity2 ) const
{
   CH_TIME("TakizukaAbe::setMeanFreeTime()");
 
   //std::array<Real,8> x0 = {0.1,0.5,1.0,3.0,6.0,9.0,10.0,100};
   //for (auto i=0; i<x0.size(); i++) {
   //   if(!procID()) cout << "JRA: x0 = " << x0[i] << endl; 
   //   if(!procID()) cout << "JRA: gammainc(x0,1.5) = " << MathUtils::gammainc(x0[i],1.5) << endl; 
   //}
  
   Real Teff1_eV, numberDensity1, energyDensity1, energy1, VT1;
   Real Teff2_eV, numberDensity2, energyDensity2, energy2, VT2;
   Real x12, psi12, nu012, nu12;
   Real x21, psi21, nu021, nu21;
         
   const Real factor = pow(Constants::QE*m_charge1*Constants::QE*m_charge2/Constants::EP0,2)/Constants::FOURPI; // [Joules-m]

   Real tau, Clog;
   Real box_nuMax=0.0; // for scattering time step calculation
   
   Real cvacSq = Constants::CVAC*Constants::CVAC;
 
   //bool verbosity = true;
   bool verbosity = false;

   const DisjointBoxLayout& grids = a_numberDensity1.disjointBoxLayout();
   DataIterator ditg(grids);
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity1 = a_numberDensity1[ditg];
      const FArrayBox& this_energyDensity1 = a_energyDensity1[ditg];
      const FArrayBox& this_numberDensity2 = a_numberDensity2[ditg];
      const FArrayBox& this_energyDensity2 = a_energyDensity2[ditg];
     
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         const IntVect ig = gbit(); // grid index
       
         // get local density and temperature and compute local VTeff
         numberDensity1 = this_numberDensity1.get(ig,0); // [1/m^3]
         numberDensity2 = this_numberDensity2.get(ig,0); // [1/m^3]
         if(numberDensity1*numberDensity2 == 0.0) continue;
         
         energyDensity1 = 0.0;
         energyDensity2 = 0.0;
         for( int dir=0; dir<3; dir++) {
            energyDensity1 = energyDensity1 + this_energyDensity1.get(ig,dir);  
            energyDensity2 = energyDensity2 + this_energyDensity2.get(ig,dir);  
         }
         energy1 = Constants::ME*energyDensity1/numberDensity1*cvacSq; // [Joules]
         energy2 = Constants::ME*energyDensity2/numberDensity2*cvacSq; // [Joules]
         Teff1_eV = Constants::EV_PER_JOULE*2.0/3.0*energy1; // local temperature [eV]
         Teff2_eV = Constants::EV_PER_JOULE*2.0/3.0*energy2; // local temperature [eV]
         VT1 = sqrt(Constants::QE*Teff1_eV/(Constants::ME*m_mass1)); // [m/s]
         VT2 = sqrt(Constants::QE*Teff2_eV/(Constants::ME*m_mass2)); // [m/s]

         Clog = 10.0;

         x12 = (Teff1_eV/m_mass1)/(Teff2_eV/m_mass2); // ratio of velocities squared
         x21 = 1./x12;

         psi12 = 2.0/sqrt(Constants::PI)*MathUtils::gammainc(x12,1.5);
         psi21 = 2.0/sqrt(Constants::PI)*MathUtils::gammainc(x21,1.5);

         nu012 = factor*Clog*numberDensity2/(pow(energy1,2))*VT1; // [Hz]
         nu021 = factor*Clog*numberDensity1/(pow(energy2,2))*VT2; // [Hz]
         
         nu12 = (1.0 + m_mass1/m_mass2)*psi12*nu012;
         nu21 = (1.0 + m_mass2/m_mass1)*psi21*nu021;

         if(!procID() && verbosity) {
            cout << "JRA: mass1 = " << m_mass1<< endl;
            cout << "JRA: mass2 = " << m_mass2<< endl;
            cout << "JRA: Teff1_eV = " << Teff1_eV << endl;
            cout << "JRA: Teff2_eV = " << Teff2_eV << endl;
            cout << "JRA: VT1 = " << VT1<< endl;
            cout << "JRA: VT2 = " << VT2<< endl;
            cout << "JRA: x12 = " << x12<< endl;
            cout << "JRA: x21 = " << x21<< endl;
            cout << "JRA: numberDen1 = " << numberDensity1 << endl;
            cout << "JRA: numberDen2 = " << numberDensity2 << endl;
            cout << "JRA: psi12 = " << psi12 << endl;
            cout << "JRA: psi21 = " << psi21 << endl;
            cout << "JRA: nu12 = " << nu12 << endl;
            cout << "JRA: nu21 = " << nu21 << endl;
         }

         // compute nuMax [Hz]
         box_nuMax = Max(box_nuMax,nu12);
         box_nuMax = Max(box_nuMax,nu21);

      }
   
   }
   
   Real global_nuMax = box_nuMax;
#ifdef CH_MPI
   MPI_Allreduce( &box_nuMax, &global_nuMax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
#endif
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]

}
      
void TakizukaAbe::applySelfScattering( PicSpecies&  a_picSpecies, 
                                 const DomainGrid&  a_mesh,
                                 const Real         a_dt_sec ) const
{
   CH_TIME("TakizukaAbe::applySelfScattering()");
 
   // predefine some pointers to be used below
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  
 
   // define references to a_picSpecies
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab_ptr = a_picSpecies.partData_binfab();
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity = a_picSpecies.getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity = a_picSpecies.getEnergyDensity(setMoments);

   // predefine some variables
   int numCell;
   Real Teff_eV, tau, numDen, eneDen;
   Real Clog=10.0;
   Real box_nuMax=0.0;
   std::array<Real,3> deltaU;
 
   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data_binfab_ptr.disjointBoxLayout();
   DataIterator ditg(grids);
   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity = numberDensity[ditg];
      const FArrayBox& this_energyDensity = energyDensity[ditg];
     
      BinFab<JustinsParticlePtr>& thisBinFab_ptr = data_binfab_ptr[ditg];
   
      std::vector<JustinsParticlePtr> vector_part_ptrs;
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
         const IntVect ig = gbit(); // grid index
       
         // get local density and temperature and compute local tau
         numDen = this_numberDensity.get(ig,0);
         if(numDen == 0.0) continue;
         /*
         eneDen = 0.0;
         for( int dir=0; dir<3; dir++) {
            eneDen = eneDen + this_energyDensity.get(ig,dir);  
         }
         Teff_eV = Constants::ME*2.0/3.0*eneDen/numDen; // [Joules]
         Teff_eV = Constants::EV_PER_JOULE*Teff_eV; // local temperature [eV]

         tau = 3.44e5*pow(Teff_eV,1.5)/(numDen*Constants::M3_PER_CM3)/Clog; // [s]
         if(m_charge1>0) tau = tau*sqrt(m_mass1/2.0);
 
         box_nuMax = Max(box_nuMax,1.0/tau);
         if(box_nuMax*a_dt_sec>10.0) { 
            if(procID()) {
               cout << "WARNING: box_nuMaxDt = " << box_nuMax*a_dt_sec << endl;
               cout << "WARNING: Teff_eV = " << Teff_eV << endl;
               cout << "WARNING: m_charge1 = " << m_charge1 << endl;
               cout << "WARNING: m_mass1 = " << m_mass1 << endl;
               cout << "WARNING: numDen= " << numDen << endl;
            }
         }
         */

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         numCell = cell_pList.length();
         if(numCell < 2) continue;
          
         // set the start index for the loop over particle collisions
         int pstart = 3;        
         if(numCell % 2 == 0) pstart = 0;
         if(procID()==0 && verbosity) {
            cout << "numCell = " << numCell << endl;
            cout << "pstart = " << pstart << endl;
         }
    
         // copy the iterators to a vector in order to shuffle
         vector_part_ptrs.clear();
         vector_part_ptrs.reserve(numCell);
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) vector_part_ptrs.push_back(lit());
         std::shuffle(vector_part_ptrs.begin(),vector_part_ptrs.end(),global_rand_gen); 
         
         for (auto p=pstart; p<vector_part_ptrs.size(); p++) { // loop over particle scattering pairs
 
            // get particle data for first particle    
            JustinsParticlePtr& this_part1 = vector_part_ptrs[p];
            this_part1_ptr = this_part1.getPointer();
            std::array<Real,3>& this_betap1 = this_part1_ptr->velocity();
            
            // get particle data for second particle    
            p++; // advance the loop index by 1
            JustinsParticlePtr& this_part2 = vector_part_ptrs[p];
            this_part2_ptr = this_part2.getPointer();
            std::array<Real,3>& this_betap2 = this_part2_ptr->velocity();
   
            // compute deltaU
            computeDeltaU( deltaU,
                           this_betap1, numDen,
                           this_betap2, numDen,
                           Clog, a_dt_sec );     
            //deltaU = {0,0,0};
            
            // update particle velocities
            for (int dir=0; dir<3; dir++) {
               this_betap1[dir] = this_betap1[dir] + m_mu/m_mass1*deltaU[dir];
               this_betap2[dir] = this_betap2[dir] - m_mu/m_mass2*deltaU[dir];
            }

         }
         if(pstart==3) {
            int p1, p2;
            for (int p=0; p<pstart; p++) { // scatter particles 0,1, and 2
               p1 = p % 2;
               if(p==0) p2 = 1;
               if(p==1) p2 = 2;
               if(p==2) p2 = 2;
  
               // get particle data for first particle    
               JustinsParticlePtr& this_part1 = vector_part_ptrs[p1];
               this_part1_ptr = this_part1.getPointer();
               std::array<Real,3>& this_betap1 = this_part1_ptr->velocity();
            
               // get particle data for second particle    
               JustinsParticlePtr& this_part2 = vector_part_ptrs[p2];
               this_part2_ptr = this_part2.getPointer();
               std::array<Real,3>& this_betap2 = this_part2_ptr->velocity();

               // compute deltaU
               computeDeltaU( deltaU,
                              this_betap1, numDen/2.0,
                              this_betap2, numDen/2.0,
                              Clog, a_dt_sec );     
               //deltaU = {0,0,0};
            
               // update particle velocities
               for (int dir=0; dir<3; dir++) {
                  this_betap1[dir] = this_betap1[dir] + m_mu/m_mass1*deltaU[dir];
                  this_betap2[dir] = this_betap2[dir] - m_mu/m_mass2*deltaU[dir];
               }

            }
         }
         verbosity=0;

      } // end loop over cells

   } // end loop over boxes
  
 
   // don't forget to set pointers back to NULL and delete
   this_part1_ptr = NULL;
   this_part2_ptr = NULL;
   delete this_part1_ptr;
   delete this_part2_ptr;
   
   // While we are here, update the mean free time
   setMeanFreeTime(numberDensity,energyDensity);

}

void TakizukaAbe::applyInterScattering( PicSpecies&  a_picSpecies1,
                                        PicSpecies&  a_picSpecies2, 
                                  const DomainGrid&  a_mesh,
                                  const Real         a_dt_sec ) const
{
   CH_TIME("TakizukaAbe::applyInterScattering()");
   
   CH_assert(m_species1_name==a_picSpecies1.name());
   CH_assert(m_species2_name==a_picSpecies2.name());

   // predefine some pointers to be used below
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  
   
   // define references to a_picSpecies1
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab1_ptr = a_picSpecies1.partData_binfab();
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity1 = a_picSpecies1.getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity1 = a_picSpecies1.getEnergyDensity(setMoments);
   
   // define references to a_picSpecies2
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab2_ptr = a_picSpecies2.partData_binfab();
   const LevelData<FArrayBox>& numberDensity2 = a_picSpecies2.getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity2 = a_picSpecies2.getEnergyDensity(setMoments);
  
 
   // predefine some variables
   int numCell1, numCell2, pMax, pMin;
   Real Teff1_eV, numDen1, eneDen1;
   Real Teff2_eV, numDen2, eneDen2;
   Real Clog=10.0;
   Real box_nuMax=0.0;
   std::array<Real,3> deltaU;
 
   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data_binfab1_ptr.disjointBoxLayout();
   DataIterator ditg(grids);
   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity1 = numberDensity1[ditg];
      const FArrayBox& this_energyDensity1 = energyDensity1[ditg];
      
      const FArrayBox& this_numberDensity2 = numberDensity2[ditg];
      const FArrayBox& this_energyDensity2 = energyDensity2[ditg];
     
      BinFab<JustinsParticlePtr>& thisBinFab1_ptr = data_binfab1_ptr[ditg];
      BinFab<JustinsParticlePtr>& thisBinFab2_ptr = data_binfab2_ptr[ditg];
   
      std::vector<JustinsParticlePtr> vector_part1_ptrs;
      std::vector<JustinsParticlePtr> vector_part2_ptrs;
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
         const IntVect ig = gbit(); // grid index
         
         // get local density and temperature and compute local tau
         numDen1 = this_numberDensity1.get(ig,0);
         numDen2 = this_numberDensity2.get(ig,0);
         if(numDen1*numDen2 == 0.0) continue;
         /*
         eneDen1 = 0.0;
         eneDen2 = 0.0;
         for( int dir=0; dir<3; dir++) {
            eneDen1 = eneDen1 + this_energyDensity1.get(ig,dir);  
            eneDen2 = eneDen2 + this_energyDensity2.get(ig,dir);  
         }
         Teff1_eV = Constants::ME*2.0/3.0*eneDen1/numDen1; // [Joules]
         Teff1_eV = Constants::EV_PER_JOULE*Teff1_eV; // local temperature [eV]
         Teff2_eV = Constants::ME*2.0/3.0*eneDen2/numDen2; // [Joules]
         Teff2_eV = Constants::EV_PER_JOULE*Teff2_eV; // local temperature [eV]
         */

         // get the number of particles in this cell
         List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
         List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
         numCell1 = cell_pList1.length();
         numCell2 = cell_pList2.length();
         if(numCell1*numCell2 < 2) continue;
         pMin = Min(numCell1,numCell2);
         pMax = Max(numCell1,numCell2);
          
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
            if(procID() && verbosity) {
               cout << "JRA: numCell1 = " << numCell1 << endl;
               cout << "JRA: numCell2 = " << numCell2 << endl;
               cout << "JRA: p = " << p << endl;
               cout << "JRA: p1 = " << p1 << endl;
               cout << "JRA: p2 = " << p2 << endl;
            }

            // get particle data for first particle    
            JustinsParticlePtr& this_part1 = vector_part1_ptrs[p1];
            this_part1_ptr = this_part1.getPointer();
            std::array<Real,3>& this_betap1 = this_part1_ptr->velocity();
            
            // get particle data for second particle    
            JustinsParticlePtr& this_part2 = vector_part2_ptrs[p2];
            this_part2_ptr = this_part2.getPointer();
            std::array<Real,3>& this_betap2 = this_part2_ptr->velocity();
   
            // compute deltaU
            computeDeltaU( deltaU,
                           this_betap1, numDen1,
                           this_betap2, numDen2,
                           Clog, a_dt_sec );     
            //deltaU = {0,0,0};
            
            // update particle velocities
            for (int dir=0; dir<3; dir++) {
               this_betap1[dir] = this_betap1[dir] + m_mu/m_mass1*deltaU[dir];
               this_betap2[dir] = this_betap2[dir] - m_mu/m_mass2*deltaU[dir];
            }

         }
         verbosity=0;

      } // end loop over cells

   } // end loop over boxes


   // don't forget to set pointers back to NULL and delete
   this_part1_ptr = NULL;
   this_part2_ptr = NULL;
   delete this_part1_ptr;
   delete this_part2_ptr;
   
   // While we are here, update the mean free time
   setMeanFreeTime(numberDensity1,energyDensity1,numberDensity2,energyDensity2);

}

void TakizukaAbe::computeDeltaU( std::array<Real,3>&  a_deltaU,
                           const std::array<Real,3>&  a_vp1,
                           const Real                 a_den1,
                           const std::array<Real,3>&  a_vp2,
                           const Real                 a_den2,
                           const Real                 a_Clog,
                           const Real                 a_dt_sec ) const
{
   CH_TIME("TakizukaAbe::computeDeltaU()");
   
   Real ux = a_vp1[0]-a_vp2[0];
   Real uy = a_vp1[1]-a_vp2[1];
   Real uz = a_vp1[2]-a_vp2[2];
   Real u = sqrt(ux*ux + uy*uy + uz*uz);

   //Real cvacSq = Constants::CVAC*Constants::CVAC;

   // compute deltasq_var
   Real den = Min(a_den1,a_den2); 
   //Real b90 = abs(m_charge1*m_charge2)/(m_mu*u*u*cvacSq)*m_b90_codeToPhys; // 90 degree impact param [m] 
   Real b90 = m_b90_cofactor/(u*u); // 90 degree impact parameter [m]
   Real deltasq_var = Constants::TWOPI*b90*b90*den*a_Clog*u*Constants::CVAC*a_dt_sec;

   // sample from gaussian distribution  
   m_delta = sqrt(deltasq_var)*MathUtils::randn();
   m_deltasq = m_delta*m_delta;

   // set the polar scattering angle
   m_sinth = 2.0*m_delta/(1.0+m_deltasq); 
   m_costh = 1.0 - 2.0*m_deltasq/(1.0+m_deltasq);

   // set random azimuthal angle
   m_phi = Constants::TWOPI*MathUtils::rand();
   m_cosphi = cos(m_phi);
   m_sinphi = sin(m_phi);

   // define deltaU
   ParticleUtils::computeDeltaU(a_deltaU,ux,uy,uz,m_costh,m_sinth,m_cosphi,m_sinphi);
               
}

#include "NamespaceFooter.H"

