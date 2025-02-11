
#include "HardSphere.H"
#include "MathUtils.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"

#include "NamespaceHeader.H"


void HardSphere::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                             const DomainGrid&           a_mesh )
{
   CH_TIME("HardSphere::initialize()");
   
   const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
   const int num_species = a_pic_species_intf.numSpecies();   

   // get pointer to species 1
   CH_assert(m_sp1<num_species);
   const PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];
   CH_assert(this_species1->charge()==0);      

   // get pointer to species 2
   CH_assert(m_sp2<num_species);
   const PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];
   CH_assert(this_species2->charge()==0);      

   // set the species names
   m_species1_name = this_species1->name();
   m_species2_name = this_species2->name();

   // set the species masses
   m_mass1 = this_species1->mass();          // species 1 mass / me
   m_mass2 = this_species2->mass();          // species 2 mass / me
   m_mu = m_mass1*m_mass2/(m_mass1 + m_mass2); // reduced mass

   // define the radius of each species if not specified in input file
   if(m_r1==-1) {
      const Real A1 = m_mass1*Constants::ME/Constants::AMU;
      m_r1 = pow(A1,1./3.)*Constants::ABOHR;
   }
   if(m_r2==-1) {
      const Real A2 = m_mass2*Constants::ME/Constants::AMU;
      m_r2 = pow(A2,1./3.)*Constants::ABOHR;
   }
 
   // define total cross section
   m_sigmaT = Constants::PI*(m_r1 + m_r2)*(m_r1 + m_r2);
   
   if (m_verbosity)  printParameters();

}

void HardSphere::setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const 
{
   CH_TIME("HardSphere::setMeanFreeTime()");

   const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 

   const PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];
   const PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];

   if (!this_species1->scatter() || !this_species2->scatter()) { return; }

   const bool setMoments = false;
   const LevelData<FArrayBox>& numberDensity1 = this_species1->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity1 = this_species1->getEnergyDensity(setMoments);

   if (m_sp1==m_sp2) { setIntraMFT(numberDensity1,energyDensity1); }
   else {
      const LevelData<FArrayBox>& numberDensity2 = this_species2->getNumberDensity(setMoments);
      const LevelData<FArrayBox>& energyDensity2 = this_species2->getEnergyDensity(setMoments);
      setInterMFT(numberDensity1,energyDensity1,numberDensity2,energyDensity2);
   }

}

void HardSphere::setIntraMFT( const LevelData<FArrayBox>&  a_numberDensity,
                              const LevelData<FArrayBox>&  a_energyDensity ) const
{
   CH_TIME("HardSphere::setIntraMFT()");
   
   Real mass = m_mass1;
   Real cvacSq = Constants::CVAC*Constants::CVAC;
   
   // predefine some variables
   Real local_Teff, local_numberDensity, local_energyDensity, local_VTeff;
   Real local_nuMax;
   Real box_nuMax=0.0; // for scattering time step calculation
 
   // loop over lists in each cell and test shuffle
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
         local_numberDensity = this_numberDensity.get(ig,0);
         if(local_numberDensity == 0.0) continue;
         local_energyDensity = 0.0;
         for( int dir=0; dir<3; dir++) {
            local_energyDensity = local_energyDensity + this_energyDensity.get(ig,dir);  
         }
         local_Teff = 2.0/3.0*local_energyDensity/local_numberDensity*cvacSq; // M/me*(V[m/s])^2
         local_VTeff = sqrt(local_Teff/mass); // thermal speed [m/s]

         // compute local local nuMax*dt
         local_nuMax = local_numberDensity*m_sigmaT*local_VTeff;
         box_nuMax = Max(box_nuMax,local_nuMax);
      }
   
   }
   
   Real global_nuMax = box_nuMax;
//#ifdef CH_MPI
//   MPI_Allreduce( &box_nuMax, &global_nuMax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
//#endif
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]

}

void HardSphere::setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                              const LevelData<FArrayBox>&  a_energyDensity1,
                              const LevelData<FArrayBox>&  a_numberDensity2,
                              const LevelData<FArrayBox>&  a_energyDensity2 ) const
{
   CH_TIME("HardSphere::setInterMFT()");
   
   // predefine some variables
   Real local_Teff1, local_numberDensity1, local_energyDensity1, local_VTeff1;
   Real local_Teff2, local_numberDensity2, local_energyDensity2, local_VTeff2;
   Real local_nuMax1, local_nuMax2;
   Real box_nuMax=0.0;
   
   Real cvacSq = Constants::CVAC*Constants::CVAC;
 
   // loop over lists in each cell and test shuffle
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

         local_energyDensity1 = 0.0;
         local_energyDensity2 = 0.0;
         for( int dir=0; dir<3; dir++) {
            local_energyDensity1 = local_energyDensity1 + this_energyDensity1.get(ig,dir);  
            local_energyDensity2 = local_energyDensity2 + this_energyDensity2.get(ig,dir);  
         }
         local_Teff1 = 2.0/3.0*local_energyDensity1/local_numberDensity1*cvacSq; // M1/me*(V1[m/s])^2
         local_Teff2 = 2.0/3.0*local_energyDensity2/local_numberDensity2*cvacSq; // M2/me*(V2[m/s])^2
         local_VTeff1 = sqrt(local_Teff1/m_mass1); // thermal speed [m/s]
         local_VTeff2 = sqrt(local_Teff2/m_mass2); // thermal speed [m/s]

         // compute local local nuMax
         local_nuMax1 = local_numberDensity2*m_sigmaT*local_VTeff1;
         local_nuMax2 = local_numberDensity1*m_sigmaT*local_VTeff2;

         // update box_nuMax for time step
         box_nuMax = Max(box_nuMax,local_nuMax1);
         box_nuMax = Max(box_nuMax,local_nuMax2);
      }
   
   }
   
   Real global_nuMax = box_nuMax;
//#ifdef CH_MPI
//   MPI_Allreduce( &box_nuMax, &global_nuMax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
//#endif
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]

}

void HardSphere::applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                            const DomainGrid&           a_mesh,
                            const Real                  a_dt_sec ) const
{
   CH_TIME("HardSphere::applyScattering()");
   
   PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 

   PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];
   if (!this_species1->scatter()) { return; }

   if (m_sp1==m_sp2) {
      applySelfScattering( *this_species1, a_mesh, a_dt_sec );
   }
   else {
      PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];
      if (!this_species2->scatter()) { return; }
      applyInterScattering( *this_species1, *this_species2, a_mesh, a_dt_sec );
   }

}
      
      
void HardSphere::applySelfScattering( PicChargedSpecies&  a_picSpecies, 
                                const DomainGrid&  a_mesh,
                                const Real         a_dt_sec ) const
{
   CH_TIME("HardSphere::applySelfScattering()");
 
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  
   JustinsParticle* this_part3_ptr = NULL;  
 
   Real mass = m_mass1;
   Real cvacSq = Constants::CVAC*Constants::CVAC;
 
   // define reference to a_picSpcies binfab container of pointers to particle data
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab_ptr = a_picSpecies.partData_binfab();
   const LevelData<FArrayBox>& numberDensity = a_picSpecies.getNumberDensity(false);
   const LevelData<FArrayBox>& energyDensity = a_picSpecies.getEnergyDensity(true);

   // predefine some variables
   Real local_Teff, local_numberDensity, local_energyDensity, local_gmax;
   Real local_nuMaxDt, local_Nmax, local_Nmax_remainder, local_Nmax_whole;
   int local_numCell;
   Real g12, q12;
   
   std::array<Real,3> deltaU;
   Real R, costh, sinth, phi; 
 
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
       
         /////////////////////////////////////////////////////////////////
         //
         // compute the number of possible collisions for this cell
         //
          
         // get local density and temperature and compute local gmax
         local_numberDensity = this_numberDensity.get(ig,0);
         if (local_numberDensity == 0.0) { continue; }
         local_energyDensity = 0.0;
         for (int dir=0; dir<3; dir++) {
            local_energyDensity = local_energyDensity + this_energyDensity.get(ig,dir);  
         }
         local_Teff = 2.0/3.0*local_energyDensity/local_numberDensity*cvacSq; // M/me*(V[m/s])^2
         local_gmax = 5.0*sqrt(local_Teff/mass); // 5x thermal speed [m/s]

         // compute local nuMax*dt
         local_nuMaxDt = local_numberDensity*m_sigmaT*local_gmax*a_dt_sec;         
         if (local_nuMaxDt>10.0 && !m_ignore_warnings) {
            cout << "WARNING: local_nuMaxDt = " << local_nuMaxDt << endl;
         }

         // compute local maximum number of collsions: Nmax = 1/2*(N-1)*n*sigmaTmax*gmax*dt
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         local_numCell = cell_pList.length();
         if (local_numCell < 2) { continue; }
         local_Nmax = 0.5*(local_numCell-1)*std::min(local_nuMaxDt,1.0);
    
         // separate Nmax into whole integer and remainder and set integer Nmax
         local_Nmax_remainder = modf(local_Nmax, &local_Nmax_whole);
         Real rand_num = MathUtils::rand();
         int local_Nmax_integer = static_cast<int>(local_Nmax_whole);
         if (rand_num<=local_Nmax_remainder) {
            local_Nmax_integer = local_Nmax_integer + 1;
         }
         
         // copy the iterators to a vector in order to randomly select
         //std::vector<JustinsParticlePtr> vector_part_ptrs;
         vector_part_ptrs.clear();
         vector_part_ptrs.reserve(local_numCell);
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) vector_part_ptrs.push_back(lit());
        
         ////////////////////////////////////////////////////////////////////////
         //
         // loop over Nmax and randomly choose two unique elements from the vector to scatter
         //
         int random_index1, random_index2;
         int numCollisions = 0;
         for (int n=0; n<local_Nmax_integer; n++) { // loop over collision for this cell
         
            random_index1 = MathUtils::randInt(0, local_numCell-1);  
            random_index2 = MathUtils::randInt(0, local_numCell-1);  
            while (random_index2==random_index1) { random_index2 = MathUtils::randInt(0, local_numCell-1); }

            // get particle data for first particle    
            JustinsParticlePtr& this_part1 = vector_part_ptrs[random_index1];
            this_part1_ptr = this_part1.getPointer();
            std::array<Real,3>& betap1 = this_part1_ptr->velocity();
            Real& wp1 = this_part1_ptr->weight();
            
            // get particle data for second particle    
            JustinsParticlePtr& this_part2 = vector_part_ptrs[random_index2];
            this_part2_ptr = this_part2.getPointer();
            std::array<Real,3>& betap2 = this_part2_ptr->velocity();
            Real& wp2 = this_part2_ptr->weight();

            // compute relative velocity magnitude
            g12 = 0.0;
            for (int dir=0; dir<3; dir++) { g12 += pow(betap1[dir]-betap2[dir],2); }
            g12 = sqrt(g12)*Constants::CVAC; // relative velocity [m/s]
         
            // determine if the pair collide and then do the collision
            q12 = g12*m_sigmaT/(local_gmax*m_sigmaT);
            Real rand_num2 = MathUtils::rand();
            if (rand_num2>q12) { continue; }

            numCollisions = numCollisions + 1;
               
            //compute deltaU
            R = MathUtils::rand();
            costh = 1.0 - 2.0*R;
            sinth = sqrt(1.0 - costh*costh);
            phi = Constants::TWOPI*MathUtils::rand();
            ScatteringUtils::computeDeltaU(deltaU,betap1,betap2,costh,sinth,phi);
            //deltaU = {0,0,0};

            // update particle velocities
            if (m_weight_method==CONSERVATIVE && wp1!=wp2) {

               if (local_numCell<3) { continue; }
               
               // get third particle from larger weight particle list
               int random_index3 = MathUtils::randInt(0, local_numCell-1);
               while (random_index3==random_index2 || random_index3==random_index1) {
                  random_index3 = MathUtils::randInt(0,local_numCell-1);
               }
               JustinsParticlePtr& this_part3 = vector_part_ptrs[random_index3];
               this_part3_ptr = this_part3.getPointer();
               std::array<Real,3>& betap3 = this_part3_ptr->velocity();
               Real& wp3 = this_part3_ptr->weight();

               if (wp1<wp2) {
                  // scatter the lower weight particle
                  for (int dir=0; dir<3; dir++) { betap1[dir] += 0.5*deltaU[dir]; }
                  std::array<Real,3> betap2p = betap2;
                  for (int dir=0; dir<3; dir++) { betap2p[dir] -= 0.5*deltaU[dir]; }
                   
                  // transform particles 2, 2p, and 3 into two equally-weighted particles such
                  // that momentum, and direction-dependent energy are conserved
                  ScatteringUtils::collapseThreeToTwo( betap2, wp2, betap3, wp3, betap2p, wp1 );
               }
               else {
                  // scatter the lower weight particle
                  std::array<Real,3> betap1p = betap1;
                  for (int dir=0; dir<3; dir++) { betap1p[dir] += 0.5*deltaU[dir]; }
                  for (int dir=0; dir<3; dir++) { betap2[dir] -= 0.5*deltaU[dir]; }
                  
                  // transform particles 1, 1p, and 3 into two equally-weighted particles such
                  // that momentum, and direction-dependent energy are conserved
                  ScatteringUtils::collapseThreeToTwo( betap1, wp1, betap3, wp3, betap1p, wp2 );
               }

            }
            else {
               Real rand_num3 = MathUtils::rand();
               if (rand_num3<=wp2/wp1) {
                  for (int dir=0; dir<3; dir++) { betap1[dir] += 0.5*deltaU[dir]; }
               }
               if (rand_num3<=wp1/wp2) {
                  for (int dir=0; dir<3; dir++) { betap2[dir] -= 0.5*deltaU[dir]; }
               }
            }

         } // end loop over collisions for this cell
         if (procID()==0 && verbosity) { cout << "numCollisions = " << numCollisions << endl; }
         verbosity=0;

      } // end loop over cells

   } // end loop over boxes
   
   this_part1_ptr = NULL;
   this_part2_ptr = NULL;
   this_part3_ptr = NULL;
   delete this_part1_ptr;
   delete this_part2_ptr;
   delete this_part3_ptr;
   
}

void HardSphere::applyInterScattering( PicChargedSpecies&  a_picSpecies1,
                                       PicChargedSpecies&  a_picSpecies2, 
                                 const DomainGrid&  a_mesh,
                                 const Real         a_dt_sec ) const
{
   CH_TIME("HardSphere::applyInterScattering()");
   
   CH_assert(m_species1_name==a_picSpecies1.name());
   CH_assert(m_species2_name==a_picSpecies2.name());

   // predefine some pointers to be used below
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  
   JustinsParticle* this_part3_ptr = NULL;  
   
   // get the assumed fixed cell volume  
   const Real Vc = a_mesh.getMappedCellVolume();

   // define references to picSpecies1
   LevelData<BinFab<JustinsParticlePtr>>& data1_binfab_ptr = a_picSpecies1.partData_binfab();
   const LevelData<FArrayBox>& numberDensity1 = a_picSpecies1.getNumberDensity(true);
   const LevelData<FArrayBox>& energyDensity1 = a_picSpecies1.getEnergyDensity(false);
   
   // define references to picSpecies2
   LevelData<BinFab<JustinsParticlePtr>>& data2_binfab_ptr = a_picSpecies2.partData_binfab();
   const LevelData<FArrayBox>& numberDensity2 = a_picSpecies2.getNumberDensity(true);
   const LevelData<FArrayBox>& energyDensity2 = a_picSpecies2.getEnergyDensity(false);

   // predefine some variables
   Real local_Teff1, local_numberDensity1, local_energyDensity1;
   Real local_Teff2, local_numberDensity2, local_energyDensity2;
   Real local_gmax, local_Nmax, local_Nmax_remainder, local_Nmax_whole;
   int local_numCell1, local_numCell2;
   Real W1, W2, Wmax;
   Real g12, q12;
   
   std::array<Real,3> deltaU;
   Real R, costh, sinth, phi; 
 
   Real local_nuMaxDt;
   
   Real cvacSq = Constants::CVAC*Constants::CVAC;

   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data1_binfab_ptr.disjointBoxLayout();
   DataIterator ditg(grids);

   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      BinFab<JustinsParticlePtr>& thisBinFab1_ptr = data1_binfab_ptr[ditg];
      const FArrayBox& this_numberDensity1 = numberDensity1[ditg];
      const FArrayBox& this_energyDensity1 = energyDensity1[ditg];
     
      BinFab<JustinsParticlePtr>& thisBinFab2_ptr = data2_binfab_ptr[ditg];
      const FArrayBox& this_numberDensity2 = numberDensity2[ditg];
      const FArrayBox& this_energyDensity2 = energyDensity2[ditg];
   
      std::vector<JustinsParticlePtr> vector_part1_ptrs;
      std::vector<JustinsParticlePtr> vector_part2_ptrs;
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over cells
         
         const IntVect ig = gbit(); // grid index
       
         /////////////////////////////////////////////////////////////////
         //
         // compute the number of possible collisions for this cell
         //
          
         // get local density and temperature and compute local gmax
         local_numberDensity1 = this_numberDensity1.get(ig,0);
         local_numberDensity2 = this_numberDensity2.get(ig,0);
         if (local_numberDensity1*local_numberDensity2 == 0.0) { continue; }
         local_energyDensity1 = 0.0;
         local_energyDensity2 = 0.0;
         for (int dir=0; dir<3; dir++) {
            local_energyDensity1 = local_energyDensity1 + this_energyDensity1.get(ig,dir);  
            local_energyDensity2 = local_energyDensity2 + this_energyDensity2.get(ig,dir);  
         }
         local_Teff1 = 2.0/3.0*local_energyDensity1/local_numberDensity1*cvacSq; // M1/me*(V1[m/s])^2
         local_Teff2 = 2.0/3.0*local_energyDensity2/local_numberDensity2*cvacSq; // M2/me*(V2[m/s])^2
         local_gmax = 2.5*sqrt(2.0*max(local_Teff1,local_Teff2)/m_mu); // 5x thermal speed [m/s]

         // compute local nuMax*dt
         local_nuMaxDt = Max(local_numberDensity1,local_numberDensity2)*m_sigmaT*local_gmax*a_dt_sec;         
         if (local_nuMaxDt>10.0) {
            cout << "WARNING: local_nuMaxDt = " << local_nuMaxDt << endl;
         }

         // compute local maximum number of collsions: Nmax = N1*n2*sigmaTmax*gmax*dt
         List<JustinsParticlePtr>& cell_pList1 = thisBinFab1_ptr(ig,0);
         List<JustinsParticlePtr>& cell_pList2 = thisBinFab2_ptr(ig,0);
         local_numCell1 = cell_pList1.length();
         local_numCell2 = cell_pList2.length();
         if (local_numCell1 < 2 && local_numCell2 < 2) { continue; }
    
         // compute the particle weights for each species assuming all particles
         // for a given species have the same weight
         W1 = local_numberDensity1/local_numCell1*Vc;
         W2 = local_numberDensity2/local_numCell2*Vc;
         Wmax = Max(W1,W2); 
         
         // compute the maximum collision number     
         local_Nmax = Wmax*local_numCell1*local_numCell2/Vc*m_sigmaT*local_gmax*a_dt_sec;

         // separate Nmax into whole integer and remainder and set integer Nmax
         local_Nmax_remainder = modf(local_Nmax, &local_Nmax_whole);
         Real rand_num = MathUtils::rand();
         int local_Nmax_integer = static_cast<int>(local_Nmax_whole);
         if (rand_num<local_Nmax_remainder) {
            local_Nmax_integer = local_Nmax_integer + 1;
         }
         
         // copy the iterators to a vector in order to randomly select
         vector_part1_ptrs.clear();
         vector_part2_ptrs.clear();
         vector_part1_ptrs.reserve(local_numCell1);
         vector_part2_ptrs.reserve(local_numCell2);
         ListIterator<JustinsParticlePtr> lit1(cell_pList1);
         ListIterator<JustinsParticlePtr> lit2(cell_pList2);
         for (lit1.begin(); lit1.ok(); ++lit1) { vector_part1_ptrs.push_back(lit1()); }
         for (lit2.begin(); lit2.ok(); ++lit2) { vector_part2_ptrs.push_back(lit2()); }
        
         ////////////////////////////////////////////////////////////////////////
         //
         // loop over Nmax and randomly choose two unique elements from the vector to scatter
         //
         int random_index1, random_index2;
         int numCollisions = 0;
         for (int n=0; n<local_Nmax_integer; n++) { // loop over collision for this cell
         
            random_index1 = MathUtils::randInt(0, local_numCell1-1);  
            random_index2 = MathUtils::randInt(0, local_numCell2-1);  

            // get particle data for first particle    
            JustinsParticlePtr& this_part1 = vector_part1_ptrs[random_index1];
            this_part1_ptr = this_part1.getPointer();
            std::array<Real,3>& betap1 = this_part1_ptr->velocity();
            Real& wp1 = this_part1_ptr->weight();
            
            // get particle data for second particle    
            JustinsParticlePtr& this_part2 = vector_part2_ptrs[random_index2];
            this_part2_ptr = this_part2.getPointer();
            std::array<Real,3>& betap2 = this_part2_ptr->velocity();
            Real& wp2 = this_part2_ptr->weight();

            // compute relative velocity magnitude
            g12 = 0.0;
            for (int dir=0; dir<3; dir++) { g12 = g12 + pow(betap1[dir]-betap2[dir],2); }
            g12 = sqrt(g12)*Constants::CVAC; // relative velocity [m/s]
         
            // determine if the pair collide and then do the collision
            q12 = g12*m_sigmaT/(local_gmax*m_sigmaT);
            Real rand_num2 = MathUtils::rand();
            if (rand_num2>q12) { continue; }

            numCollisions = numCollisions + 1;
               
            //compute deltaU
            R = MathUtils::rand();
            costh = 1.0 - 2.0*R;
            sinth = sqrt(1.0 - costh*costh);
            phi = Constants::TWOPI*MathUtils::rand();
            ScatteringUtils::computeDeltaU(deltaU,betap1,betap2,costh,sinth,phi);
            //deltaU = {0,0,0};
            
            // update particle velocities
            if (m_weight_method==CONSERVATIVE && wp1!=wp2) {

               if (wp1<wp2) {

                  if (local_numCell2<2) { continue; }

                  // scatter the lower weight particle
                  for (int dir=0; dir<3; dir++) { betap1[dir] += 0.5*deltaU[dir]; }
                  std::array<Real,3> betap2p = betap2;
                  for (int dir=0; dir<3; dir++) { betap2p[dir] -= 0.5*deltaU[dir]; }
               
                  // get third particle from larger weight particle list
                  int random_index3 = MathUtils::randInt(0, local_numCell2-1);
                  while (random_index3==random_index2) {
                     random_index3 = MathUtils::randInt(0,local_numCell2-1);
                  }
                  JustinsParticlePtr& this_part3 = vector_part2_ptrs[random_index3];
                  this_part3_ptr = this_part3.getPointer();
                  std::array<Real,3>& betap3 = this_part3_ptr->velocity();
                  Real& wp3 = this_part3_ptr->weight();
             
                  // transform particles 2, 2p, and 3 into two equally-weighted particles such
                  // that momentum, and direction-dependent energy are conserved
                  ScatteringUtils::collapseThreeToTwo( betap2, wp2, betap3, wp3, betap2p, wp1 );

               }
               else {

                  if (local_numCell1<2) { continue; }

                  // scatter the lower weight particle
                  std::array<Real,3> betap1p = betap1;
                  for (int dir=0; dir<3; dir++) { betap1p[dir] += 0.5*deltaU[dir]; }
                  for (int dir=0; dir<3; dir++) { betap2[dir] -= 0.5*deltaU[dir]; }
               
                  // get third particle from larger weight particle list
                  int random_index3 = MathUtils::randInt(0, local_numCell1-1);
                  while (random_index3==random_index1) {
                     random_index3 = MathUtils::randInt(0,local_numCell1-1);
                  }
                  JustinsParticlePtr& this_part3 = vector_part1_ptrs[random_index3];
                  this_part3_ptr = this_part3.getPointer();
                  std::array<Real,3>& betap3 = this_part3_ptr->velocity();
                  Real& wp3 = this_part3_ptr->weight();
                  
                  // transform particles 1, 1p, and 3 into two equally-weighted particles such
                  // that momentum, and direction-dependent energy are conserved
                  ScatteringUtils::collapseThreeToTwo( betap1, wp1, betap3, wp3, betap1p, wp2 );

               }
               
            }
            else {
               Real rand_num3 = MathUtils::rand();
               if (rand_num3<=wp2/wp1) {
                  for (int dir=0; dir<3; dir++) { betap1[dir] +=  m_mu/m_mass1*deltaU[dir]; }
               }
               if (rand_num3<=wp1/wp2) {
                  for (int dir=0; dir<3; dir++) { betap2[dir] -= m_mu/m_mass2*deltaU[dir]; }
               }
            }

         } // end loop over collisions for this cell
         if (procID()==0 && verbosity) { cout << "numCollisions = " << numCollisions << endl; }
         verbosity=0;

      } // end loop over cells

   } // end loop over boxes
   
   this_part1_ptr = NULL;
   this_part2_ptr = NULL;
   this_part3_ptr = NULL;
   delete this_part1_ptr;
   delete this_part2_ptr;
   delete this_part3_ptr;
   
}

#include "NamespaceFooter.H"

