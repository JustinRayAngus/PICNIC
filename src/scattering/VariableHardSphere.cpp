
#include "VariableHardSphere.H"
#include "Constants.H"
#include "CodeUnits.H"
#include "MathUtils.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"

#include "NamespaceHeader.H"


void VariableHardSphere::initialize( const DomainGrid&         a_mesh,
                                     const PicSpeciesPtrVect&  a_picSpeciesPtrVect )
{
   CH_TIME("VariableHardSphere::initialize()");
      
   int this_sp = m_sp1;
   CH_assert(m_sp1<a_picSpeciesPtrVect.size());
   PicSpeciesPtr this_picSpecies(a_picSpeciesPtrVect[this_sp]);
   CH_assert(this_picSpecies->scatter()); // assert that collisions are allowed for this species
   CH_assert(this_picSpecies->charge()==0);      
   
   if(m_sp1==m_sp2) { // self scattering

      m_species1_name = this_picSpecies->name();
      m_species2_name = this_picSpecies->name();
      
      m_mass1 = this_picSpecies->mass(); // normalized mass
      m_mass2 = this_picSpecies->mass(); // normalized mass
 
      // define alpha and Aconst
      //
      m_alpha = 4./(2.*m_eta-1.);
      const Real Mass_kg = m_mass1*Constants::ME; // species mass [kg]
      const Real VT0 = sqrt(Constants::KB*m_T0/Mass_kg);     // reference thermal speed [m/s]
      const Real Gamma0 = tgamma(4.-2./m_alpha);
      m_Aconst = 15./32./Gamma0/m_mu0*Mass_kg/sqrt(Constants::PI)*VT0*pow(4.*VT0*VT0,2./m_alpha);

   } 
   else { // inter-species scattering
      cout << " inter-species VHS collisions not yet implemented !!!!!!!!!! " << endl;
   }
   
   //
   // set the mean free time
   //
   
   // define references to picSpecies1
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity1 = this_picSpecies->getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity1 = this_picSpecies->getEnergyDensity(setMoments);
   
   if(m_sp1==m_sp2) {
      setMeanFreeTime(numberDensity1,energyDensity1);
   }
   else {
      // define references to picSpecies2
      //const LevelData<FArrayBox>& numberDensity2 = this_picSpecies2->getNumberDensity(setMoments);
      //const LevelData<FArrayBox>& energyDensity2 = this_picSpecies2->getEnergyDensity(setMoments);
      //setMeanFreeTime(numberDensity1,energyDensity1,numberDensity2,energyDensity2);
   }

   if (m_verbosity)  printParameters();

}


void VariableHardSphere::setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity,
                                          const LevelData<FArrayBox>&  a_energyDensity ) const
{
   CH_TIME("VariableHardSphere::setMeanFreeTime()");
   
   Real mass = m_mass1;
   
   // predefine some variables
   Real local_Teff, local_numberDensity, local_energyDensity, local_VTeff;
   Real local_sigmaTmax, local_nuMax;
   Real fourPiA = 4.*Constants::PI*m_Aconst;
   Real fourOverAlpha = 4.0/m_alpha;   
   Real box_nuMax=0.0; // for scattering time step calculation
 
   Real cvacSq = Constants::CVAC*Constants::CVAC;

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

         // compute local sigmaTmax = sigmaT(VTeff) local nuMax*dt
         local_sigmaTmax = fourPiA*pow(local_VTeff,-fourOverAlpha);
         local_nuMax = local_numberDensity*local_sigmaTmax*local_VTeff;
         box_nuMax = Max(box_nuMax,local_nuMax);
      }
   
   }
   
   // compute mean free time in physical units [s]
   Real global_nuMax = box_nuMax;
#ifdef CH_MPI
   MPI_Allreduce( &box_nuMax, &global_nuMax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
#endif
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]
   
   //const Real tscale = a_units.getScale(a_units.TIME);
   //m_scatter_dt = m_scatter_dt/tscale;

}
/*
void VariableHardSphere::setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity1,
                                          const LevelData<FArrayBox>&  a_energyDensity1,
                                          const LevelData<FArrayBox>&  a_numberDensity2,
                                          const LevelData<FArrayBox>&  a_energyDensity2 )
{
   CH_TIME("VariableHardSphere::setMeanFreeTime()");
   
   Real mass = m_mass1;
   
   // predefine some variables
   Real local_Teff, local_numberDensity, local_energyDensity, local_VTeff;
   Real local_sigmaTmax, local_nuMax;
   Real fourPiA = 4.*Constants::PI*m_Aconst;
   Real fourOverAlpha = 4.0/m_alpha;   
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
         local_Teff = 2.0/3.0*local_energyDensity/local_numberDensity; // M/me*(V[m/s])^2
         local_VTeff = sqrt(local_Teff/mass); // thermal speed [m/s]

         // compute local sigmaTmax = sigmaT(VTeff) local nuMax*dt
         local_sigmaTmax = fourPiA*pow(local_VTeff,-fourOverAlpha);
         local_nuMax = local_numberDensity*local_sigmaTmax*local_VTeff;
         box_nuMax = Max(box_nuMax,local_nuMax);
      }
   
   }
   
   Real global_nuMax = box_nuMax;
#ifdef CH_MPI
   MPI_Allreduce( &box_nuMax, &global_nuMax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
#endif
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]

}
*/
      
void VariableHardSphere::applySelfScattering( PicSpecies&  a_picSpecies, 
                                        const DomainGrid&  a_mesh,
                                        const Real         a_dt_sec ) const
{
   CH_TIME("VariableHardSphere::applySelfScattering()");
 
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  
 
   Real mass = m_mass1;
   Real cvacSq = Constants::CVAC*Constants::CVAC;
   
   // define reference to a_picSpcies binfab container of pointers to particle data
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab_ptr = a_picSpecies.partData_binfab();
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity = a_picSpecies.getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity = a_picSpecies.getEnergyDensity(setMoments);

   // predefine some variables
   Real local_Teff, local_numberDensity, local_energyDensity, local_gmax;
   Real local_sigmaTmax, local_nuMaxDt, local_Nmax, local_Nmax_remainder, local_Nmax_whole;
   int local_numCell;
   Real local_sigmaT, local_nuDt, g12, q12;
   Real fourPiA = 4.*Constants::PI*m_Aconst;
   Real fourOverAlpha = 4.0/m_alpha;   
   
   std::array<Real,3> deltaU;
   Real theta; 
 
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
         if(local_numberDensity == 0.0) continue;
         local_energyDensity = 0.0;
         for( int dir=0; dir<3; dir++) {
            local_energyDensity = local_energyDensity + this_energyDensity.get(ig,dir);  
         }
         local_Teff = 2.0/3.0*local_energyDensity/local_numberDensity*cvacSq; // M/me*(V[m/s])^2
         local_gmax = 5.0*sqrt(local_Teff/mass); // 5x thermal speed [m/s]

         // compute local sigmaTmax = sigmaT(gmax) local nuMax*dt
         local_sigmaTmax = fourPiA*pow(local_gmax,-fourOverAlpha);
         local_nuMaxDt = local_numberDensity*local_sigmaTmax*local_gmax*a_dt_sec;         
         if(local_nuMaxDt>10.0) 
            cout << "WARNING: local_nuMaxDt = " << local_nuMaxDt << endl;

         // compute local maximum number of collsions: Nmax = 1/2*(N-1)*n*sigmaTmax*gmax*dt
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         local_numCell = cell_pList.length();
         if(local_numCell < 2) continue;
         local_Nmax = 0.5*(local_numCell-1)*local_nuMaxDt;
    
         // separate Nmax into whole integer and remainder and set integer Nmax
         local_Nmax_remainder = modf(local_Nmax, &local_Nmax_whole);
         Real rand_num = MathUtils::rand();
         int local_Nmax_integer = static_cast<int>(local_Nmax_whole);
         if(rand_num<local_Nmax_remainder) {
            local_Nmax_integer = local_Nmax_integer + 1;
         }
         
         if(!procID() && verbosity) {
            cout << "JRA: local_numCell = " << local_numCell << endl;
            cout << "JRA: local_nuMaxDt = " << local_nuMaxDt << endl;
            cout << "JRA: local_Nmax    = " << local_Nmax << endl;
            cout << "JRA: local_Nmax_whole = " << local_Nmax_whole << endl;
            cout << "JRA: local_Nmax_remainder = " << local_Nmax_remainder << endl;
            cout << "JRA: local_Nmax_integer = " << local_Nmax_integer << endl;
            cout << "JRA: local_numberDensity = " << local_numberDensity << endl;
            cout << "JRA: local_energyDensity = " << local_energyDensity << endl;
            cout << "JRA: local_gmax [m/s] = " << local_gmax << endl;
         }
        
         // copy the iterators to a vector in order to randomly selection
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
         for (auto n=0; n<local_Nmax_integer; n++) { // loop over collision for this cell
         
            random_index1 = MathUtils::randInt(0, local_numCell-1);  
            random_index2 = MathUtils::randInt(0, local_numCell-1);  
            while(random_index2==random_index1) random_index2 = MathUtils::randInt(0, local_numCell-1);  

            // get particle data for first particle    
            JustinsParticlePtr& this_part1 = vector_part_ptrs[random_index1];
            this_part1_ptr = this_part1.getPointer();
            const uint64_t& this_ID1 = this_part1_ptr->ID();
            const RealVect& this_xp1 = this_part1_ptr->position();
            std::array<Real,3>& this_betap1 = this_part1_ptr->velocity();
            
            // get particle data for second particle    
            JustinsParticlePtr& this_part2 = vector_part_ptrs[random_index2];
            this_part2_ptr = this_part2.getPointer();
            const uint64_t& this_ID2 = this_part2_ptr->ID();
            const RealVect& this_xp2 = this_part2_ptr->position();
            std::array<Real,3>& this_betap2 = this_part2_ptr->velocity();

            // compute relative velocity magnitude
            g12 = 0.0;
            for (int dir=0; dir<3; dir++) {
               g12 = g12 + pow(this_betap1[dir]-this_betap2[dir],2);
            }
            g12 = sqrt(g12)*Constants::CVAC; // relavive velocity [m/s]
         
            // compute local sigmaT = sigmaT(g12) local nu*dt
            local_sigmaT = fourPiA*pow(g12,-fourOverAlpha);
            local_nuDt = local_numberDensity*local_sigmaT*g12*a_dt_sec;         

            // determine if the pair collide and then do the collision
            q12 = g12*local_sigmaT/(local_gmax*local_sigmaTmax);
            rand_num = MathUtils::rand();
            if(rand_num<=q12) { // this pair collides

               numCollisions = numCollisions + 1;
               
               //compute deltaU
               theta = Constants::TWOPI*MathUtils::rand();
               ScatteringUtils::computeDeltaU(deltaU,this_betap1,this_betap2,theta);
               //deltaU = {0,0,0};

               // update particle velocities
               for (int dir=0; dir<3; dir++) {
                  this_betap1[dir] = this_betap1[dir] + 0.5*deltaU[dir];
                  this_betap2[dir] = this_betap2[dir] - 0.5*deltaU[dir];
               }
               if(procID()==0 && verbosity) {
                  cout << "JRA random_index1 = " << random_index1 << endl;
                  cout << "JRA random_index2 = " << random_index2 << endl;
                  cout << "JRA: ID1 = " << this_ID1 << endl;
                  cout << "JRA: ID2 = " << this_ID2 << endl;
                  cout << "JRA: position1 = " << this_xp1 << endl;
                  cout << "JRA: position2 = " << this_xp2 << endl;
               }
            }

         } // end loop over collisions for this cell
         if(procID()==0 && verbosity) cout << "numCollisions = " << numCollisions << endl;
         verbosity=0;

      } // end loop over cells

   } // end loop over boxes
   
   this_part1_ptr = NULL;
   this_part2_ptr = NULL;
   delete this_part1_ptr;
   delete this_part2_ptr;
   
   // While we are here, update the mean free time
   setMeanFreeTime(numberDensity,energyDensity);
   
}

void VariableHardSphere::applyInterScattering( PicSpecies&  a_picSpecies1,
                                               PicSpecies&  a_picSpecies2, 
                                         const DomainGrid&  a_mesh,
                                         const Real         a_dt_sec ) const
{
   CH_TIME("VariableHardSphere::applyInterScattering()");

   cout << "VariableHardSphere::applyInterScattering() NOT YET IMPLEMENTED" << endl;
 
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  
 
   // define references to picSpecies1
   //
   LevelData<BinFab<JustinsParticlePtr>>& data1_binfab_ptr = a_picSpecies1.partData_binfab();
   const bool setMoments = false; // It is the job of the caller to make sure the moments are pre-computed
   const LevelData<FArrayBox>& numberDensity1 = a_picSpecies1.getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity1 = a_picSpecies1.getEnergyDensity(setMoments);
   
   // define references to picSpecies2
   //
   LevelData<BinFab<JustinsParticlePtr>>& data2_binfab_ptr = a_picSpecies2.partData_binfab();
   const LevelData<FArrayBox>& numberDensity2 = a_picSpecies2.getNumberDensity(setMoments);
   const LevelData<FArrayBox>& energyDensity2 = a_picSpecies2.getEnergyDensity(setMoments);

  
}

#include "NamespaceFooter.H"

