
#include "VariableHardSphere.H"
#include "Constants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"

#include "NamespaceHeader.H"

      
void VariableHardSphere::applySelfScattering( PicSpecies&             a_picSpecies, 
                                        const DomainGrid&             a_mesh,
                                        const LevelData<FArrayBox>&   a_numberDensity,
                                        const LevelData<FArrayBox>&   a_energyDensity,
                                        const Real                    a_dt ) const
{
   CH_TIME("VariableHardSphere::applySelfScattering()");
 
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  
    
   // define refence to a_picSpcies binfab container of pointers to particle data
   LevelData<BinFab<JustinsParticlePtr>>& data_binfab_ptr = a_picSpecies.partData_binfab();

   // loop over lists in each cell and test shuffle
   //
   const DisjointBoxLayout& grids = data_binfab_ptr.disjointBoxLayout();
   DataIterator ditg(grids);
 
   const Real Mass = a_picSpecies.mass(); // species mass in units of electron mass
   Real local_Teff, local_numberDensity, local_energyDensity, local_gmax;
   int local_numCell;
   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes

      const FArrayBox& this_numberDensity = a_numberDensity[ditg];
      const FArrayBox& this_energyDensity = a_energyDensity[ditg];
     
      BinFab<JustinsParticlePtr>& thisBinFab_ptr = data_binfab_ptr[ditg];
   
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices
         
         const IntVect ig = gbit(); // grid index
        
         // get local density and temperature and compute local gmax
         local_numberDensity = this_numberDensity.get(ig,0); 
         local_energyDensity = this_energyDensity.get(ig,0); 
         local_Teff = 2.0/3.0*local_energyDensity/local_numberDensity; // M/me*(V[m/s])^2
         local_gmax = 5.0*sqrt(local_Teff/Mass); // 5x thermal speed [m/s]

         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         local_numCell = cell_pList.length();
         if(local_numCell < 2) break;        
         if(!procID() && verbosity) {
            cout << "JRA: local_numCell = " << local_numCell << endl;
            cout << "JRA: local_numberDensity = " << local_numberDensity << endl;
            cout << "JRA: local_energyDensity = " << local_energyDensity << endl;
            cout << "JRA: local_gmax [m/s] = " << local_gmax << endl;
         }
        
         // copy the iterators to a vector in order to randomly selection
         std::vector<JustinsParticlePtr> vector_part_ptrs;
         vector_part_ptrs.reserve(local_numCell);
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) vector_part_ptrs.push_back(lit());
        
         // loop over number of collisions to perform and randomly choose two unique elements from the vector
         int random_index1, random_index2; 
         for (auto n=0; n<local_numCell; n++) {
            random_index1 = std::rand() % local_numCell;  
            random_index2 = std::rand() % local_numCell;
            while(random_index2==random_index1) random_index2 = std::rand() % local_numCell;  

            // get particle data for first particle    
            JustinsParticlePtr& this_part1 = vector_part_ptrs[random_index1];
            this_part1_ptr = this_part1.getPointer();
            const uint64_t& this_ID1 = this_part1_ptr->ID();
            const RealVect& this_xp1 = this_part1_ptr->position();
            const RealVect& this_vp1 = this_part1_ptr->velocity();
            
            // get particle data for second particle    
            JustinsParticlePtr& this_part2 = vector_part_ptrs[random_index2];
            this_part2_ptr = this_part2.getPointer();
            const uint64_t& this_ID2 = this_part2_ptr->ID();
            const RealVect& this_xp2 = this_part2_ptr->position();
            const RealVect& this_vp2 = this_part2_ptr->velocity();

            if(procID()==0 && verbosity) {
               cout << "JRA random_index1 = " << random_index1 << endl;
               cout << "JRA random_index2 = " << random_index2 << endl;
               cout << "JRA: ID1 = " << this_ID1 << endl;
               cout << "JRA: ID2 = " << this_ID2 << endl;
               cout << "JRA: position1 = " << this_xp1 << endl;
               cout << "JRA: position2 = " << this_xp2 << endl;
            }

         }
         verbosity=0;
      }

   }
   
   this_part1_ptr = NULL;
   this_part2_ptr = NULL;
   delete this_part1_ptr;
   delete this_part2_ptr;

}

#include "NamespaceFooter.H"

