
#include "Elastic.H"
#include "Constants.H"
#include "MathUtils.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"

#include <iostream>
#include <fstream>

#include "NamespaceHeader.H"


void Elastic::initialize( const PicSpeciesPtrVect&  a_picSpeciesPtrVect,
                             const DomainGrid&         a_mesh )
{
   CH_TIME("Elastic::initialize()");
   
   // get pointer to species 1
   CH_assert(m_sp1<a_picSpeciesPtrVect.size());
   PicSpeciesPtr this_species1(a_picSpeciesPtrVect[m_sp1]);

   // get pointer to species 2
   CH_assert(m_sp2<a_picSpeciesPtrVect.size());
   PicSpeciesPtr this_species2(a_picSpeciesPtrVect[m_sp2]);
   
   // set the species names
   m_species1_name = this_species1->name();
   m_species2_name = this_species2->name();

   // set the species masses
   m_mass1 = this_species1->mass();   // species 1 mass / me
   m_mass2 = this_species2->mass();   // species 2 mass / me
   m_m1om2 = m_mass1/m_mass2;         // mass ratio
   m_mu = m_mass1*m_mass2/(m_mass1 + m_mass2); // reduced mass

   // pre-define some constants
   m_mcSq = Constants::ME*Constants::CVAC*Constants::CVAC/Constants::QE; // 0.511e6 
         
   // load cross section from file
   if(m_xsec_type == TEXT_FILE) {
      readCrossSectionFile();
      setOkhrimovskyyXI();
   }

   //
   // set the mean free time
   //
   
   if(this_species1->scatter() && this_species2->scatter()) {
   
      const bool setMoments = true;
      const LevelData<FArrayBox>& numberDensity1 = this_species1->getNumberDensity(setMoments);
      const LevelData<FArrayBox>& energyDensity1 = this_species1->getEnergyDensity(setMoments);
   
      const LevelData<FArrayBox>& numberDensity2 = this_species2->getNumberDensity(setMoments);
      const LevelData<FArrayBox>& energyDensity2 = this_species2->getEnergyDensity(setMoments);
 
      setMeanFreeTime(numberDensity1,energyDensity1,numberDensity2,energyDensity2);

   }

   if (m_verbosity)  printParameters();

}

void Elastic::readCrossSectionFile()
{

   std::ifstream xsecFile;
   std::string this_line;
   bool data_start_line_found=false;
   int this_line_number=0;
   Real this_m1om2;
   Real this_energy;
   Real this_xsec;
      
   xsecFile.open(m_xsec_fname,std::ifstream::in);
   if(xsecFile.is_open()) {
      while(getline(xsecFile, this_line)) {
         stringstream ss;
         ss << this_line;
         if(this_line_number==0) CH_assert(this_line=="elastic" || this_line=="ELASTIC")
         if(this_line_number==1) m_reaction_label = this_line;
         if(this_line_number==2) {
            stringstream(this_line) >> this_m1om2;
            //CH_assert(m_m1om2==this_mM); // don't compare doubles !!!
         }
         if(data_start_line_found && this_line.compare(1,5,"-----")==0) {
            break;
         }
         if(data_start_line_found) {
            string temp;
            ss >> temp; // first element is energy [eV]
            stringstream(temp) >> this_energy;
            m_Evec.push_back(this_energy);
            ss >> temp; // second element is total elastic xsec [m^2]
            stringstream(temp) >> this_xsec;
            m_Qela.push_back(this_xsec);
            ss >> temp; // third element is elastic momentum-tansfer xsec [m^2]
            stringstream(temp) >> this_xsec;
            m_Qelm.push_back(this_xsec);
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
   CH_assert(m_Evec.front()==0.0);
   
}

void Elastic::printCrossSections() const
{
   std::ios_base::fmtflags f(cout.flags());
   cout << std::setprecision(5) << std::scientific << endl;
   cout << "energy [eV]  Qela[m^2]    Qelm[m^2]    scattering xi" << endl;
   for (int i = 0; i < m_Evec.size(); i++) {
      cout << m_Evec[i] << "  ";
      cout << m_Qela[i] << "  ";
      cout << m_Qelm[i] << "  ";
      cout << m_xi[i] << "  ";
      cout << endl;
   }
   cout << endl;  
   cout.flags(f);
}

void Elastic::setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity1,
                               const LevelData<FArrayBox>&  a_energyDensity1,
                               const LevelData<FArrayBox>&  a_numberDensity2,
                               const LevelData<FArrayBox>&  a_energyDensity2 ) const
{
   CH_TIME("Elastic::setMeanFreeTime()");
   
   Real local_numberDensity1, local_betaSqEff1;
   Real local_numberDensity2, local_betaSqEff2;
   Real g12, xi, sigma;
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
         sigma = getSigma(xi, g12);
         
         // compute local local nuMax
         local_nuMax = local_numberDensity2*sigma*g12*Constants::CVAC;

         // update box_nuMax for time step
         box_nuMax = Max(box_nuMax,local_nuMax);
         
      }
   
   }
   
   Real global_nuMax = box_nuMax;
   m_scatter_dt = 1.0/global_nuMax; // mean free time [s]

}

void Elastic::applyScattering( PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                            const DomainGrid&         a_mesh,
                            const Real                a_dt_sec ) const
{
   CH_TIME("Elastic::applyScattering()");
      
   PicSpeciesPtr this_species1(a_pic_species_ptr_vect[m_sp1]);
   PicSpeciesPtr this_species2(a_pic_species_ptr_vect[m_sp2]);
   if(!this_species1->scatter()) return;
   if(!this_species2->scatter()) return;
   
   // electron impact elastic: e + A => e + A
   electronImpact( *this_species1, *this_species2, a_mesh, a_dt_sec );

}
      
void Elastic::electronImpact( PicSpecies&  a_picSpecies1,
                              PicSpecies&  a_picSpecies2, 
                        const DomainGrid&  a_mesh,
                        const Real         a_dt_sec ) const
{
   CH_TIME("Elastic::electronImpact()");
   
   // predefine some pointers to be used below
   JustinsParticle* this_part1_ptr = NULL;  
   JustinsParticle* this_part2_ptr = NULL;  
   JustinsParticle* this_part3_ptr = NULL;  

   // define references to picSpecies1 (primary species)
   LevelData<BinFab<JustinsParticlePtr>>& data1_binfab_ptr = a_picSpecies1.partData_binfab();
   
   // define references to picSpecies2
   LevelData<BinFab<JustinsParticlePtr>>& data2_binfab_ptr = a_picSpecies2.partData_binfab();
   // get/set density here from binfab_ptr, which  doesn't contain particles created
   // during inelastic events, but does know about those that are killed ... 
   const LevelData<FArrayBox>& numberDensity2 = a_picSpecies2.getNumberDensityFromBinFab();
  
   // predefine some variables
   Real local_numberDensity2;
   int local_numCell1, local_numCell2;
   Real g12, arg, q12;
   
   std::array<Real,3> deltaU;
   Real costh, sinth, phi, rand_num, R; 
   Real xi, sigma;
   
   // loop over lists in each cell and test shuffle
   const DisjointBoxLayout& grids = data1_binfab_ptr.disjointBoxLayout();
   DataIterator ditg(grids);

   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes
      
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
            std::array<Real,3>& betap2 = this_part2_ptr->velocity();
            Real& wp2 = this_part2_ptr->weight();

            // compute relative beta magnitude (non-relativistic)
            g12 = 0.0;
            for (int dir=0; dir<3; dir++) g12 = g12 + pow(betap1[dir]-betap2[dir],2);
            g12 = sqrt(g12);

            // get the cross section
            sigma = getSigma(xi, g12);
            if(sigma==0.0) continue;
         
            // determine if the pair collide and then do the collision
            arg = g12*Constants::CVAC*sigma*local_numberDensity2*a_dt_sec;
            q12 = 1.0 - exp(-arg);
            rand_num = MathUtils::rand();
            
            if(rand_num<=q12) { // this pair collides

               numCollisions = numCollisions + 1;

               //compute deltaU
               phi = Constants::TWOPI*MathUtils::rand();
               R = MathUtils::rand();
               //costh = 1.0 - 2.0*R;
               //costh = 1.0 - 2.0*R*(1.0-xi)/(1.0 + xi*(1.0 - 2.0*R));
               costh = ScatteringUtils::getScatteringCos(R,xi);
               sinth = sqrt(1.0 - costh*costh);
               ScatteringUtils::computeDeltaU(deltaU,betap1,betap2,costh,sinth,phi);

               // update particle velocity vectors
               if (m_weight_method==CONSERVATIVE && wp1<wp2) {

                  if(local_numCell2<2) continue;

                  // scatter the lower weight particle
                  for (int dir=0; dir<3; dir++) betap1[dir] += m_mu/m_mass1*deltaU[dir];
                  
                  std::array<Real,3> betap2p = betap2;
                  for (int dir=0; dir<3; dir++) betap2p[dir] -= m_mu/m_mass2*deltaU[dir];
             
                  // get third particle from larger weight particle list
                  int random_index3 = MathUtils::randInt(0, local_numCell2-1);
                  while(random_index3==random_index2) {
                     random_index3 = MathUtils::randInt(0,local_numCell2-1);
                  }
                  JustinsParticlePtr& this_part3 = vector_part2_ptrs[random_index3];
                  this_part3_ptr = this_part3.getPointer();
                  std::array<Real,3>& betap3 = this_part3_ptr->velocity();
                  Real& wp3 = this_part3_ptr->weight();

                  // transform particles 2, 2p, and 3 into two equally-weighted particles such
                  // that momentum, and direction-dependent energy are conserved
                  ScatteringUtils::collapseThreeToTwo( betap2, wp2, betap3, wp3, betap2p, wp1 );
                  
                  // what about when wp1>wp2? Need to randomly get another lit().. 

               }
               else {
                  Real rand_num2 = MathUtils::rand();
                  if(rand_num2<=wp2/wp1) {
                     for (int dir=0; dir<3; dir++) betap1[dir] += m_mu/m_mass1*deltaU[dir];
                  }
                  if(rand_num2<=wp1/wp2) {
                     for (int dir=0; dir<3; dir++) betap2[dir] -= m_mu/m_mass2*deltaU[dir];
                  }
               }

            }

         } // end loop over collisions for this cell
         if(procID()==0 && verbosity) cout << "numCollisions = " << numCollisions << endl;
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

Real Elastic::getSigma( Real&  a_xi, const Real  a_beta ) const
{
   Real sigma = 0.0;
   a_xi = 0;
   if(m_xsec_type==CONSTANT) sigma = m_const_sigma;
   if(m_xsec_type==TEXT_FILE) sigma = getTextSigma(a_xi, a_beta);

   return sigma;
}

Real Elastic::getTextSigma( Real&  a_xi, const Real  a_beta ) const
{
   CH_TIME("Elastic::getTextSigma()");
 
   const Real betaSq = a_beta*a_beta;
   const Real KE_cm = m_mu*m_mcSq*betaSq/2.0; // relative energy in cm frame [eV]

   Real sigma=0.0, xi=0.0;
   const int N = m_Evec.size();
   int index = 0;
   
   if(m_angular_scattering == ISOTROPIC) {

      // use momentum-transfer cross section for isotropic scattering
      if(KE_cm >= m_Evec.back()) { // SigM ~ ln(E)/E
         //sigma = m_Qelm.back(); // JRA, testing E/N with elastic only
         sigma = m_Qelm.back()*log(KE_cm)/log(m_Evec.back())*m_Evec.back()/KE_cm;
      }
      else {

         // find the index
         index = N/2;
         while (KE_cm < m_Evec[index]) index--;
         while (KE_cm > m_Evec[index+1]) index++;

         // get sigma
         if(m_use_loglog_interp && m_Qelm[index]*m_Evec[index]>0.0) {
            sigma = ScatteringUtils::loglogInterp(m_Evec,m_Qelm,KE_cm,index);
         }
         else {
            sigma = ScatteringUtils::linearInterp(m_Evec,m_Qelm,KE_cm,index);
         }

      }
      
   }
   else if(m_angular_scattering == OKHRIMOVSKYY) {

      if(KE_cm >= m_Evec.back()) { // SigT ~ 1/E
         //sigma = m_Qela.back(); // JRA, testing E/N with elastic only
         sigma = m_Qela.back()*m_Evec.back()/KE_cm;
         xi = m_xi.back();
      }
      else {

         // find the index
         index = N/2;
         while (KE_cm < m_Evec[index]) index--;
         while (KE_cm > m_Evec[index+1]) index++;

         // get sigma
         if(m_use_loglog_interp && m_Qela[index]*m_Evec[index]>0.0) {
            sigma = ScatteringUtils::loglogInterp(m_Evec,m_Qela,KE_cm,index);
         }
         else {
            sigma = ScatteringUtils::linearInterp(m_Evec,m_Qela,KE_cm,index);
         }

         // get xi
         if(m_Evec[index]*KE_cm>0.0) {
            xi = ScatteringUtils::semilogInterp(m_Evec,m_xi,KE_cm,index);
         }
         else {
            xi = ScatteringUtils::linearInterp(m_Evec,m_xi,KE_cm,index);
         }

      }

   }
   else {
      MayDay::Error("non-implemented scattering model in Elastic class");
   }
   
   a_xi = xi;
   return sigma;   

}  

void Elastic::setOkhrimovskyyXI()
{
   //  compute Okhrimovskyy xi
   //  sigmaM/sigmaT = (1-xi)/(2*xi^2)*( (1+xi)*ln((1+xi)/(1-xi)) - 2*xi )
   //  used in cos(theta) = 1 - 2*R*(1-xi)/(1+xi*(1-2*R))
 
   Real x, z, RHS, ratio, f, dfdx;
   Real sigmaT, sigmaM, rel_diff;
   int iter;
   const Real rel_tol = 1.0e-10;
   const int iter_max = 100;

   // initialize xi to zeros
   const int N = m_Evec.size();
   m_xi.resize(N,0.0);

   // solve for xi using cross sections
   x = 0.1;
   for(int n=0; n<N; n++) {

      sigmaT = m_Qela.at(n);
      sigmaM = m_Qelm.at(n);
      ratio = sigmaM/sigmaT;
      if (ratio==1.0) continue;

      // initial guess
      if(n>0) x = m_xi.at(n-1);
      if(x==0.0) x = 0.1;

      // use simple Newton method to get xi
      rel_diff = 1.0;
      iter = 0;
      while(rel_diff>rel_tol) {
         z = (1.0 + x)/(1.0 - x);
         RHS = (1.0 - x)/2.0/x/x*( (1.0 + x)*log(z) - 2.0*x );
         f = RHS - ratio;
         dfdx = (1.0 - x)/(2.0*x*x)*(log(z) + 2.0/(1.0 - x) - 2.0) 
              + ((1.0 + x)*log(z) - 2.0*x)*(x - 2.0)/2.0/x/x/x;
         x = x - f/dfdx;
         rel_diff = abs(f/dfdx);
         iter = iter + 1;
         if(iter==iter_max) {
            cout << "WARNING: rel_diff = " << rel_diff << endl;
            cout << "WARNING: iter = iter_max = " << iter_max << endl;
            MayDay::Error( "Elastic scattering xi calc: iter = iter_max" );
            break;
         }
      }
      m_xi.at(n) = x;

   }

}

#include "NamespaceFooter.H"

