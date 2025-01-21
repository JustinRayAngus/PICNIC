#include "MonteCarloNull.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ParticleData.H"
#include "BinFab.H"
#include "ScatteringUtils.H"

#include <iostream>
#include <fstream>

#include "NamespaceHeader.H"


void MonteCarloNull::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                                 const DomainGrid&           a_mesh )
{
   CH_TIME("MonteCarloNull::initialize()");
   
   const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
   const int num_species = a_pic_species_intf.numSpecies();   

   // get pointer to species 1
   CH_assert(m_sp1<num_species);
   const PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];

   // get pointer to species 2
   CH_assert(m_sp2<num_species);
   const PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];

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
   readCrossSectionFile();
   setOkhrimovskyyXI();

   // scattering test
   bool scatteringTest = false;
   if(scatteringTest && !procID()) {

      Real FE_norm = 12.444/m_mcSq;
      Real mass_df1 = 3720.0;
      Real mass_df2 = 3720.0;
      Real mu_df = mass_df1*mass_df2/(mass_df1 + mass_df2);
    

      std::array<Real,3> beta_df1 = {3.0e-5,1.0e-5,-2.0e-5};
      std::array<Real,3> beta_df2 = {2.0e-5,-1.0e-5,1.0e-5};
      //std::array<Real,3> beta_df2 = beta_df1;
      std::array<Real,3> deltaUdf;

      Real energy_before = 0.0;
      for (int dir=0; dir<3; dir++) {
         energy_before += 0.5*beta_df1[dir]*beta_df1[dir]*mass_df1;
         energy_before += 0.5*beta_df2[dir]*beta_df2[dir]*mass_df2;
      }
      cout << std::setprecision(6) << std::scientific << endl;
      cout << " Frag energy [eV] = " << FE_norm*m_mcSq << endl;
      cout << " before: momX = " << beta_df1[0]*mass_df1 + beta_df2[0]*mass_df2 << endl;
      cout << " before: momY = " << beta_df1[1]*mass_df1 + beta_df2[1]*mass_df2 << endl;
      cout << " before: momZ = " << beta_df1[2]*mass_df1 + beta_df2[2]*mass_df2 << endl;
      cout << " before: energy [eV] = " << energy_before*m_mcSq << endl;
      std::cout << setprecision(2) << fixed;
                  
      //compute scattering angles
      Real phi_df = Constants::TWOPI*MathUtils::rand();
      Real costh_df = 1.0 - 2.0*MathUtils::rand();
      Real sinth_df = sqrt(1.0 - costh_df*costh_df);
      ScatteringUtils::computeDeltaU(deltaUdf,beta_df1,beta_df2,
                                    -FE_norm,mu_df,costh_df,sinth_df,phi_df);
      for (int dir=0; dir<3; dir++) beta_df1[dir] += mu_df/mass_df1*deltaUdf[dir];
      for (int dir=0; dir<3; dir++) beta_df2[dir] -= mu_df/mass_df2*deltaUdf[dir];
      
      Real energy_after = 0.0;
      for (int dir=0; dir<3; dir++) {
         energy_after += 0.5*beta_df1[dir]*beta_df1[dir]*mass_df1;
         energy_after += 0.5*beta_df2[dir]*beta_df2[dir]*mass_df2;
      }
      cout << " after: momX = " << beta_df1[0]*mass_df1 + beta_df2[0]*mass_df2 << endl;
      cout << " after: momY = " << beta_df1[1]*mass_df1 + beta_df2[1]*mass_df2 << endl;
      cout << " after: momZ = " << beta_df1[2]*mass_df1 + beta_df2[2]*mass_df2 << endl;
      cout << " after: energy [eV] = " << energy_after*m_mcSq << endl;
      cout << " energy diff [eV] = " << (energy_after - energy_before)*m_mcSq << endl;
      cout << endl;

   }

   // check mass and charge conservation for inelastic reactions
   for (int n=0; n<m_num_exc; n++) {
      int spc = m_exc_species[n];
      if (spc<0) { continue; }
      const PicChargedSpeciesPtr exc_species = pic_species_ptr_vect[species_map[spc]];
      Real this_mass = exc_species->mass();
      int this_q = exc_species->charge();
      CH_assert(this_mass==m_mass2);
      CH_assert(this_q==this_species2->charge());
   }
   
   for (int n=0; n<m_num_dis; n++) {
      int spc1 = m_dis_species1[n];
      int spc2 = m_dis_species2[n];
      if (spc1<0 || spc2<0) { continue; }
      const PicChargedSpeciesPtr dis_species1 = pic_species_ptr_vect[species_map[spc1]];
      const PicChargedSpeciesPtr dis_species2 = pic_species_ptr_vect[species_map[spc2]];
      Real this_mass1 = dis_species1->mass();
      Real this_mass2 = dis_species2->mass();
      Real mass_tot = this_mass1 + this_mass2;
      int this_q1 = dis_species1->charge();
      int this_q2 = dis_species2->charge();
      int q_tot = this_q1 + this_q2;
      int spc3 = m_dis_species3[n];
      if (spc3>=0) {
         const PicChargedSpeciesPtr dis_species3 = pic_species_ptr_vect[species_map[spc3]];
         Real this_mass3 = dis_species3->mass();
         int this_q3 = dis_species3->charge();
         q_tot += this_q3;
         mass_tot += this_mass3;
      }
      if(m_dis_type[n]==DIS_RECOMBINATION) {
         q_tot -= this_species1->charge();
         mass_tot -= this_species1->mass();
      }
      if(mass_tot!=this_species2->mass()) {
         cout << "JRA: this_species1_mass = " << this_species1->mass() << endl;
         cout << "JRA: this_species2_mass = " << this_species2->mass() << endl;
         cout << "JRA: dis_m1 = " << this_mass1 << endl;
         cout << "JRA: dis_m2 = " << this_mass2 << endl;
         cout << "JRA: m_tot = " << mass_tot << endl;
         cout << "JRA: spc3 = " << spc3 << endl;
         CH_assert(mass_tot==m_mass2);
      }
      if(q_tot!=this_species2->charge()) {
         cout << "JRA: this_species1_charge= " << this_species1->charge() << endl;
         cout << "JRA: this_species2_charge= " << this_species2->charge() << endl;
         cout << "JRA: dis_q1 = " << this_q1 << endl;
         cout << "JRA: dis_q2 = " << this_q2 << endl;
         cout << "JRA: q_tot = " << q_tot << endl;
         cout << "JRA: spc3 = " << spc3 << endl;
         CH_assert(q_tot==this_species2->charge()); 
      }
   }

   for (int n=0; n<m_num_izn; n++) {
      int spcE = m_izn_speciesE[n];
      int spcI = m_izn_speciesI[n];
      if (spcE<0 || spcI<0) { continue; }
      const PicChargedSpeciesPtr izn_speciesE = pic_species_ptr_vect[species_map[spcE]];
      const PicChargedSpeciesPtr izn_speciesI = pic_species_ptr_vect[species_map[spcI]];
      Real this_massE = izn_speciesE->mass();
      Real this_massI = izn_speciesI->mass();
      int this_qE = izn_speciesE->charge();
      int this_qI = izn_speciesI->charge();
      CH_assert(this_massE+this_massI==m_mass2);
      CH_assert(this_massE+this_massI==m_mass2);
      CH_assert(this_qE+this_qI==this_species2->charge());
   }

   if (m_verbosity)  printParameters();

}

void MonteCarloNull::readCrossSectionFile()
{

   std::ifstream xsecFile;
   std::string this_line;
   bool data_start_line_found=false;
   int this_line_number=0;
   Real this_m1om2;
   Real this_energy;
   Real this_xsec;

   int num_exc, num_dis, num_izn;
   int this_exc_type;
   int this_dis_type;
   double this_Uexc;     
   double this_Udis;     
   double this_Uizn;     

   xsecFile.open(m_xsec_fname,std::ifstream::in);
   if(xsecFile.is_open()) {
      while(getline(xsecFile, this_line)) {
         stringstream ss;
         ss << this_line;
         if(this_line_number==0) CH_assert(this_line=="monte carlo null" 
                                        || this_line=="MONTE CARLO NULL")
         if(this_line_number==1) m_reaction_label = this_line;
         if(this_line.compare(0,3,"m/M")==0) {
            string temp;
            while (ss.good()) ss >> temp;
            stringstream(temp) >> this_m1om2;
            //cout << "JRA: this_m1om2 = " << this_m1om2 << endl;
         }
         
         if(this_line.compare(0,17,"elastic count = 0")==0) {
            m_no_elastic = true;
         }

         if(this_line.compare(0,16,"excitation count")==0) {
            string temp;
            while (ss.good()) ss >> temp;
            stringstream(temp) >> num_exc;
            CH_assert(num_exc==m_num_exc);
            m_Qexc.resize(num_exc);
            m_Qexc0.resize(num_exc);
            //if(!procID()) cout << "JRA: num_exc = " << num_exc << endl;
         }
         if(this_line.compare(0,15,"excitation type")==0) {
            string temp;
            while (temp.compare("=")!=0) ss >> temp;
            for (int n=0; n<num_exc; n++) {
               ss >> temp;
               stringstream(temp) >> this_exc_type;
               m_exc_type.push_back(static_cast<EXCITATION_TYPE>(this_exc_type));
               //if(!procID()) cout << "JRA: exc_type = " << m_exc_type.back() << endl;
            }
         }
         if(this_line.compare(0,20,"excitation potential")==0) {
            string temp;
            while (temp.compare("=")!=0) ss >> temp;
            for (int n=0; n<num_exc; n++) {
               ss >> temp;
               stringstream(temp) >> this_Uexc;
               m_Uexc.push_back(this_Uexc);
               //if(!procID()) cout << "JRA: Uexc = " << m_Uexc.back() << endl;
            }
         }

         if(this_line.compare(0,18,"dissociation count")==0) {
            string temp;
            while (ss.good()) ss >> temp;
            stringstream(temp) >> num_dis;
            CH_assert(num_dis==m_num_dis);
            m_Qdis.resize(num_dis);
            m_Qdis0.resize(num_dis);
            //if(!procID()) cout << "JRA: num_dis = " << num_dis << endl;
         }
         if(this_line.compare(0,17,"dissociation type")==0) {
            string temp;
            while (temp.compare("=")!=0) ss >> temp;
            for (int n=0; n<num_dis; n++) {
               ss >> temp;
               stringstream(temp) >> this_dis_type;
               m_dis_type.push_back(static_cast<DISSOCIATION_TYPE>(this_dis_type));
               //if(!procID()) cout << "JRA: dis_type = " << m_dis_type.back() << endl;
            }
         }
         if(this_line.compare(0,22,"dissociation potential")==0) {
            string temp;
            while (temp.compare("=")!=0) ss >> temp;
            for (int n=0; n<num_dis; n++) {
               ss >> temp;
               stringstream(temp) >> this_Udis;
               m_Udis.push_back(this_Udis);
               //if(!procID()) cout << "JRA: Udis = " << m_Udis.back() << endl;
            }
         }
         if(this_line.compare(0,20,"diss-fragment energy")==0) {
            string temp;
            while (temp.compare("=")!=0) ss >> temp;
            for (int n=0; n<num_dis; n++) {
               ss >> temp;
               Real this_FEdis;
               stringstream(temp) >> this_FEdis;
               CH_assert(this_FEdis>=0.0);
               m_FEdis.push_back(this_FEdis);
               //if(!procID()) cout << "JRA: FEdis = " << m_FEdis.back() << endl;
            }
         }

         if(this_line.compare(0,16,"ionization count")==0) {
            string temp;
            while (ss.good()) ss >> temp;
            stringstream(temp) >> num_izn;
            CH_assert(num_izn==m_num_izn);
            m_Qizn.resize(num_izn);
            m_Qizn0.resize(num_izn);
            //if(!procID()) cout << "JRA: num_izn = " << num_izn << endl;
         }
         if(this_line.compare(0,20,"ionization potential")==0) {
            string temp;
            while (temp.compare("=")!=0) ss >> temp;
            for (int n=0; n<num_izn; n++) {
               ss >> temp;
               stringstream(temp) >> this_Uizn;
               m_Uizn.push_back(this_Uizn);
               //if(!procID()) cout << "JRA: Uizn = " << m_Uizn.back() << endl;
            }
         }
         
         if(data_start_line_found && this_line.compare(1,5,"-----")==0) {
            break;
         }
         if(data_start_line_found) {
            string temp;
            ss >> temp; // first element is energy [eV]
            stringstream(temp) >> this_energy;
            m_Evec.push_back(this_energy);

            if(m_no_elastic) {
               m_Qela.push_back(0.0);
               m_Qelm.push_back(0.0);
            }
            else {
               ss >> temp; // second element is total elastic xsec [m^2]
               stringstream(temp) >> this_xsec;
               m_Qela.push_back(this_xsec);
               ss >> temp; // third element is elastic-momentum transfer xsec [m^2]
               stringstream(temp) >> this_xsec;
               m_Qelm.push_back(this_xsec);
            } 

            for (int n=0; n<num_exc; n++) {
               ss >> temp;
               stringstream(temp) >> this_xsec;
               m_Qexc.at(n).push_back(this_xsec);
            }
            
            for (int n=0; n<num_dis; n++) {
               ss >> temp;
               stringstream(temp) >> this_xsec;
               m_Qdis.at(n).push_back(this_xsec);
            }
            
            for (int n=0; n<num_izn; n++) {
               ss >> temp;
               stringstream(temp) >> this_xsec;
               m_Qizn.at(n).push_back(this_xsec);
            }
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

   // set the total cross section with ground state
   double sigma_tot;
   for (int i = 0; i < m_Evec.size(); i++) {
      if(m_angular_scattering == ISOTROPIC) sigma_tot = m_Qelm[i];
      else sigma_tot = m_Qela[i];
      for (int n = 0; n < m_Qexc.size(); n++) sigma_tot += m_Qexc[n][i];
      sigma_tot +=  m_Qizn[0][i];
      //m_Qtot.push_back(sigma_tot);
   }

}

void MonteCarloNull::printCrossSections() const
{
   std::ios_base::fmtflags f(cout.flags());
   cout << std::setprecision(5) << std::scientific << endl;
   cout << "energy [eV]  Qela[m^2]    Qelm[m^2]    scattering xi" << endl;
   for (int i = 0; i < m_Evec.size(); i++) {
      cout << m_Evec[i] << "  ";
      cout << m_Qela[i] << "  ";
      cout << m_Qelm[i] << "  ";
      cout << m_xi[i] << "  ";
      for (int n=0; n<m_Qexc.size(); n++) cout << m_Qexc[n][i] << "  ";
      for (int n=0; n<m_Qdis.size(); n++) cout << m_Qdis[n][i] << "  ";
      for (int n=0; n<m_Qizn.size(); n++) cout << m_Qizn[n][i] << "  ";
      cout << endl;
   }
   cout << endl;
   cout.flags(f);
}

void MonteCarloNull::setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const 
{
   CH_TIME("MonteCarloNull::setMeanFreeTime()");
   
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

void MonteCarloNull::setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                                  const LevelData<FArrayBox>&  a_energyDensity1,
                                  const LevelData<FArrayBox>&  a_numberDensity2,
                                  const LevelData<FArrayBox>&  a_energyDensity2 ) const
{
   CH_TIME("MonteCarloNull::setInterMFT()");
  
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

void MonteCarloNull::applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                                const DomainGrid&           a_mesh,
                                const Real                  a_dt_sec ) const
{
   CH_TIME("MonteCarloNull::applyScattering()");
   
   PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getChargedPtrVect();
   const std::vector<int>& species_map = a_pic_species_intf.getSpeciesMap(); 
   const int num_species = a_pic_species_intf.numSpecies();   

   PicChargedSpeciesPtr this_species1 = pic_species_ptr_vect[species_map[m_sp1]];
   PicChargedSpeciesPtr this_species2 = pic_species_ptr_vect[species_map[m_sp2]];
   if (!this_species1->scatter()) { return; }
   if (!this_species2->scatter()) { return; }

   // set vector of pointers to species created by excitation
   std::vector<PicChargedSpeciesPtr> exc_species(m_num_exc);
   for (int n=0; n<m_num_exc; n++) {
      if (m_exc_species.at(n)<num_species && m_exc_species.at(n)>=0) {
         exc_species.at(n) = pic_species_ptr_vect[species_map[m_exc_species.at(n)]];
      }
   }   
   
   // set vector of pointers to species created by dissociation
   std::vector<PicChargedSpeciesPtr> dis_species1(m_num_dis);
   std::vector<PicChargedSpeciesPtr> dis_species2(m_num_dis);
   std::vector<PicChargedSpeciesPtr> dis_species3(m_num_dis);
   for (int n=0; n<m_num_dis; n++) {
      if (m_dis_species1.at(n)<num_species && m_dis_species1.at(n)>=0) {
         dis_species1.at(n) = pic_species_ptr_vect[species_map[m_dis_species1.at(n)]];
         dis_species2.at(n) = pic_species_ptr_vect[species_map[m_dis_species2.at(n)]];
      }
      if (m_dis_species3.at(n)<num_species && m_dis_species3.at(n)>=0) {
         dis_species3.at(n) = pic_species_ptr_vect[species_map[m_dis_species3.at(n)]];
      }
   }

   // set vector of pointers to species created by ionization
   std::vector<PicChargedSpeciesPtr> izn_speciesE(m_num_izn);
   std::vector<PicChargedSpeciesPtr> izn_speciesI(m_num_izn);
   for (int n=0; n<m_num_izn; n++) {
      if (m_izn_speciesI.at(n)<num_species && m_izn_speciesI.at(n)>=0) {
         izn_speciesE.at(n) = pic_species_ptr_vect[species_map[m_izn_speciesE.at(n)]];
         izn_speciesI.at(n) = pic_species_ptr_vect[species_map[m_izn_speciesI.at(n)]];
      }
   }   
   
   // electron impact
   electronImpact( *this_species1, *this_species2, exc_species,
                   dis_species1, dis_species2, dis_species3,
                   izn_speciesE, izn_speciesI, a_mesh, a_dt_sec );

}
      
void MonteCarloNull::electronImpact( PicChargedSpecies&  a_picSpecies1,
                                     PicChargedSpecies&  a_picSpecies2, 
                                     std::vector<PicChargedSpeciesPtr>&  a_exc_species, 
                                     std::vector<PicChargedSpeciesPtr>&  a_dis_species1, 
                                     std::vector<PicChargedSpeciesPtr>&  a_dis_species2, 
                                     std::vector<PicChargedSpeciesPtr>&  a_dis_species3, 
                                     std::vector<PicChargedSpeciesPtr>&  a_izn_speciesE, 
                                     std::vector<PicChargedSpeciesPtr>&  a_izn_speciesI, 
                               const DomainGrid&  a_mesh,
                               const Real         a_dt_sec ) const
{
   CH_TIME("MonteCarloNull::electronImpact()");
   
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
   const LevelData<FArrayBox>& numberDensity2 = a_picSpecies2.getNumberDensityFromBinFab(true);

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

      // initilize lists needed to create new particles
      List<JustinsParticle> new_pList1;
      std::vector<List<JustinsParticle>> new_exc_pList(m_num_exc);
      std::vector<List<JustinsParticle>> new_dis_pList1(m_num_dis);
      std::vector<List<JustinsParticle>> new_dis_pList2(m_num_dis);
      std::vector<List<JustinsParticle>> new_dis_pList3(m_num_dis);
      std::vector<List<JustinsParticle>> new_izn_pListE(m_num_izn);
      std::vector<List<JustinsParticle>> new_izn_pListI(m_num_izn);

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
            
            Real wpNew;

            // compute relative beta magnitude (non-relativistic)
            g12 = 0.0;
            for (int dir=0; dir<3; dir++) g12 += pow(betap1[dir]-betap2[dir],2);
            g12 = sqrt(g12);

            // get the total cross section
            sigma = getSigma(xi, g12);
            arg = g12*Constants::CVAC*sigma*local_numberDensity2*a_dt_sec;
            if(arg==0.0) continue;
         
            // determine if the pair collide
            q12 = 1.0 - exp(-arg);
            rand_num = MathUtils::rand();
            if(rand_num>q12) continue;

            numCollisions = numCollisions + 1;

            // deterime what type of collision occurs
            collisionInfo info;
            info = getCollisionInfo( sigma );
            CH_assert(info.nR>=0);

            //compute scattering angles
            phi = Constants::TWOPI*MathUtils::rand();
            R = MathUtils::rand();
            if(m_angular_scattering == ISOTROPIC) costh = 1.0 - 2.0*R;
            else costh = ScatteringUtils::getScatteringCos(R,xi);
            sinth = sqrt(1.0 - costh*costh);

            // compute deltaU
            if(info.type == NULL_ELASTIC) {
               ScatteringUtils::computeDeltaU(deltaU,betap1,betap2,costh,sinth,phi);
            }
            else {
               ScatteringUtils::computeDeltaU(deltaU,betap1,betap2,
                                              info.Unorm,m_mu,costh,sinth,phi);
            }
            
            //
            //  excitation collision: e + A ==> e + A*
            //

            if(info.type == NULL_EXCITATION && !a_exc_species[info.nR].isNull()) {

               // adjust weights to account for new excitation species
               if(wp1<=wp2) {
                  wpNew = wp1;
                  wp2 = wp2 - wp1;
               }
               else { // split primary and kill secondary
                  wpNew = wp2;
                  Real wpne = wp1 - wp2; // new electron to create
                  wp1 = wp2;
                  wp2 = 0.0; // will be tagged to kill below
                  if(wpne>=1.0e-8*wpNew) {
                     const RealVect& Xpe = this_part1_ptr->position();
                     JustinsParticle particleNE(wpne, Xpe, betap1);
                     new_pList1.append(particleNE);
                  }
               }

               // update velocity of primary
               for (int dir=0; dir<3; dir++) betap1[dir] += m_mu/m_mass1*deltaU[dir];
                  
               // update velocity of new excitation particle
               std::array<Real,3> betapExc = betap2;
               for (int dir=0; dir<3; dir++) betapExc[dir] -= m_mu/m_mass2*deltaU[dir]; 

               // create new excitation particle
               const RealVect& Xp = this_part2_ptr->position();
               JustinsParticle particleExc(wpNew, Xp, betapExc);

               // append new particles to the appropriate lists
               new_exc_pList[info.nR].append(particleExc);

               // lower the density by wpNew for next collision probability
               // local_numberDensity2 -= wpNew/cell_volume; // norm? geom?

            }
            
            //
            //  dissociation collision: e + A  ==> e + A1 + A2      ( DIS_EXCITATION )
            //                       or e + A  ==> e + A1 + A2+ + e ( DIS_IONIZATION )
            //                       or e + A+ ==> A1 + A2  ( DIS_RECOMBINATION )
            //

            else if(info.type == NULL_DISSOCIATION && !a_dis_species1[info.nR].isNull()) {

            if (m_dis_type[info.nR] == DIS_EXCITATION || 
                m_dis_type[info.nR] == DIS_IONIZATION) {

               // adjust weights to account for new dissociation species
               if(wp1<=wp2) {
                  wpNew = wp1;
                  wp2 = wp2 - wp1;
               }
               else { // split primary and kill secondary
                  wpNew = wp2;
                  Real wpne = wp1 - wp2; // new electron to create
                  wp1 = wp2;
                  wp2 = 0.0; // will be tagged to kill below
                  if(wpne>=1.0e-8*wpNew) {
                     const RealVect& Xpe = this_part1_ptr->position();
                     JustinsParticle particleNE(wpne, Xpe, betap1);
                     new_pList1.append(particleNE);
                  }
               }

               // update velocity of primary
               for (int dir=0; dir<3; dir++) betap1[dir] += m_mu/m_mass1*deltaU[dir];
                  
               // update velocity of new dissociation pair as in excitation
               std::array<Real,3> betapDis = betap2;
               for (int dir=0; dir<3; dir++) betapDis[dir] -= m_mu/m_mass2*deltaU[dir]; 

               if(m_dis_type[info.nR] == DIS_EXCITATION) {

                  // create new dissociation particles
                  const RealVect& Xp = this_part2_ptr->position();
                  JustinsParticle particleDis1(wpNew, Xp, betapDis);
                  JustinsParticle particleDis2(wpNew, Xp, betapDis);

                  // append new particles to the appropriate lists
                  new_dis_pList1[info.nR].append(particleDis1);
                  new_dis_pList2[info.nR].append(particleDis2);

               }

               // set beta vector of diss fragments 1 and 2   
               std::array<Real,3> beta_df1 = betapDis;
               std::array<Real,3> beta_df2 = betapDis;

               if(info.FEnorm>0.0) { // add energy and scatter fragments 1 and 2
         
                  Real mass_df1 = a_dis_species1[info.nR]->mass();
                  Real mass_df2 = a_dis_species2[info.nR]->mass();
                  Real mu_df = mass_df1*mass_df2/(mass_df1 + mass_df2);

                  std::array<Real,3> deltaUdf = {0.0,0.0,0.0};
                  
                  //compute scattering angles
                  Real phi_df = Constants::TWOPI*MathUtils::rand();
                  Real costh_df = 1.0 - 2.0*MathUtils::rand();
                  Real sinth_df = sqrt(1.0 - costh_df*costh_df);
                  ScatteringUtils::computeDeltaU(deltaUdf,beta_df1,beta_df2,
                                                -info.FEnorm,mu_df,costh_df,sinth_df,phi_df);
                 for (int dir=0; dir<3; dir++) beta_df1[dir] += mu_df/mass_df1*deltaUdf[dir];
                 for (int dir=0; dir<3; dir++) beta_df2[dir] -= mu_df/mass_df2*deltaUdf[dir];

               }
                  
               // create new dissociation particles
               const RealVect& Xp = this_part2_ptr->position();
               JustinsParticle particleDis1(wpNew, Xp, beta_df1);
               JustinsParticle particleDis2(wpNew, Xp, beta_df2);
                  
               // append new particles to the appropriate lists
               new_dis_pList1[info.nR].append(particleDis1);
               new_dis_pList2[info.nR].append(particleDis2);
               
               if(m_dis_type[info.nR] == DIS_IONIZATION) {
                  JustinsParticle particleDis3(wpNew, Xp, betapDis);
                  new_dis_pList3[info.nR].append(particleDis3);
               }

               // lower the density by wpNew for next collision probability
               // local_numberDensity2 -= wpNew/cell_volume; // norm? geom?

            }

            if (m_dis_type[info.nR] == DIS_RECOMBINATION ) { 
               //cout << "JRA: DIS_RECOMBINATION...." << endl;
               //need to think about charge conservation......
            }

            }

            //
            //  ionization collision: e + A ==> e + A+ + e
            //
 
            else if(info.type == NULL_IONIZATION && !a_izn_speciesI[info.nR].isNull()) {

               // adjust weights to account for ionization
               if(wp1<=wp2) {
                  wpNew = wp1;
                  wp2 = wp2 - wp1;
               }
               else { // split primary and kill secondary
                  wpNew = wp2;
                  Real wpne = wp1 - wp2; // new electron to create
                  wp1 = wp2;
                  wp2 = 0.0; // will be tagged to kill below
                  if(wpne>=1.0e-8*wpNew) {
                     const RealVect& Xpe = this_part1_ptr->position();
                     JustinsParticle particleNE(wpne, Xpe, betap1);
                     new_pList1.append(particleNE);
                  }
               }

               // update velocity of primary
               for (int dir=0; dir<3; dir++) betap1[dir] += m_mu/m_mass1*deltaU[dir];
                  
               // update velocity of new ele/ion pair as in excitation
               std::array<Real,3> betapEI = betap2;
               for (int dir=0; dir<3; dir++) betapEI[dir] -= m_mu/m_mass2*deltaU[dir]; 

               // create new ele/ion pair
               const RealVect& Xp = this_part2_ptr->position();
               JustinsParticle particleE(wpNew, Xp, betapEI);
               JustinsParticle particleI(wpNew, Xp, betapEI);

               // append new particles to the appropriate lists
               new_izn_pListE[info.nR].append(particleE);
               new_izn_pListI[info.nR].append(particleI);

               // lower the density by wpNew for next collision probability
               // local_numberDensity2 -= wpNew/cell_volume; // norm? geom?

            }

            //
            //  elastic collision: e + A ==> e + A
            //

            else {

               if (m_weight_method==CONSERVATIVE && wp1<wp2 && local_numCell2>1) {

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
               
            // update inelastic energy probes
            if(info.type==NULL_EXCITATION) { 
               m_deltaE_exc += m_Uexc.at(info.nR)*wp1*Constants::JOULE_PER_EV; // Joules
            }
            if(info.type==NULL_DISSOCIATION && m_dis_type.at(info.nR)!=DIS_RECOMBINATION) { 
               m_deltaE_exc += (m_Udis.at(info.nR)-m_FEdis.at(info.nR))*wp1*Constants::JOULE_PER_EV; // Joules
            }
            if(info.type==NULL_IONIZATION) { 
               m_deltaE_izn += m_Uizn.at(info.nR)*wp1*Constants::JOULE_PER_EV; // Joules
            }
               
            // remove particles from list that are to be killed
            if(wp1<1.0e-8*wpNew) {
               this_part1_ptr->setKillTag();
               cell_pList1.remove(this_part1_ptr); // JRA, not sure this works !
               local_numCell1 = local_numCell1 - 1;
            }
            if(wp2<1.0e-8*wpNew) {
               vector_part2_ptrs.erase(vector_part2_ptrs.begin()+random_index2);
               this_part2_ptr->setKillTag();
               cell_pList2.remove(this_part2_ptr); // JRA, not sure this works !
               local_numCell2 = local_numCell2 - 1;
            }
            if(local_numCell1*local_numCell2==0) break;

         } // end loop over collisions for this cell
         if(procID()==0 && verbosity) cout << "numCollisions = " << numCollisions << endl;
         verbosity=0;

      } // end loop over cells
      
      // add the new particles from excitation to the appropriate containers
      for (int n=0; n<m_num_exc; n++) {
         if(a_exc_species.at(n).isNull() || new_exc_pList.at(n).length()==0) continue;
         ParticleData<JustinsParticle>& pDataExc = a_exc_species.at(n)->partData();
         pDataExc[ditg].addItemsDestructive(new_exc_pList.at(n));
         CH_assert(new_exc_pList.at(n).length()==0);
      }
      
      // add the new particles from dissociation to the appropriate containers
      for (int n=0; n<m_num_dis; n++) {
         if(a_dis_species1.at(n).isNull() || new_dis_pList1.at(n).length()==0) continue;
         ParticleData<JustinsParticle>& pDataDis1 = a_dis_species1.at(n)->partData();
         ParticleData<JustinsParticle>& pDataDis2 = a_dis_species2.at(n)->partData();
         pDataDis1[ditg].addItemsDestructive(new_dis_pList1.at(n));
         pDataDis2[ditg].addItemsDestructive(new_dis_pList2.at(n));
         CH_assert(new_dis_pList1.at(n).length()==0);
         CH_assert(new_dis_pList2.at(n).length()==0);
         if(a_dis_species3.at(n).isNull() || new_dis_pList3.at(n).length()==0) continue;
         ParticleData<JustinsParticle>& pDataDis3 = a_dis_species3.at(n)->partData();
         pDataDis3[ditg].addItemsDestructive(new_dis_pList3.at(n));
         CH_assert(new_dis_pList3.at(n).length()==0);
      }
      
      // add the new particles from ionization to the appropriate containers
      for (int n=0; n<m_num_izn; n++) {
         if(a_izn_speciesI.at(n).isNull() || new_izn_pListI.at(n).length()==0) continue;
         ParticleData<JustinsParticle>& pDataE = a_izn_speciesE.at(n)->partData();
         ParticleData<JustinsParticle>& pDataI = a_izn_speciesI.at(n)->partData();
         pDataE[ditg].addItemsDestructive(new_izn_pListE.at(n));
         pDataI[ditg].addItemsDestructive(new_izn_pListI.at(n));
         CH_assert(new_izn_pListE.at(n).length()==0);
         CH_assert(new_izn_pListI.at(n).length()==0);
      }

      // add or remove particles from species 1 if tagged to kill
      ParticleData<JustinsParticle>& pData1 = a_picSpecies1.partData();
      pData1[ditg].addItemsDestructive(new_pList1);
      CH_assert(new_pList1.length()==0);

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
   this_part3_ptr = NULL;
   delete this_part1_ptr;
   delete this_part2_ptr;
   delete this_part3_ptr;
   
}

MonteCarloNull::collisionInfo MonteCarloNull::getCollisionInfo( const Real  a_sigma ) const
{ 

   // deterime what type of collision occurs and set the 
   // reaction number nR and normalized potential Unorm
   COLLISION_TYPE type;
   int nR = -1;
   Real Unorm = 0.0;
   Real FEnorm = 0.0;
   Real rand_num = MathUtils::rand();
   collisionInfo info;
   
   Real Qsum = m_Qela0;
   if(rand_num<=Qsum/a_sigma) {
      type = NULL_ELASTIC;
      nR = 0;
   }
   if(nR<0) { 
      for (int n=0; n<m_Qexc0.size(); n++) {
         Qsum += m_Qexc0.at(n);
         if(rand_num<=Qsum/a_sigma) {
            type = NULL_EXCITATION;
            nR = n;
            Unorm = m_Uexc.at(n)/m_mcSq;
            break;
         }
      }
   }
   if(nR<0) { 
      for (int n=0; n<m_Qdis0.size(); n++) {
         Qsum += m_Qdis0.at(n);
         if(rand_num<=Qsum/a_sigma) {
            type = NULL_DISSOCIATION;
            nR = n;
            Unorm = m_Udis.at(n)/m_mcSq;
            FEnorm = m_FEdis.at(n)/m_mcSq;
            break;
         }
      }
   }
   if(nR<0) {
      for (int n=0; n<m_Qizn0.size(); n++) {
         Qsum += m_Qizn0.at(n);
         if(rand_num<=Qsum/a_sigma) {
            type = NULL_IONIZATION;
            nR = n;
            Unorm = m_Uizn.at(n)/m_mcSq;
            break;
         }
      } 
   }
   CH_assert(nR>=0);

   info.type  = type;
   info.nR    = nR;
   info.Unorm = Unorm;
   info.FEnorm = FEnorm;

   return info;

}

Real MonteCarloNull::getSigma( Real&  a_xi, const Real  a_beta ) const
{
   Real sigma_tot = 0.0;
   a_xi = 0.0;
   
   // compute the relative center-of-mass energy in eV
   const Real betaSq = a_beta*a_beta;
   const Real KE_cm = m_mu*m_mcSq*betaSq/2.0;

   //if(KE_cm>10.0) cout << "JRA KE_cm = " << KE_cm << endl;

   // find the index for interpolation
   int index = 0;
   const int N = m_Evec.size();
   if(KE_cm>=m_Evec.back()) {
      index = m_Evec.size();
      if(m_angular_scattering == OKHRIMOVSKYY) a_xi = m_xi.back();
   }
   else {

      index = N/2;
      while (KE_cm < m_Evec[index]) index--;
      while (KE_cm > m_Evec[index+1]) index++;
         
      // get xi
      if(m_angular_scattering == OKHRIMOVSKYY) {
         if(m_Evec[index]*KE_cm>0.0) {
            a_xi = ScatteringUtils::semilogInterp(m_Evec,m_xi,KE_cm,index);
         }
         else {
            a_xi = MathUtils::linearInterp(m_Evec,m_xi,KE_cm,index);
         }
      }

   }

   // get cross sections
   m_Qela0 = 0.0;
   std::fill(m_Qexc0.begin(), m_Qexc0.end(), 0.0 );
   std::fill(m_Qdis0.begin(), m_Qdis0.end(), 0.0 );
   std::fill(m_Qizn0.begin(), m_Qizn0.end(), 0.0 );

   m_Qela0 = getElasticSigma(index, KE_cm);
   sigma_tot = m_Qela0;
 
   for (int n=0; n<m_Qexc0.size(); n++) {
      m_Qexc0[n] = getExcitationSigma(index, n, KE_cm);
      sigma_tot += m_Qexc0[n];
   }

   for (int n=0; n<m_Qdis0.size(); n++) {
      m_Qdis0[n] = getDissociationSigma(index, n, KE_cm);
      sigma_tot += m_Qdis0[n];
   }

   for (int n=0; n<m_Qizn0.size(); n++) {
      m_Qizn0[n] = getIonizationSigma(index, n, KE_cm);
      sigma_tot += m_Qizn0[n];
   }
   
   return sigma_tot;
}

Real MonteCarloNull::getElasticSigma( const int   a_index, 
                                      const Real  a_KE ) const
{
 
   Real sigma=0.0;
   
   if(m_angular_scattering == ISOTROPIC) {

      // use momentum-transfer cross section for isotropic scattering
      if(a_KE >= m_Evec.back()) { // SigM ~ ln(E)/E
         sigma = m_Qelm.back();
         //sigma = m_Qelm.back()*log(a_KE)/log(m_Evec.back())*m_Evec.back()/a_KE;
      }
      else {
         if(m_use_loglog_interp && m_Qelm[a_index]*m_Evec[a_index]>0.0) {
            sigma = ScatteringUtils::loglogInterp(m_Evec,m_Qelm,a_KE,a_index);
         }
         else {
            sigma = MathUtils::linearInterp(m_Evec,m_Qelm,a_KE,a_index);
         }
      }
      
   }
   else if(m_angular_scattering == OKHRIMOVSKYY) {

      if(a_KE >= m_Evec.back()) { // SigT ~ 1/E
         sigma = m_Qela.back();
         //sigma = m_Qela.back()*m_Evec.back()/a_KE;
      }
      else {
         if(m_use_loglog_interp && m_Qela[a_index]*m_Evec[a_index]>0.0) {
            sigma = ScatteringUtils::loglogInterp(m_Evec,m_Qela,a_KE,a_index);
         }
         else {
            sigma = MathUtils::linearInterp(m_Evec,m_Qela,a_KE,a_index);
         }
      }

   }
   else {
      MayDay::Error("non-implemented scattering model in MonteCarloNull class");
   }
   
   return sigma;   

}  

Real MonteCarloNull::getExcitationSigma( const int   a_index, 
                                         const int   a_comp,
                                         const Real  a_KE ) const
{
 
   Real sigma=0.0;
   
   if(a_KE >= m_Evec.back()) {
      if(m_exc_type[a_comp]==FORBIDDEN) { // SigExc ~ 1/E^3
         sigma = m_Qexc[a_comp].back()*pow(m_Evec.back()/a_KE,3);
      }
      if(m_exc_type[a_comp]==ALLOWED) { // SigExc ~ log(E)/E
         sigma = m_Qexc[a_comp].back();
         //sigma = m_Qexc[a_comp].back()*log(a_KE)/log(m_Evec.back())*m_Evec.back()/a_KE;
      }
   }
   else if(a_KE > m_Uexc.at(a_comp)) {
      if(m_use_loglog_interp && m_Qexc[a_comp][a_index]*m_Evec[a_index]>0.0) {
         sigma = ScatteringUtils::loglogInterp(m_Evec,m_Qexc[a_comp],a_KE,a_index);
      }
      else {
         sigma = MathUtils::linearInterp(m_Evec,m_Qexc[a_comp],a_KE,a_index);
      }
   }

   return sigma;   

}  

Real MonteCarloNull::getDissociationSigma( const int   a_index, 
                                           const int   a_comp,
                                           const Real  a_KE ) const
{
 
   Real sigma=0.0;
   
   if(a_KE >= m_Evec.back()) {
      if(m_dis_type[a_comp]==DIS_EXCITATION || 
         m_dis_type[a_comp]==DIS_RECOMBINATION) { // forbidden SigExc ~ 1/E^3
         sigma = m_Qdis[a_comp].back()*pow(m_Evec.back()/a_KE,3);
      }
      if(m_dis_type[a_comp]==DIS_IONIZATION) { // SigIzn ~ log(E)/E
         sigma = m_Qdis[a_comp].back();
         //sigma = m_Qdis[a_comp].back()*log(a_KE)/log(m_Evec.back())*m_Evec.back()/a_KE;
      }
   }
   else if(a_KE > m_Udis.at(a_comp)) {
      if(m_use_loglog_interp && m_Qdis[a_comp][a_index]*m_Evec[a_index]>0.0) {
         sigma = ScatteringUtils::loglogInterp(m_Evec,m_Qdis[a_comp],a_KE,a_index);
      }
      else {
         sigma = MathUtils::linearInterp(m_Evec,m_Qdis[a_comp],a_KE,a_index);
      }
   }

   return sigma;   

}  

Real MonteCarloNull::getIonizationSigma( const int   a_index, 
                                         const int   a_comp,
                                         const Real  a_KE ) const
{
 
   Real sigma=0.0;
   
   if(a_KE >= m_Evec.back()) { // ionization, SigIzn ~ log(E)/E
      sigma = m_Qizn[a_comp].back();
      ///sigma = m_Qizn[a_comp].back()*log(a_KE)/log(m_Evec.back())*m_Evec.back()/a_KE;
   }
   else if(a_KE > m_Uizn.at(a_comp)) {
      if(m_use_loglog_interp && m_Qizn[a_comp][a_index]*m_Evec[a_index]>0.0) {
         sigma = ScatteringUtils::loglogInterp(m_Evec,m_Qizn[a_comp],a_KE,a_index);
      }
      else {
         sigma = MathUtils::linearInterp(m_Evec,m_Qizn[a_comp],a_KE,a_index);
      }
   }

   return sigma;   

}  

void MonteCarloNull::setOkhrimovskyyXI()
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
   if(!m_no_elastic) {
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
            MayDay::Error( "MonteCarloNull scattering xi calc: iter = iter_max" );
            break;
         }
      }
      m_xi.at(n) = x;

   }
   }

}

#include "NamespaceFooter.H"

