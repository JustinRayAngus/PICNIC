#ifndef _FUSION_H_
#define _FUSION_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * Fusion scattering class:
 * D + D   ==> T(1.01 MeV)   + p(3.02 MeV)
 * D + D   ==> He3(0.82 MeV) + n(2.45 MeV)
 * D + T   ==> He4(3.5 MeV)  + n(14.1 MeV)
 * D + He3 ==> He4(3.6 MeV)  + p(14.7 MeV)
*/

class Fusion : public Scattering
{
   public:

      Fusion( ParmParse&  a_pp,
        const int&        a_verbosity )
         : Scattering(a_verbosity)
      {
         m_scatter_dt = DBL_MAX;
         m_deltaE_fusion = 0.0;
         parseParameters( a_pp );
      }

      virtual ~Fusion() = default;

      virtual ScatteringType getScatteringType() const { return FUSION; }

      virtual std::string getScatteringSubTypeName() const {
          if      (m_fusion_type==DDa)  { return "DDa"; }
          else if (m_fusion_type==DDb)  { return "DDb"; }
          else if (m_fusion_type==DDab) { return "DDab"; }
          else if (m_fusion_type==DT)   { return "DT"; }
          else if (m_fusion_type==DHe3) { return "DHe3"; }
          else { return "invalid"; }
      }

      virtual int species1() const { return m_sp1; }

      virtual int species2() const { return m_sp2; }

      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const;

      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }

      virtual Real getFusionQ() const { return m_Q; }
      virtual Real getFusionQb() const { return m_Qb; }
      virtual Real getDeltaEfusion() const { return m_deltaE_fusion; }
      virtual void zeroDeltaEfusion() const { m_deltaE_fusion=0.0; }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Fusion scattering parameters:" << std::endl;
            std::cout << "  species 1 name:   " << m_species1_name << std::endl;
            std::cout << "  species 2 name:   " << m_species2_name << std::endl;
            std::cout << "  species 3 name:   " << m_species3_name << std::endl;
            std::cout << "  species 4 name:   " << m_species4_name << std::endl;
            if (m_fusion_type==DDab) {
                std::cout << "  species 3b name:  " << m_species3b_name << std::endl;
                std::cout << "  species 4b name:  " << m_species4b_name << std::endl;
                std::cout << "  fusion type:  D + D (a) ==>   T(1.01 MeV) + p(3.02 MeV)" << std::endl;
                std::cout << "                      (b) ==> He3(0.82 MeV) + n(2.45 MeV)" << std::endl;
            }
            else if (m_fusion_type==DDa) {
                std::cout << "  fusion type:  D + D ==>   T(1.01 MeV) + p(3.02 MeV)" << std::endl;
            }
            else if (m_fusion_type==DDb) {
                std::cout << "  fusion type:  D + D ==> He3(0.82 MeV) + n(2.45 MeV)" << std::endl;
            }
            else if (m_fusion_type==DT) {
                std::cout << "  fusion type:  D + T ==> He4(3.5 MeV) + n(14.1 MeV)" << std::endl;
            }
            else if (m_fusion_type==DHe3) {
                std::cout << "  fusion type:  D + He3 ==> He4(3.6 MeV) + p(14.7 MeV)" << std::endl;
            }
            if (m_fusion_type==DDab) { 
                std::cout << "  fusion Qa = " << m_Q  << " eV" << std::endl;
                std::cout << "  fusion Qb = " << m_Qb << " eV" << std::endl;
            }
            else { std::cout << "  fusion Q = " << m_Q << " eV" << std::endl; }
            std::cout << "  fmulti = " << m_fmulti << std::endl;
            std::cout << std::endl;
         }
      }

      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const;

   private:

      void setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                        const LevelData<FArrayBox>&  a_energyDensity1,
                        const LevelData<FArrayBox>&  a_numberDensity2,
                        const LevelData<FArrayBox>&  a_energyDensity2 ) const;

      void readCrossSectionFile();

      void interSpeciesFusion( PicSpecies&  a_picSpecies1,
                               PicSpecies&  a_picSpecies2,
                               PicSpecies&  a_picSpecies3,
                               PicSpecies&  a_picSpecies4,
                         const DomainGrid&  a_mesh,
                         const Real         a_dt_sec ) const;

      void intraSpeciesFusion( PicSpecies&  a_picSpecies1,
                               PicSpecies&  a_picSpecies3,
                               PicSpecies&  a_picSpecies4,
                               PicSpecies&  a_picSpecies3b,
                               PicSpecies&  a_picSpecies4b,
                         const DomainGrid&  a_mesh,
                         const Real         a_dt_sec ) const;

      enum XSEC_TYPE {
         LOTZ,
         TEXT_FILE
      };

      enum FUSION_TYPE {
         DDa,  // D + D   ==> T   + p
         DDb,  // D + D   ==> He3 + n
         DDab, // D + D   ==> T  + p and ==> He3 + n
         DT,   // D + T   ==> He4 + n
         DHe3  // D + He3 ==> He4 + p
      };

      Real getSigma( Real& a_ratio_b, const Real  a_Erel_eV ) const;
      Real getDDaSigma( const Real  a_Erel_eV ) const;
      Real getDDbSigma( const Real  a_Erel_eV ) const;
      Real getDTSigma( const Real  a_Erel_eV ) const;
      Real getDHe3Sigma( const Real  a_Erel_eV ) const;

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp )
      {
          // get the fusion type
          std::string fusion_type;
          a_pp.get( "fusion_type", fusion_type );
          if      (fusion_type=="DDa")  { m_fusion_type = DDa; }
          else if (fusion_type=="DDb")  { m_fusion_type = DDb; }
          else if (fusion_type=="DDab") { m_fusion_type = DDab; }
          else if (fusion_type=="DT")   { m_fusion_type = DT; }
          else if (fusion_type=="DHe3") { m_fusion_type = DHe3; }
          else {
              MayDay::Error( "Fusion: invalid fusion_type. Valid types are DDa, DDb, DDab, DT, and DHe3" );
          }
          a_pp.query( "fmulti", m_fmulti );

          // get number associated with each species
          a_pp.get( "species1", m_sp1 );
          a_pp.get( "species2", m_sp2 );

          if (m_fusion_type==DDab) {
              std::vector<int>  species3_vec, species4_vec;
              species3_vec.resize(2);
              species4_vec.resize(2);
              a_pp.getarr( "species3", species3_vec, 0, 2 );
              a_pp.getarr( "species4", species4_vec, 0, 2 );
              m_sp3 = species3_vec.at(0);
              m_sp4 = species4_vec.at(0);
              m_sp3b = species3_vec.at(1);
              m_sp4b = species4_vec.at(1);
          }
          else {
              a_pp.get( "species3", m_sp3 );
              a_pp.get( "species4", m_sp4 );
              m_sp3b = m_sp3;
              m_sp4b = m_sp4;
          }

      }

      std::string m_species1_name;
      std::string m_species2_name;
      std::string m_species3_name;
      std::string m_species4_name;
      std::string m_species3b_name;
      std::string m_species4b_name;
      mutable Real m_scatter_dt;    // (mean free time)
      mutable Real m_deltaE_fusion; // running sum of energy gained from fusion
      int m_sp1;      // species 1 identifier
      int m_sp2;      // species 2 identifier
      int m_sp3;      // species 3 identifier
      int m_sp4;      // species 4 identifier
      int m_sp3b;     // species 3b identifier
      int m_sp4b;     // species 4b identifier
      Real m_mass1;   // species 1 mass / electron mass
      Real m_mass2;   // species 2 mass / electron mass
      Real m_mass3;   // species 3 mass / electron mass
      Real m_mass4;   // species 4 mass / electron mass
      Real m_mass3b;  // species 3b mass / electron mass
      Real m_mass4b;  // species 4b mass / electron mass
      Real m_mu;      // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
      Real m_Q  = 0.0;     // Reaction energy [eV]
      Real m_Qb = 0.0;     // Reaction energy for channel b of DDab [eV]
      Real m_fmulti = 1.0; // production factor (for better statistics)

      FUSION_TYPE m_fusion_type;
      Real m_mcSq;

};

#include "NamespaceFooter.H"

#endif

