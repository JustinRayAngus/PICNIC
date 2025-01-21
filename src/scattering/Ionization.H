#ifndef _IONIZATION_H_
#define _IONIZATION_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 *
 * generic ionization class: A + B ==> A + e + B+
 *
*/

class Ionization : public Scattering
{
   public:
      
      Ionization( ParmParse&  a_pp, 
            const int&        a_verbosity )
         : Scattering(a_verbosity),
           m_use_loglog_interp(false),
           m_angular_scattering(ISOTROPIC),
           m_energy_sharing(ZERO)
      {
         m_scatter_dt = DBL_MAX;
         m_deltaE_izn = 0.0;
         parseParameters( a_pp );
      }

      virtual ~Ionization() {;}
      
      virtual ScatteringType getScatteringType() const { return IONIZATION; }
      
      virtual int species1() const { return m_sp1; }
      
      virtual int species2() const { return m_sp2; }
      
      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const;
      
      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }
      
      virtual Real getDeltaEizn() const { return m_deltaE_izn; }
      virtual void zeroDeltaEizn() const { m_deltaE_izn=0.0; } 

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Ionization scattering parameters:" << std::endl;
            std::cout << "  species A name:        " << m_species1_name << std::endl;
            std::cout << "  species B name:        " << m_species2_name << std::endl;
            std::cout << "  ele species name:      " << m_speciesE_name << std::endl;
            std::cout << "  ion species name:      " << m_speciesI_name << std::endl;
            std::cout << "  Izn potential [eV]:    " << m_Uizn << std::endl;
            std::cout << "  angular scattering:    " << m_angular_scattering << std::endl;
            std::cout << "  energy sharing:        " << m_energy_sharing << std::endl;
            std::cout << "  cross section type:    " << m_xsec_type << std::endl;
            if(m_xsec_type==LOTZ) {
               std::cout << "  lotz a [m^2]:    " << m_lotz_a << std::endl;
               std::cout << "  lotz xi:         " << m_lotz_xi << std::endl;
            }
            if(m_xsec_type==TEXT_FILE) {
               std::cout << "  file name:       " << m_xsec_fname << std::endl;
               std::cout << "  reaction:        " << m_reaction_label << std::endl;
               std::cout << "  loglog interp:   " << m_use_loglog_interp << std::endl;
            }
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

      void electronImpact( PicChargedSpecies&  a_picSpecies1,
                           PicChargedSpecies&  a_picSpecies2,
                           PicChargedSpecies&  a_picSpeciesE,
                           PicChargedSpecies&  a_picSpeciesI,
                     const DomainGrid&  a_mesh,
                     const Real         a_dt_sec ) const;
      
      enum XSEC_TYPE {
         LOTZ,
         TEXT_FILE
      };
      
      enum ENERGY_SHARING {
         ZERO,
         EQUAL,
         OPAL
      };
      
      enum ANGULAR_SCATTERING {
         ISOTROPIC,
         SURENDRA,
         OKHRIMOVSKYY
      };

      // prohibit copying
      //Ionization( const Constant& );
      //Ionization& operator=( const Constant& );

      Real getSigma( const Real  a_beta ) const;
      Real getLotzSigma( const Real  a_beta ) const;
      Real getTextSigma( const Real  a_beta ) const;

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp ) 
      {
        a_pp.get( "speciesA", m_sp1 ); 
        a_pp.get( "speciesB", m_sp2 );
        a_pp.get( "speciesE", m_spE );
        a_pp.get( "speciesI", m_spI );
        a_pp.get( "potential", m_Uizn );
        
        // look for specified angular scattering model
        std::string angular_scattering;
        a_pp.query("angular_scattering",angular_scattering);
        if(!angular_scattering.empty()) {
          if(angular_scattering=="ISOTROPIC" || angular_scattering=="isotropic") {
            m_angular_scattering = ISOTROPIC;
          }
          else if(angular_scattering=="SURENDRA" || angular_scattering=="surendra") {
            m_angular_scattering = SURENDRA;
          }
          else if(angular_scattering=="OKHRIMOVSKYY" || angular_scattering=="okhrimovskyy") {
            m_angular_scattering = OKHRIMOVSKYY;
          }
          else {
             cout << "angular_scattering = " << angular_scattering << endl;
             MayDay::Error( "Ionization: invalid angular_scattering" );
          }
        }

        // look for specified secondary energy sharing model
        std::string energy_sharing;
        a_pp.query("energy_sharing",energy_sharing);
        if(!energy_sharing.empty()) {
          if(energy_sharing=="ZERO" || energy_sharing=="zero") {
            m_energy_sharing = ZERO;
          }
          else if(energy_sharing=="EQUAL" || energy_sharing=="equal") {
            m_energy_sharing = EQUAL;
          }
          else if(energy_sharing=="OPAL" || energy_sharing=="opal") {
            m_energy_sharing = OPAL;
            a_pp.get( "opal.bi", m_opal_bi );
            a_pp.get( "opal.ki", m_opal_ki );
          }
          else {
             cout << "energy_sharing = " << energy_sharing << endl;
             MayDay::Error( "Ionization: invalid energy_sharing" );
          }
        }

        // get the cross section type
        std::string xsec_type;
        a_pp.get( "cross_section_type", xsec_type );
        if(xsec_type=="lotz" || xsec_type=="Lotz") {
          m_xsec_type = LOTZ;
          a_pp.get( "Lotz.a", m_lotz_a );   // fitting coeff a [m^2]
          a_pp.get( "Lotz.b", m_lotz_b );   // fitting coeff b
          a_pp.get( "Lotz.c", m_lotz_c );   // fitting coeff c
          a_pp.get( "Lotz.xi", m_lotz_xi ); // # ele in outer shell          
        }
        else if(xsec_type=="text file") {
          m_xsec_type = TEXT_FILE;
          a_pp.get("file_name",m_xsec_fname);
          a_pp.query("use_loglog_interpolation",m_use_loglog_interp);
        }
        else {
          MayDay::Error( "Ionization: invalid cross_section_type" );
        }

      }

      std::string m_species1_name;
      std::string m_species2_name;
      std::string m_speciesE_name;
      std::string m_speciesI_name;
      mutable Real m_scatter_dt; // (mean free time)
      mutable Real m_deltaE_izn; // running sum of energy spent ionizing
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      int m_spE;     // electron species identifier
      int m_spI;     // ion species identifier
      Real m_mass1;  // species 1 mass / electron mass
      Real m_mass2;  // species 1 mass / electron mass
      Real m_massE;  // species E mass / electron mass
      Real m_massI;  // species I mass / electron mass
      Real m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
      Real m_Uizn;   // ionization potential [eV]
      //
      XSEC_TYPE m_xsec_type;
      bool m_use_loglog_interp;
      ANGULAR_SCATTERING m_angular_scattering;
      ENERGY_SHARING m_energy_sharing;
      //
      Real m_opal_bi, m_opal_ki;
      Real m_lotz_a, m_lotz_b, m_lotz_c;
      int m_lotz_xi;
      //
      // cross section file members
      std::string m_xsec_fname;
      std::string m_reaction_label;
      std::vector<double> m_Evec;   
      std::vector<double> m_Qvec;   
   

      Real m_mcSq;

};

#include "NamespaceFooter.H"

#endif

