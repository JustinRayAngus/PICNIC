#ifndef _ELASTIC_H_
#define _ELASTIC_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 *
 * generic elastic scattering class: A + B ==> A + B
 *
*/

class Elastic : public Scattering
{
   public:
      
      Elastic( ParmParse&    a_pp, 
         const WeightMethod  a_weight_method,
         const int           a_verbosity )
         : Scattering(a_verbosity),
           m_print_cross_sections(false),
           m_use_loglog_interp(false),
           m_angular_scattering(ISOTROPIC),
           m_weight_method(a_weight_method)
      {
         m_scatter_dt = DBL_MAX;
         parseParameters( a_pp );
      }

      virtual ~Elastic() {;}
      
      virtual ScatteringType getScatteringType() const { return ELASTIC; }
      
      virtual int species1() const { return m_sp1; }
      
      virtual int species2() const { return m_sp2; }
      
      virtual void setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity1,
                                    const LevelData<FArrayBox>&  a_energyDensity1,
                                    const LevelData<FArrayBox>&  a_numberDensity2,
                                    const LevelData<FArrayBox>&  a_energyDensity2 ) const;
    
      virtual void initialize( const PicSpeciesPtrVect&  a_picSpeciesPtrVect,
                               const DomainGrid&         a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }
      
      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Elastic scattering parameters:" << std::endl;
            std::cout << "  species A name:        " << m_species1_name << std::endl;
            std::cout << "  species B name:        " << m_species2_name << std::endl;
            std::cout << "  angular scattering:    " << m_angular_scattering << std::endl;
            std::cout << "  weight method:         " << m_weight_method << std::endl;
            std::cout << "  cross section type:    " << m_xsec_type << std::endl;
            if(m_xsec_type==CONSTANT) {
               std::cout << "  sigma [m^2] = :         " << m_const_sigma << std::endl;
            }
            if(m_xsec_type==TEXT_FILE) {
               std::cout << "  file name:       " << m_xsec_fname << std::endl;
               std::cout << "  reaction:        " << m_reaction_label << std::endl;
               std::cout << "  loglog interp:   " << m_use_loglog_interp << std::endl;
               std::cout << std::endl;
               if(m_print_cross_sections) printCrossSections();
            }
         }
      }
      
      virtual void applyScattering( PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                              const DomainGrid&         a_mesh,
                              const Real                a_dt_sec ) const;
      
   private:

      void readCrossSectionFile();
      
      void setOkhrimovskyyXI();
      
      void printCrossSections() const;

      void electronImpact( PicSpecies&  a_picSpecies1,
                           PicSpecies&  a_picSpecies2,
                     const DomainGrid&  a_mesh,
                     const Real         a_dt_sec ) const;
      
      enum XSEC_TYPE {
         CONSTANT,
         TEXT_FILE
      };
      
      enum ANGULAR_SCATTERING {
         ISOTROPIC,
         SURENDRA,
         OKHRIMOVSKYY
      };
      
      Real getSigma( Real&  a_xi, const Real  a_beta ) const;
      Real getTextSigma( Real&  a_xi, const Real  a_beta ) const;

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp ) 
      {
        a_pp.get( "speciesA", m_sp1 ); 
        a_pp.get( "speciesB", m_sp2 );
        
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
             MayDay::Error( "Elastic: invalid angular_scattering" );
          }
        }
        
        // look for specified method for A and B having different weights
        std::string weight_method;
        a_pp.query("weight_method",weight_method);
        if(!weight_method.empty()) {
          if(weight_method=="PROBABILISTIC" || weight_method=="probabilistic") {
            m_weight_method = PROBABILISTIC;
          }
          else if(weight_method=="CONSERVATIVE" || weight_method=="conservative") {
            m_weight_method = CONSERVATIVE;
          }
          else {
             cout << "weight_method = " << weight_method << endl;
             MayDay::Error( "Elastic: invalid weight_method" );
          }
        }

        // get the cross section type
        std::string xsec_type;
        a_pp.get( "cross_section_type", xsec_type );
        if(xsec_type=="constant") {
          m_xsec_type = CONSTANT;
          a_pp.get( "constant.sigma", m_const_sigma ); // [m^2]          
        }
        else if(xsec_type=="text file") {
          m_xsec_type = TEXT_FILE;
          a_pp.get("file_name",m_xsec_fname);
          a_pp.query("print_cross_sections",m_print_cross_sections);
          a_pp.query("use_loglog_interpolation",m_use_loglog_interp);
        }
        else {
          MayDay::Error( "Elastic: invalid cross_section_type" );
        }

      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      Real m_mass1;  // species 1 mass / electron mass
      Real m_mass2;  // species 1 mass / electron mass
      Real m_m1om2;  // mass ratio
      Real m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
      //
      XSEC_TYPE m_xsec_type;
      bool m_print_cross_sections;
      bool m_use_loglog_interp;
      ANGULAR_SCATTERING m_angular_scattering;
      WeightMethod m_weight_method;
      //
      Real m_const_sigma;          
      
      // cross section file members
      std::string m_xsec_fname;
      std::string m_reaction_label;
      std::vector<double> m_Evec;   
      std::vector<double> m_Qela; // total cross section   
      std::vector<double> m_Qelm; // momentum-transfer cross section
      std::vector<double> m_xi;   // Okhrimovskyy xi

      Real m_mcSq;

};

#include "NamespaceFooter.H"

#endif

