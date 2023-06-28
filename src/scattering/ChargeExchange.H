#ifndef _CHARGEEXCHANGE_H_
#define _CHARGEEXCHANGE_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * Charge exchange scattering: A+ + B => A + B+
 *
*/
class ChargeExchange : public Scattering
{
   public:

      ChargeExchange( ParmParse&    a_pp, 
                const WeightMethod  a_weight_method,
                const int           a_verbosity )
         : Scattering(a_verbosity),
           m_use_loglog_interp(false),
           m_weight_method(a_weight_method)
      {
         m_scatter_dt = DBL_MAX;
         parseParameters( a_pp );
      }

      virtual ~ChargeExchange() {;}
      
      virtual ScatteringType getScatteringType() const { return CHARGE_EXCHANGE; }
      
      virtual int species1() const { return m_sp1; }
      
      virtual int species2() const { return m_sp2; }
      
      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const; 
      
      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "ChargeExchange scattering parameters:" << std::endl;
            std::cout << "  species A name:          " << m_species1_name << std::endl;
            std::cout << "  species B name:          " << m_species2_name << std::endl;
            std::cout << "  weight method:         " << m_weight_method << std::endl;
            std::cout << "  cross section type:    " << m_xsec_type << std::endl;
            if(m_xsec_type==TEXT_FILE) {
               std::cout << "  file name:       " << m_xsec_fname << std::endl;
               std::cout << "  reaction:        " << m_reaction_label << std::endl;
               std::cout << "  loglog interp:   " << m_use_loglog_interp << std::endl;
            }
            std::cout << std::endl;
         }
      }
      
      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const;
      
   private:
      
      void ionImpact( PicSpecies&  a_picSpecies1,
                      PicSpecies&  a_picSpecies2,
                const DomainGrid&  a_mesh,
                const Real         a_dt_sec ) const;
      
      void setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                        const LevelData<FArrayBox>&  a_energyDensity1,
                        const LevelData<FArrayBox>&  a_numberDensity2,
                        const LevelData<FArrayBox>&  a_energyDensity2 ) const;
      
      void readCrossSectionFile();
      
      enum XSEC_TYPE {
         ANALYTIC,
         TEXT_FILE
      };
      
      enum ANALYTIC_MODEL {
         HpH,
         HepHe,
         H2pH2
      };

      // prohibit copying
      //ChargeExchange( const Constant& );
      //ChargeExchange& operator=( const Constant& );
      
      Real getSigmaIso( const Real  a_beta ) const;
      Real getSigmaBack( const Real  a_beta ) const;
      Real getTextSigmaIso( const Real  a_beta ) const;
      Real getTextSigmaBack( const Real  a_beta ) const;
      void setgSigmaMax();

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp ) 
      {
        a_pp.get( "speciesA", m_sp1 ); 
        a_pp.get( "speciesB", m_sp2 );

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
             MayDay::Error( "ChargeExchange: invalid weight_method" );
          }
        }

        // get the cross section type
        std::string xsec_type;
        a_pp.get( "cross_section_type", xsec_type );
        if(xsec_type=="analytic") {
          m_xsec_type = ANALYTIC;
          //a_pp.get( "ANALYTIC.model", m_analytic_model ); // [m^2]          
        }
        else if(xsec_type=="text file") {
          m_xsec_type = TEXT_FILE;
          a_pp.get("file_name",m_xsec_fname);
          a_pp.query("use_loglog_interpolation",m_use_loglog_interp);
        }
        else {
          MayDay::Error( "ChargeExchange: invalid cross_section_type" );
        }
      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      Real m_mass1;  // species 1 mass / electron mass
      Real m_mass2;  // species 1 mass / electron mass
      Real m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
      
      Real m_gSigmaMax; // max(vrel*QTot(vrel)) 
      XSEC_TYPE m_xsec_type;
      bool m_use_loglog_interp;
      ANALYTIC_MODEL m_xsec_model;
      WeightMethod m_weight_method;

      // cross section file members
      std::string m_xsec_fname;
      std::string m_reaction_label;
      std::vector<double> m_Evec;
      std::vector<double> m_Qiso;
      std::vector<double> m_Qback;
      
      Real m_mcSq;

};

#include "NamespaceFooter.H"

#endif

