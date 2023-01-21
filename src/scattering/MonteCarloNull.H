#ifndef _MONTECARLONULL_H_
#define _MOTECARLONULL_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 *
 * monte carlo null scattering class: e + B ==> e + B
 * with elastic, excitation, and ionization collisions.
 *
 * B can have multiple internal states
 *
*/

class MonteCarloNull : public Scattering
{
   public:
      
      MonteCarloNull( ParmParse&    a_pp, 
                const WeightMethod  a_weight_method,
                const int           a_verbosity )
         : Scattering(a_verbosity),
           m_use_loglog_interp(false),
           m_angular_scattering(ISOTROPIC),
           m_weight_method(a_weight_method),
           m_no_elastic(false),
           m_print_cross_sections(false)
      {
         m_scatter_dt = DBL_MAX;
         m_deltaE_izn = 0.0;
         m_deltaE_exc = 0.0;
         parseParameters( a_pp );
      }

      virtual ~MonteCarloNull() {;}
      
      virtual ScatteringType getScatteringType() const { return MONTE_CARLO_NULL; }
      
      virtual int species1() const { return m_sp1; }
      
      virtual int species2() const { return m_sp2; }
      
      virtual void setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity1,
                                    const LevelData<FArrayBox>&  a_energyDensity1,
                                    const LevelData<FArrayBox>&  a_numberDensity2,
                                    const LevelData<FArrayBox>&  a_energyDensity2 ) const;
    
      virtual void initialize( const PicSpeciesPtrVect&  a_picSpeciesPtrVect,
                               const DomainGrid&         a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }

      virtual Real getDeltaEizn() const { return m_deltaE_izn; }
      virtual void zeroDeltaEizn() const { m_deltaE_izn=0.0; }      
      
      virtual Real getDeltaEexc() const { return m_deltaE_exc; }
      virtual void zeroDeltaEexc() const { m_deltaE_exc=0.0; }      

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "MonteCarloNull scattering parameters:" << std::endl;
            std::cout << "  file name:       " << m_xsec_fname << std::endl;
            std::cout << "  reaction:        " << m_reaction_label << std::endl;
            std::cout << "  species A name:       " << m_species1_name << std::endl;
            std::cout << "  species B name:       " << m_species2_name << std::endl;
            std::cout << "  num excitations    =   " << m_num_exc << std::endl;
            std::cout << "  exc potential [eV] =   ";
            for (int n=0; n<m_num_exc; n++) cout << m_Uexc[n] << "  ";
            cout << std::endl;
            std::cout << "  num dissociations  =   " << m_num_dis << std::endl;
            std::cout << "  dis potential [eV]   =   ";
            for (int n=0; n<m_num_dis; n++) cout << m_Udis[n] << "  ";
            cout << std::endl;
            std::cout << "  dis-frag energy [eV] =   ";
            for (int n=0; n<m_num_dis; n++) cout << m_FEdis[n] << "  ";
            cout << std::endl;
            std::cout << "  num ionizations    =   " << m_num_izn << std::endl;
            std::cout << "  izn potential [eV] =   ";
            for (int n=0; n<m_num_izn; n++) cout << m_Uizn[n] << "  ";
            cout << std::endl;
            std::cout << "  angular scattering:   " << m_angular_scattering << std::endl;
            std::cout << "  weight method:        " << m_weight_method << std::endl;
            std::cout << "  loglog interp:        " << m_use_loglog_interp << std::endl;
            std::cout << std::endl;
            if(m_print_cross_sections) printCrossSections();
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
                           std::vector<PicSpeciesPtr>&  a_exc_species, 
                           std::vector<PicSpeciesPtr>&  a_dis_species1, 
                           std::vector<PicSpeciesPtr>&  a_dis_species2, 
                           std::vector<PicSpeciesPtr>&  a_dis_species3, 
                           std::vector<PicSpeciesPtr>&  a_izn_speciesE, 
                           std::vector<PicSpeciesPtr>&  a_izn_speciesI, 
                     const DomainGrid&  a_mesh,
                     const Real         a_dt_sec ) const;
      
      enum ANGULAR_SCATTERING {
         ISOTROPIC,
         //SURENDRA,
         OKHRIMOVSKYY
      };
      
      enum EXCITATION_TYPE {
         FORBIDDEN,
         ALLOWED
      };
      
      enum DISSOCIATION_TYPE {
         DIS_EXCITATION,
         DIS_IONIZATION,
         DIS_RECOMBINATION
      };
      
      enum COLLISION_TYPE {
         NULL_ELASTIC,
         NULL_EXCITATION,
         NULL_DISSOCIATION,
         NULL_IONIZATION,
      };
      
      struct collisionInfo {
         COLLISION_TYPE type;
         int nR;     // sub-reaction number specifier
         Real Unorm;  // normalized reaction potential energy
         Real FEnorm; // normalized reaction fragmentation energy (only for dissociation)
      };

      collisionInfo getCollisionInfo( const Real  a_sigma ) const;

      Real getSigma( Real&  a_xi, const Real  a_beta ) const;
      Real getElasticSigma( const int  a_index, const Real  a_KE ) const;
      Real getExcitationSigma( const int   a_index, 
                               const int   a_comp,
                               const Real  a_KE ) const;
      Real getDissociationSigma( const int   a_index, 
                                 const int   a_comp,
                                 const Real  a_KE ) const;
      Real getIonizationSigma( const int   a_index, 
                               const int   a_comp,
                               const Real  a_KE ) const;

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
          else if(angular_scattering=="OKHRIMOVSKYY" || angular_scattering=="okhrimovskyy") {
            m_angular_scattering = OKHRIMOVSKYY;
          }
          else {
             cout << "angular_scattering = " << angular_scattering << endl;
             MayDay::Error( "MonteCarloNull: invalid angular_scattering" );
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
             MayDay::Error( "MonteCarloNull: invalid weight_method" );
          }
        }
        
        a_pp.query("use_loglog_interpolation",m_use_loglog_interp);

        a_pp.get("file_name",m_xsec_fname);
        
        a_pp.query("no_elastic",m_no_elastic);
        a_pp.query("print_cross_sections",m_print_cross_sections);
        
        a_pp.get("num_exc",m_num_exc); // must match that in cross section file
        m_exc_species.resize(m_num_exc,-1);
        if(m_num_exc>0 && a_pp.contains("exc_species")) {
           a_pp.getarr("exc_species",m_exc_species,0,m_num_exc);
        }

        a_pp.get("num_dis",m_num_dis); // must match that in cross section file
        m_dis_species1.resize(m_num_dis,-1);
        m_dis_species2.resize(m_num_dis,-1);
        m_dis_species3.resize(m_num_dis,-1);
        if(m_num_dis>0) {
           a_pp.getarr("dis_species1",m_dis_species1,0,m_num_dis);
           a_pp.getarr("dis_species2",m_dis_species2,0,m_num_dis);
           a_pp.getarr("dis_species3",m_dis_species3,0,m_num_dis);
        }

        a_pp.get("num_izn",m_num_izn); // must match that in cross section file
        m_izn_speciesE.resize(m_num_izn,-1);
        m_izn_speciesI.resize(m_num_izn,-1);
        if(m_num_izn>0 && a_pp.contains("izn_speciesI")) {
           a_pp.getarr("izn_speciesE",m_izn_speciesE,0,m_num_izn);
           a_pp.getarr("izn_speciesI",m_izn_speciesI,0,m_num_izn);
        }

      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // mean free time
      mutable Real m_deltaE_izn; // running sum of energy spent ionizing
      mutable Real m_deltaE_exc; // running sum of energy spent exciting
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      Real m_mass1;  // species 1 mass / electron mass
      Real m_mass2;  // species 1 mass / electron mass
      Real m_m1om2;  // mass ratio
      Real m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
      //
      bool m_use_loglog_interp;
      ANGULAR_SCATTERING m_angular_scattering;
      WeightMethod m_weight_method;
      
      // cross section file members
      std::string m_xsec_fname;
      std::string m_reaction_label;
      std::vector<double> m_Evec;   
      std::vector<double> m_xi;    // Okhrimovskyy xi
      std::vector<double> m_Qela;  // elastic cross section
      std::vector<double> m_Qelm;  // elastic-momentum cross section
      std::vector<vector<double>> m_Qexc;  // excitation cross sections
      std::vector<vector<double>> m_Qdis;  // dissociation cross sections
      std::vector<vector<double>> m_Qizn;  // ionization cross sections
     
      std::vector<EXCITATION_TYPE> m_exc_type;
      std::vector<DISSOCIATION_TYPE> m_dis_type;
      std::vector<double> m_Uexc;   
      std::vector<double> m_Udis;   
      std::vector<double> m_FEdis;   
      std::vector<double> m_Uizn;

      bool m_no_elastic;
      bool m_print_cross_sections;
      int m_num_exc;
      int m_num_dis;
      int m_num_izn;
      std::vector<int> m_exc_species;
      std::vector<int> m_dis_species1;
      std::vector<int> m_dis_species2;
      std::vector<int> m_dis_species3;
      std::vector<int> m_izn_speciesI;
      std::vector<int> m_izn_speciesE;

      mutable double m_Qela0;
      mutable std::vector<double> m_Qexc0;
      mutable std::vector<double> m_Qdis0;
      mutable std::vector<double> m_Qizn0;

      Real m_mcSq;

};

#include "NamespaceFooter.H"

#endif

