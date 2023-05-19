#ifndef _COULOMB_H_
#define _COULOMB_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * Monte-Carlo scattering model for Coulomb collisions in a plasma
 * See Takizuka and Abe, JCP, 25 (1977)
 * See Nanbu, PRE, 55 (1997)
 * See Nanbu and Yonemura, JCP, 145 (1998)
 * See Sentoku and Kemp, JCP, 227 (2008)
 * See Perez et al, PoP, 19 (2012)
 * See Higginson et al, JCP 413 (2020)
*/
class Coulomb : public Scattering
{
   public:

      Coulomb( ParmParse&  a_pp, 
               const CodeUnits&    a_units,
               const int&        a_verbosity )
         : Scattering(a_verbosity),
           m_angular_scattering(TAKIZUKA)
      {
         m_scatter_dt = DBL_MAX;
         m_Clog = 10.0;
         m_b90_codeToPhys = Constants::QE*Constants::QE/(2.0*Constants::PI*Constants::EP0*Constants::ME);
         m_volume_scale = a_units.getScale(a_units.VOLUME);
         parseParameters( a_pp );
      }
      
      Coulomb( const int   a_sp1,
               const int   a_sp2,
               const Real  a_Clog, 
               const Real  a_volume_scale,
               const std::string&  a_scatter_type, 
               const int   a_verbosity )
         : Scattering(a_verbosity),
           m_angular_scattering(TAKIZUKA)
      {
         m_scatter_dt = DBL_MAX;
         m_b90_codeToPhys = Constants::QE*Constants::QE/(2.0*Constants::PI*Constants::EP0*Constants::ME); 
         m_volume_scale = a_volume_scale;
         m_sp1 = a_sp1;
         m_sp2 = a_sp2;
         m_Clog = a_Clog; CH_assert(m_Clog>=3.0);
         if(a_scatter_type=="TAKIZUKA" || a_scatter_type=="takizuka") {
            m_angular_scattering = TAKIZUKA;
         }
         else if(a_scatter_type=="NANBU" || a_scatter_type=="nanbu") {
            m_angular_scattering = NANBU;
         }
         else {
            cout << "angular_scattering = " << a_scatter_type << endl;
            cout << " valid options: 1) takizuka, 2) nanbu " << endl;
            MayDay::Error( "Coulomb: invalid angular_scattering" );
         }
      }

      virtual ~Coulomb() {;}
      
      virtual ScatteringType getScatteringType() const { return COULOMB; }
      
      virtual int species1() const { return m_sp1; }
      
      virtual int species2() const { return m_sp2; }
      
      virtual void setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity,
                                    const LevelData<FArrayBox>&  a_energyDensity ) const;

      virtual void setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity1,
                                    const LevelData<FArrayBox>&  a_energyDensity1,
                                    const LevelData<FArrayBox>&  a_numberDensity2,
                                    const LevelData<FArrayBox>&  a_energyDensity2 ) const;
    
      virtual void initialize( const PicSpeciesPtrVect&  a_picSpeciesPtrVect,
                               const DomainGrid&         a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Coulomb scattering parameters:" << std::endl;
            std::cout << " species A name:         " << m_species1_name << std::endl;
            std::cout << " species B name:         " << m_species2_name << std::endl;
            if(m_angular_scattering==TAKIZUKA) {
            std::cout << " angular scattering:     TAKIZUKA " << std::endl;
            }
            if(m_angular_scattering==NANBU) {
            std::cout << " angular scattering:     NANBU " << std::endl;
            }
            std::cout << " Clog =                  " << m_Clog << std::endl << std::endl;
         }
      }
      
      virtual void applyScattering( PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                              const DomainGrid&         a_mesh,
                              const Real                a_dt_sec ) const;
      
   private:

      enum ANGULAR_SCATTERING {
         TAKIZUKA,
         NANBU 
      };

      // scatter like species
      virtual void applySelfScattering( PicSpecies&  a_picSpecies, 
                                  const DomainGrid&  a_mesh,
                                  const Real         a_dt_sec ) const;
      
      // scatter unlike species
      virtual void applyInterScattering( PicSpecies&  a_picSpecies1,
                                         PicSpecies&  a_picSpecies2,
                                   const DomainGrid&  a_mesh,
                                   const Real         a_dt_sec ) const;
 
      void GalileanScatter( std::array<Real,3>&  a_deltaU,
                      const std::array<Real,3>&  a_vp1,
                      const std::array<Real,3>&  a_vp2,
                      const Real                 a_den12,
                      const Real                 a_Clog,
                      const Real                 a_dt_sec ) const;
     
      void LorentzScatter( std::array<Real,3>&  a_up1,
                           std::array<Real,3>&  a_up2,
                     const bool                 a_scatter2,
                     const long double          a_mass1,
                     const long double          a_mass2,
                     const Real                 a_den12,
                     const Real                 a_dt_sec ) const;

      // prohibit copying
      //Coulomb( const Constant& );
      //Coulomb& operator=( const Constant& );

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp ) 
      {
         a_pp.get( "speciesA", m_sp1 ); 
         a_pp.get( "speciesB", m_sp2 );
         a_pp.query( "coulomb_logarithm", m_Clog ); CH_assert(m_Clog>=3.0);
        
         // look for specified scattering model
         std::string scatter_type_str;
         a_pp.query("angular_scattering",scatter_type_str);
         if(!scatter_type_str.empty()) {
            if(scatter_type_str=="TAKIZUKA" || scatter_type_str=="takizuka") {
               m_angular_scattering = TAKIZUKA;
            }
            else if(scatter_type_str=="NANBU" || scatter_type_str=="nanbu") {
               m_angular_scattering = NANBU;
            }
            else {
               cout << "angular_scattering = " << scatter_type_str << endl;
               cout << " valid options: 1) takizuka, 2) nanbu " << endl;
               MayDay::Error( "Coulomb: invalid angular_scattering" );
            }
         }

      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      mutable Real m_Clog;       // Coulomb Logarithm
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      signed int m_charge1; // species 1 charge / |electron charge|
      signed int m_charge2; // species 2 charge / |electron charge|
      long double m_mass1;  // species 1 mass / electron mass
      long double m_mass2;  // species 2 mass / electron mass
      long double m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
 
      ANGULAR_SCATTERING m_angular_scattering;
   
      // precompute a conversion factor used for each scattering event
      Real m_volume_scale, m_b90_codeToPhys, m_b90_fact; 

      // preallocate some constants to avoid redeclaring many times
      mutable Real m_s12, m_costh, m_sinth, m_phi, m_cosphi, m_sinphi;

      mutable Real m_A12, m_s12sq, m_s12cu;
      inline void setNANBUcosthsinth(const Real a_s12) const 
      {

         // costh = 1 + 1/A*ln[ 1 - rand()*(1-exp(-2*A)) ]
         // and A is found from coth(A) - 1/A = exp(-s12)
         if(a_s12<0.1466) { 
            m_A12 = 1.0/(a_s12*(1.0 - a_s12/2.0 + a_s12*a_s12/6.0));
            m_costh = 1.0 + 1.0/m_A12*std::log(1.0 - MathUtils::rand()*(1.0-std::exp(-2.0*m_A12)));
         }
         else if(a_s12<3.0) {
            m_s12sq = a_s12*a_s12;
            m_s12cu = a_s12*m_s12sq;
            m_A12 = 1.0/(0.0056958 + 0.9560202*a_s12 - 0.508139*m_s12sq
                       + 0.47913906*m_s12cu - 0.12788975*m_s12sq*m_s12sq 
                       + 0.02389567*m_s12cu*m_s12sq);
            m_costh = 1.0 + 1.0/m_A12*std::log(1.0 - MathUtils::rand()*(1.0-std::exp(-2.0*m_A12)));
         }
         else if(a_s12<6.0) {
            m_A12 = 3.0*std::exp(-a_s12);
            m_costh = 1.0 + 1.0/m_A12*std::log(1.0 - MathUtils::rand()*(1.0-std::exp(-2.0*m_A12)));
         }
         else m_costh = 2.0*MathUtils::rand() - 1.0;
         m_sinth = std::sqrt(1.0 - m_costh*m_costh);

      };

};

#include "NamespaceFooter.H"

#endif

