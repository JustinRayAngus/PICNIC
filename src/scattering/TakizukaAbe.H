#ifndef _TAKIZUKAABE_H_
#define _TAKIZUKAABE_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * Takizuka and Abe scattering model for Coulomb collisions in a plasma amongst
 * equally-weighted macro particles
*/
class TakizukaAbe : public Scattering
{
   public:

      TakizukaAbe( ParmParse&  a_pp, 
                    const int&        a_verbosity )
         : Scattering(a_verbosity)
      {
         m_scatter_dt = DBL_MAX;
         m_Clog = 10.0;
#ifdef RELATIVISTIC_PARTICLES
         m_b90_codeToPhys = Constants::QE*Constants::QE/(2.0*Constants::PI*Constants::EP0*Constants::ME);
#else
         m_b90_codeToPhys = Constants::QE*Constants::QE/(4.0*Constants::PI*Constants::EP0*Constants::ME); 
#endif
         parseParameters( a_pp );
      }
      
      TakizukaAbe( const int   a_sp1,
                   const int   a_sp2,
                   const Real  a_Clog, 
                   const int   a_verbosity )
         : Scattering(a_verbosity)
      {
         m_scatter_dt = DBL_MAX;
#ifdef RELATIVISTIC_PARTICLES
         m_b90_codeToPhys = Constants::QE*Constants::QE/(2.0*Constants::PI*Constants::EP0*Constants::ME); 
#else
         m_b90_codeToPhys = Constants::QE*Constants::QE/(4.0*Constants::PI*Constants::EP0*Constants::ME); 
#endif
         //
         m_sp1 = a_sp1;
         m_sp2 = a_sp2;
         m_Clog = a_Clog; CH_assert(m_Clog>=3.0);
      }

      virtual ~TakizukaAbe() {;}
      
      virtual ScatteringType getScatteringType() const { return COULOMB; }
      
      virtual int species1() const { return m_sp1; }
      
      virtual int species2() const { return m_sp2; }
      
      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const;
      
      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "TakizukaAbe scattering parameters:" << std::endl;
            std::cout << "  species A name:      " << m_species1_name << std::endl;
            std::cout << "  species B name:      " << m_species2_name << std::endl;
            std::cout << "  Clog =               " << m_Clog << std::endl << std::endl;
         }
      }
      
      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const;
      
   private:

      // scatter like species
      void applySelfScattering( PicSpecies&  a_picSpecies, 
                          const DomainGrid&  a_mesh,
                          const Real         a_dt_sec ) const;
      
      // scatter unlike species
      void applyInterScattering( PicSpecies&  a_picSpecies1,
                                 PicSpecies&  a_picSpecies2,
                           const DomainGrid&  a_mesh,
                           const Real         a_dt_sec ) const;
      
      void setIntraMFT( const LevelData<FArrayBox>&  a_numberDensity,
                        const LevelData<FArrayBox>&  a_energyDensity ) const;

      void setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                        const LevelData<FArrayBox>&  a_energyDensity1,
                        const LevelData<FArrayBox>&  a_numberDensity2,
                        const LevelData<FArrayBox>&  a_energyDensity2 ) const;
      
      void computeDeltaU( std::array<Real,3>&  a_deltaU,
                    const std::array<Real,3>&  a_vp1,
                    const Real                 a_den1,
                    const std::array<Real,3>&  a_vp2,
                    const Real                 a_den2,
                    const Real                 a_Clog,
                    const Real                 a_dt_sec ) const;
     
      void LorentzScatter( std::array<Real,3>&  a_up1,
                           std::array<Real,3>&  a_up2,
                     const long double          a_mass1,
                     const long double          a_mass2,
                     const Real                 a_den2,
                     const Real                 a_dt_sec ) const;

      // prohibit copying
      //TakizukaAbe( const Constant& );
      //TakizukaAbe& operator=( const Constant& );

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp ) 
      {
         a_pp.get( "speciesA", m_sp1 ); 
         a_pp.get( "speciesB", m_sp2 );
         a_pp.query( "coulomb_logarithm", m_Clog ); CH_assert(m_Clog>=3.0);
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
   
      // precompute a conversion factor used for each scattering event
      Real m_b90_codeToPhys, m_b90_fact; 

      // preallocate some constants to avoid redeclaring many times
      mutable Real m_delta, m_deltasq;
      mutable Real m_costh, m_sinth, m_phi, m_cosphi, m_sinphi;

};

#include "NamespaceFooter.H"

#endif

