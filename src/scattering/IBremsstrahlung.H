#ifndef _INVERSE_BREMSSTRAHLUNG_H_
#define _INVERSE_BREMSSTRAHLUNG_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * IBremsstrahlung scattering class:
 * photon + e  ==> e
*/

class IBremsstrahlung : public Scattering
{
   public:

      IBremsstrahlung( ParmParse&  a_pp,
        const int&        a_verbosity )
         : Scattering(a_verbosity)
      {
         m_scatter_dt = DBL_MAX;
         m_deltaE_IBremsstrahlung = 0.0;
         parseParameters( a_pp );
      }

      virtual ~IBremsstrahlung() = default;

      virtual ScatteringType getScatteringType() const { return IBREMSSTRAHLUNG; }

      virtual int species1() const { return m_sp1; }

      virtual int species2() const { return m_sp2; }

      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const;

      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }

      virtual Real getDeltaEIBremsstrahlung() const { return m_deltaE_IBremsstrahlung; }
      virtual void zeroDeltaEIBremsstrahlung() const { m_deltaE_IBremsstrahlung=0.0; }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "IBremsstrahlung scattering parameters:" << std::endl;
            std::cout << "  species 1 name:   " << m_species1_name << std::endl;
            std::cout << "  species 2 name:   " << m_species2_name << std::endl;
            std::cout << std::endl;
         }
      }

      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const;

   private:

      void applyIBremsstrahlung( PicPhotonSpecies&   a_picSpecies1,
                                 PicChargedSpecies&  a_picSpecies2,
                           const LevelData<FArrayBox>&  a_nuei,
                           const DomainGrid&  a_mesh,
                           const Real         a_dt_sec ) const;

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp )
      {

          // get number associated with each species
          a_pp.get( "species1", m_sp1 );
          a_pp.get( "species2", m_sp2 );
#ifdef RELATIVISTIC_PARTICLES
#else
          MayDay::Error( "IBremsstrahlung scattering model requires RELATIVISTIC_PARTICLES" );
#endif

      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt;    // (mean free time)
      mutable Real m_deltaE_IBremsstrahlung; // running sum of absorbed bremss energy
      int m_sp1;     // species 1 identifier (photon)
      int m_sp2;     // species 2 identifier (electron)
      Real m_mass1;  // photon mass / electron mass
      Real m_mass2;  // species 2 mass / electron mass
      Real m_energy_fraction = 1.0e-3; // this parameter is used for the electron energy absorbing method

      Real m_mcSq_eV;

};

#include "NamespaceFooter.H"

#endif

