#ifndef _HARDSPHERE_H_
#define _HARDSPHERE_H_

#include "Scattering.H"
#include "Constants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * standard hard sphere scattering class for elastic neutral-neutral binary scattering
 * for species1 with radius r1 and species 2 with radius r2. r1 and r2 can be specified
 * from the input file or by default is given as r=(A)^(1/3)*a_Bohr with A the mass of 
 * the species in units of atomic mass unit and a_Bohr = 5.2918e-11 m the Bohr radius.
 *
 * total cross section [m^2]: sigmaT = pi*(r1 + r2)^2
*/
class HardSphere : public Scattering
{
   public:

      // Constructor
      //
      HardSphere( ParmParse&  a_pp, 
                    const int&        a_verbosity )
         : Scattering(a_verbosity)
      {
         m_scatter_dt = DBL_MAX;
         parseParameters( a_pp );
      }

      // Destructor
      //
      virtual ~HardSphere() {;}
      
      // return integer identifier of species 1
      virtual int species1() const { return m_sp1; }
      
      // return integer identifier of species 1
      virtual int species2() const { return m_sp2; }
      
      // set mean free time
      //
      virtual void setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity,
                                    const LevelData<FArrayBox>&  a_energyDensity ) const; 
      
      virtual void setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity1,
                                    const LevelData<FArrayBox>&  a_energyDensity1,
                                    const LevelData<FArrayBox>&  a_numberDensity2,
                                    const LevelData<FArrayBox>&  a_energyDensity2 ) const;
    
      // initialization
      //
      virtual void initialize( const DomainGrid&         a_mesh,
                               const PicSpeciesPtrVect&  a_picSpeciesPtrVect ); 

      virtual Real scatterDt() const { return m_scatter_dt; }

      // Print object parameters.
      //
      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "HardSphere scattering parameters:" << std::endl;
            std::cout << "  species A name:          " << m_species1_name << std::endl;
            std::cout << "  species B name:          " << m_species2_name << std::endl;
            std::cout << "  species 1 radius [m]:    " << m_r1 << std::endl;
            std::cout << "  species 2 radius [m]:    " << m_r2 << std::endl;
            std::cout << "  sigmaT = pi(r1+r2)^2:    " << m_sigmaT << std::endl << std::endl;
         }
      }
      
      // scatter like species
      virtual void applySelfScattering( PicSpecies&  a_picSpecies, 
                                  const DomainGrid&  a_mesh,
                                  const Real         a_dt_sec ) const;
      
      // scatter unlike species
      virtual void applyInterScattering( PicSpecies&  a_picSpecies1,
                                         PicSpecies&  a_picSpecies2,
                                   const DomainGrid&  a_mesh,
                                   const Real         a_dt_sec ) const;
 
 
   private:

      // prohibit copying
      //HardSphere( const Constant& );
      //HardSphere& operator=( const Constant& );

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp ) 
      {
         a_pp.get( "speciesA", m_sp1 ); 
         a_pp.get( "speciesB", m_sp2 );
         a_pp.query( "radius1", m_r1 );
         a_pp.query( "radius2", m_r2 );
      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      Real m_mass1;  // species 1 mass / electron mass
      Real m_mass2;  // species 1 mass / electron mass
      Real m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
      Real m_r1=-1;  // species 1 radius [m]
      Real m_r2=-1;  // species 2 radius [m]
      Real m_sigmaT; // total cross section [m^2] = pi*(r1+r2)^2

};

#include "NamespaceFooter.H"

#endif

