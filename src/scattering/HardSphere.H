#ifndef _HARDSPHERE_H_
#define _HARDSPHERE_H_

#include "Scattering.H"
#include "PicnicConstants.H"
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

      HardSphere( ParmParse&    a_pp, 
            const WeightMethod  a_weight_method,
            const int           a_verbosity )
         : Scattering(a_verbosity),
           m_weight_method(a_weight_method)
      {
         m_scatter_dt = DBL_MAX;
         parseParameters( a_pp );
      }

      virtual ~HardSphere() {;}
      
      virtual ScatteringType getScatteringType() const { return ELASTIC; }
      
      virtual int species1() const { return m_sp1; }
      
      virtual int species2() const { return m_sp2; }
      
      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const;
      
      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "HardSphere scattering parameters:" << std::endl;
            std::cout << "  species A name:          " << m_species1_name << std::endl;
            std::cout << "  species B name:          " << m_species2_name << std::endl;
            std::cout << "  species 1 radius [m]:    " << m_r1 << std::endl;
            std::cout << "  species 2 radius [m]:    " << m_r2 << std::endl;
            std::cout << "  sigmaT = pi(r1+r2)^2:    " << m_sigmaT << std::endl;
            std::cout << "  weight method:           " << m_weight_method << std::endl;
            std::cout << endl;
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
              MayDay::Error( "HardSphere: invalid weight_method" );
           }
         }

      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      WeightMethod m_weight_method;
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

