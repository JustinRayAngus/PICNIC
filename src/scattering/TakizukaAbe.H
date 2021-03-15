#ifndef _TAKIZUKAABE_H_
#define _TAKIZUKAABE_H_

#include "Scattering.H"
#include "Constants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * Takizuka and Abe scattering model for Coulomb collisions in a plasma amongtst
 * equally-weighted macro particles
 *
 * total cross section [m^2]: sigmaT = pi*(r1 + r2)^2
*/
class TakizukaAbe : public Scattering
{
   public:

      // Constructor
      //
      TakizukaAbe( ParmParse&  a_pp, 
                    const int&        a_verbosity )
         : Scattering(a_verbosity)
      {
         m_scatter_dt = DBL_MAX;
          m_b90_codeToPhys = Constants::QE*Constants::QE/(4.0*Constants::PI*Constants::EP0*Constants::ME); 
         parseParameters( a_pp );
      }

      // Destructor
      //
      virtual ~TakizukaAbe() {;}
      
      // return integer identifier of species 1
      virtual int species1() const { return m_sp1; }
      
      // return integer identifier of species 1
      virtual int species2() const { return m_sp2; }
      
      // set mean free time
      //
      virtual void setMeanFreeTime( const DomainGrid&            a_mesh,
                                    const LevelData<FArrayBox>&  a_numberDensity,
                                    const LevelData<FArrayBox>&  a_energyDensity ) ;
      
      // initialization
      //
      virtual void initialize( const DomainGrid&         a_mesh,
                               const PicSpeciesPtrVect&  a_picSpeciesPtrVect ); 

      virtual Real scatterDt() const { return m_scatter_dt; }

      // Print object parameters.
      //
      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "TakizukaAbe scattering parameters:" << std::endl;
            std::cout << "  species 1 name:          " << m_species1_name << std::endl;
            std::cout << "  species 2 name:          " << m_species2_name << std::endl;
         }
      }
      
      // scatter like species
      virtual void applySelfScattering( PicSpecies&            a_picSpecies, 
                                  const DomainGrid&            a_mesh,
                                  const LevelData<FArrayBox>&  a_numberDensity,
                                  const LevelData<FArrayBox>&  a_energyDensity,
                                  const Real                   a_dt ) const;
      
      // scatter unlike species
      virtual void applyInterScattering( PicSpecies&            a_picSpecies1,
                                         PicSpecies&            a_picSpecies2,
                                   const DomainGrid&            a_mesh,
                                   const Real                   a_dt ) const;
 
 
   private:


      void computeDeltaU( std::array<Real,3>&  a_deltaU,
                    const std::array<Real,3>&  a_vp1,
                    const Real                 a_den1,
                    const std::array<Real,3>&  a_vp2,
                    const Real                 a_den2,
                    const Real                 a_Clog,
                    const Real                 a_dt ) const;

      // prohibit copying
      //TakizukaAbe( const Constant& );
      //TakizukaAbe& operator=( const Constant& );

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp ) 
      {
         a_pp.get( "species1", m_sp1 ); m_sp1 = m_sp1-1; 
         a_pp.get( "species2", m_sp2 ); m_sp2 = m_sp2-1;
      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      signed int m_charge1; // species 1 charge / |electron charge|
      signed int m_charge2; // species 2 charge / |electron charge|
      Real m_mass1;  // species 1 mass / electron mass
      Real m_mass2;  // species 2 mass / electron mass
      Real m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
   
      // precompute a conversion factor used for each scattering event
      Real m_b90_codeToPhys; 

};

#include "NamespaceFooter.H"

#endif

