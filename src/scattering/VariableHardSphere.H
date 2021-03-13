#ifndef _VARIABLEHARDSPHERE_H_
#define _VARIALBEHARDSPHERE_H_

#include "Scattering.H"
#include "Constants.H"
#include "PicSpecies.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * variable hard sphere scattering class for elastic neutral-neutral binary scattering
 * See Viscosity of Gases by marcia L. Huber for temp dependent viscosity of common gases
 * See Nanbu 2000 pg 982 and Hong Xiao 2014 for VHS model
 *
 * total cross section: sigmaT = 4*pi*Aconst*g^(-4/alpha), g is relative speed [m/s]
*/
class VariableHardSphere : public Scattering
{
   public:

      // Constructor
      //
      VariableHardSphere( ParmParse&  a_pp, 
                    const int&        a_verbosity )
         : Scattering(a_verbosity)
      {
         m_scatter_dt = DBL_MAX;
         parseParameters( a_pp );
      }

      // Destructor
      //
      virtual ~VariableHardSphere() {;}
      
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
            std::cout << "VariableHardSphere scattering parameters:" << std::endl;
            std::cout << "  species 1 name:              " << m_species1_name << std::endl;
            std::cout << "  species 2 name:              " << m_species2_name << std::endl;
            std::cout << "  reference temperature [K]  : " << m_T0 << std::endl;
            std::cout << "  reference viscosity [Pa-s] : " << m_mu0 << std::endl;
            std::cout << "  mu=mu0*(T/T0)^eta with eta : " << m_eta << std::endl;
            std::cout << "  Aconst :                     " << m_Aconst << std::endl;
         }
      }
      
      virtual void applySelfScattering( PicSpecies&            a_picSpecies, 
                                  const DomainGrid&            a_mesh,
                                  const LevelData<FArrayBox>&  a_numberDensity,
                                  const LevelData<FArrayBox>&  a_energyDensity,
                                  const Real                   a_dt ) const;
 
   private:

      // prohibit copying
      //VariableHardSphere( const Constant& );
      //VariableHardSphere& operator=( const Constant& );

      /// Parse the input database for parameters.
      /**
       */
      void parseParameters( ParmParse&  a_pp ) 
      {
         a_pp.get( "species1", m_sp1 ); m_sp1 = m_sp1-1; 
         a_pp.get( "species2", m_sp2 ); m_sp2 = m_sp2-1;
         a_pp.get( "reference_temperature", m_T0 );
         a_pp.get( "reference_viscosity", m_mu0 );
         a_pp.get( "eta", m_eta );
      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      Real m_sp1;    // species 1
      Real m_sp2;    // species 2
      Real m_T0;     // reference temperature [K]
      Real m_mu0;    // reference viscosity [Pa-s]
      Real m_eta;    // mu = mu0*(T/T0)^eta
      Real m_alpha;  // alpha = 4/(2*eta-1)
      Real m_mass1;  // neutral mass / electron mass
      Real m_mass2;  // neutral mass / electron mass
      Real m_Aconst; // sigmaT = 4*pi*Aconst*g^(-4/alpha), g is relative speed [m/s]

};

#include "NamespaceFooter.H"

#endif

