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
               const WeightMethod  a_weight_method,
               const int&        a_verbosity )
         : Scattering(a_verbosity),
           m_Nsubcycle_max(1),
           m_angular_scattering(TAKIZUKA),
           m_weight_method(PROBABILISTIC),
           m_enforce_conservations(false),
           m_sort_weighted_particles(false),
           m_energy_fraction(0.05)
      {
         m_scatter_dt = DBL_MAX;
         m_Clog = 0.0;
         parseSpecies( a_pp );
         parseParameters( a_pp );
      }
      
      Coulomb( const int   a_sp1,
               const int   a_sp2,
               const ParmParse&  a_pp_scatterC,
               const int   a_verbosity )
         : Scattering(a_verbosity),
           m_Nsubcycle_max(1),
           m_angular_scattering(TAKIZUKA),
           m_weight_method(PROBABILISTIC),
           m_enforce_conservations(false),
           m_sort_weighted_particles(false),
           m_energy_fraction(0.05)
      {
         m_scatter_dt = DBL_MAX;
         m_sp1 = a_sp1;
         m_sp2 = a_sp2;
         parseParameters( a_pp_scatterC );
      }

      virtual ~Coulomb() {;}
      
      virtual ScatteringType getScatteringType() const { return COULOMB; }
      
      virtual int species1() const { return m_sp1; }
      
      virtual int species2() const { return m_sp2; }
      
      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const; 
      
      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh );

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
            if(m_weight_method==PROBABILISTIC) {
              std::cout << " weight method:          PROBABILISTIC " << std::endl;
              std::cout << " NxN pairings:           " << (m_NxN?"true":"false") << std::endl;
              if (!m_NxN) {
                std::cout << " NxN Nthresh =           " << m_NxN_Nthresh << std::endl;
              }
              if(m_enforce_conservations) {
                std::cout << " enforce conservations:  true" << std::endl;
                std::cout << "    sort_weighted_particles " << m_sort_weighted_particles << std::endl;
                std::cout << "    energy_fraction         " << m_energy_fraction << std::endl;
              }
              else std::cout << " enforce conservations:  false" << std::endl;
            }
            if(m_weight_method==CONSERVATIVE) {
              std::cout << " weight method:          CONSERVATIVE " << std::endl;
            }
            if(m_Clog>0.0) std::cout << " fixed Clog =            " << m_Clog << std::endl;
            std::cout << " max subcycles           " << m_Nsubcycle_max << std::endl << std::endl;
         }
      }
      
      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const;
      
   private:

      enum ANGULAR_SCATTERING {
         TAKIZUKA,
         NANBU 
      };

      // scatter like species
      void applyIntraScattering_PROB( PicSpecies&            a_picSpecies, 
                                const DomainGrid&            a_mesh,
                                const LevelData<FArrayBox>&  a_LDe_m,
                                const Real                   a_dt_sec ) const;
      
      void applyIntraScattering_PROB_alt( PicSpecies&        a_picSpecies, 
                                    const DomainGrid&            a_mesh,
                                    const LevelData<FArrayBox>&  a_LDe_m,
                                    const Real                   a_dt_sec ) const;
      
      void applyIntraScattering_SK08( PicSpecies&            a_picSpecies, 
                                const DomainGrid&            a_mesh,
                                const LevelData<FArrayBox>&  a_LDe_m,
                                const Real                   a_dt_sec ) const;
      
      // scatter unlike species
      void applyInterScattering_PROB( PicSpecies&            a_picSpecies1,
                                      PicSpecies&            a_picSpecies2,
                                const DomainGrid&            a_mesh,
                                const LevelData<FArrayBox>&  a_LDe_m,
                                const Real                   a_dt_sec ) const;
      
      void applyInterScattering_SK08( PicSpecies&            a_picSpecies1,
                                      PicSpecies&            a_picSpecies2,
                                const DomainGrid&            a_mesh,
                                const LevelData<FArrayBox>&  a_LDe_m,
                                const Real                   a_dt_sec ) const;
      
      void setIntraMFT( const LevelData<FArrayBox>&  a_debyeLength,
                        const LevelData<FArrayBox>&  a_numDensity,
                        const LevelData<FArrayBox>&  a_momDensity,
                        const LevelData<FArrayBox>&  a_eneDensity ) const;

      void setInterMFT( const LevelData<FArrayBox>&  a_debyeLength,
                        const LevelData<FArrayBox>&  a_numDensity1,
                        const LevelData<FArrayBox>&  a_momDensity1,
                        const LevelData<FArrayBox>&  a_eneDensity1,
                        const LevelData<FArrayBox>&  a_numDensity2,
                        const LevelData<FArrayBox>&  a_momDensity2,
                        const LevelData<FArrayBox>&  a_eneDensity2 ) const;
      
      void GalileanScatter( std::array<Real,3>&  a_deltaU,
                      const std::array<Real,3>&  a_vp1,
                      const std::array<Real,3>&  a_vp2,
                      const Real                 a_den12,
                      const Real                 a_bmax,
                      const Real                 a_dt_sec ) const;
     
      void LorentzScatter( std::array<Real,3>&  a_up1,
                           std::array<Real,3>&  a_up2,
                     const bool                 a_scatter2,
                     const long double          a_mass1,
                     const long double          a_mass2,
                     const Real                 a_den12,
                     const Real                 a_bmax,
                     const Real                 a_dt_sec ) const;

      // prohibit copying
      //Coulomb( const Constant& );
      //Coulomb& operator=( const Constant& );
      
      // Parse the species
      void parseSpecies( const ParmParse&  a_pp ) 
      {
         a_pp.get( "speciesA", m_sp1 ); 
         a_pp.get( "speciesB", m_sp2 );
      }

      // Parse the input database for parameters.
      void parseParameters( const ParmParse&  a_pp ) 
      {
        // a_pp.get( "speciesA", m_sp1 ); 
        // a_pp.get( "speciesB", m_sp2 );
         if(a_pp.contains("coulomb_logarithm")) {
            a_pp.get( "coulomb_logarithm", m_Clog ); 
            CH_assert(m_Clog>=2.0);
         }
         a_pp.query( "max_subcycles", m_Nsubcycle_max ); CH_assert(m_Nsubcycle_max>0);
         a_pp.query( "NxN", m_NxN );
         a_pp.query( "NxN_Nthresh", m_NxN_Nthresh );
        
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
         
         // look for specified weighted particle method
         std::string weight_method_str;
         a_pp.query("weight_method",weight_method_str);
         if(!weight_method_str.empty()) {
            if(weight_method_str=="PROBABILISTIC" || weight_method_str=="probabilistic") {
               m_weight_method = PROBABILISTIC;
            }
            else if(weight_method_str=="CONSERVATIVE" || weight_method_str=="conservative") {
               m_weight_method = CONSERVATIVE;
            }
            else {
               cout << "weight_method = " << weight_method_str << endl;
               cout << " valid options: 1) PROBABILISTIC, 2) CONSERVATIVE " << endl;
               MayDay::Error( "Coulomb: invalid weight_method" );
            }
         }

         // these parameters are used for energy-conserving weighted particle method 
         a_pp.query("enforce_conservations",m_enforce_conservations);
         if(m_enforce_conservations) {
            a_pp.query("sort_weighted_particles",m_sort_weighted_particles);
            a_pp.query("energy_fraction",m_energy_fraction);
         }

      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      Real m_Clog;       // Coulomb Logarithm
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      signed int m_charge1; // species 1 charge / |electron charge|
      signed int m_charge2; // species 2 charge / |electron charge|
      long double m_mass1;  // species 1 mass / electron mass
      long double m_mass2;  // species 2 mass / electron mass
      long double m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
      int m_Nsubcycle_max;  // maximum number of subcycling for weighted particle pairs
      
      bool m_NxN = false;
      int m_NxN_Nthresh = 11;

      ANGULAR_SCATTERING m_angular_scattering;
      WeightMethod m_weight_method;
   
      // precompute a conversion factor used for each scattering event
      Real m_b90_fact, m_bqm_fact; 

      // member variables for energy-conserving method for weighted particle collisions
      bool m_enforce_conservations;
      bool m_sort_weighted_particles;
      Real m_energy_fraction;
      int m_loop_count_max = 20; // hard coded for now
      int m_conservation_Nmin_save = 20;

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
      
      inline void modEnergyPairwise( std::array<Real,3>&  a_betap1,
                                     std::array<Real,3>&  a_betap2,
                               const long double          a_wpmp1,
                               const long double          a_wpmp2,
                                    const double          a_mult,
                                     long double&         a_deltaE ) const
      { 
         // method to add/subtract a small percent of energy from a 
         // pair of particles such that momentum is not disturbed
         //
         // This is done by treating it as an inelastic scattering event
         // with potential U = a_deltaE and zero scattering angle
         // a_deltaE > 0 ==> need to take energy away
         // a_deltaE < 0 ==> need to add energy

         const long double muR = a_wpmp1*a_wpmp2/(a_wpmp1 + a_wpmp2);
         int sign = 1;
         if (a_deltaE<0.0) { sign = -1; }

         const Real ux = a_betap1[0] - a_betap2[0];
         const Real uy = a_betap1[1] - a_betap2[1];
         const Real uz = a_betap1[2] - a_betap2[2];
         const long double usq = ux*ux + uy*uy + uz*uz;
         const long double Erel = muR/2.0*usq;
         if (Erel==0.0) return;

         long double deltaE = sign*m_energy_fraction*a_mult*Erel;
         if (std::abs(deltaE) > std::abs(a_deltaE)) {
            deltaE = a_deltaE;
            a_deltaE = 0.0;
         }
         else { a_deltaE -= deltaE; }
         const long double uprime_over_u = std::sqrt(1.0 - deltaE/Erel);

         std::array<Real,3> deltaU;
         deltaU[0] = uprime_over_u*ux - ux;
         deltaU[1] = uprime_over_u*uy - uy;
         deltaU[2] = uprime_over_u*uz - uz;

         for (int n=0; n<3; n++) {
            a_betap1[n] += muR/a_wpmp1*deltaU[n];
            a_betap2[n] -= muR/a_wpmp2*deltaU[n];
         }

      }
                     
      inline void enforceEnergyConservation( std::array<Real,3>& a_betap, 
                                       const Real a_mass, 
                                       const Real a_Eafter ) const
      { 

         // See Sentoku and Kemp JCP 227, 6846-6861 (2008)              
         // add correction to betap with direction randomly oriented in plane normal to betap
         // and with magnitude such that energy is conserved

         std::array<Real,3> dbeta_perp;
         Real betap_r = a_betap[0]*a_betap[0] + a_betap[1]*a_betap[1];
         Real betap_mag = betap_r + a_betap[2]*a_betap[2];
         const Real Eafter2 = a_mass/2.0*betap_mag;
         betap_mag = std::sqrt(betap_mag);
         betap_r = std::sqrt(betap_r);
         if(a_Eafter<Eafter2) return; // can happen for particles where weights are off to machine precision 
         const Real dbeta_perp_mag = std::sqrt(2.0/a_mass*(a_Eafter - Eafter2));
         m_phi = Constants::TWOPI*MathUtils::rand();
         m_cosphi = cos(m_phi);
         m_sinphi = sin(m_phi);
         dbeta_perp[0] = (a_betap[2]*a_betap[0]*m_cosphi - betap_mag*a_betap[1]*m_sinphi)/betap_r*dbeta_perp_mag/betap_mag;
         dbeta_perp[1] = (a_betap[2]*a_betap[1]*m_cosphi + betap_mag*a_betap[0]*m_sinphi)/betap_r*dbeta_perp_mag/betap_mag;
         dbeta_perp[2] = -betap_r*m_cosphi*dbeta_perp_mag/betap_mag;

         // add delta beta correction
         for (int n=0; n<3; n++) { a_betap[n] += dbeta_perp[n]; }

      }

};

#include "NamespaceFooter.H"

#endif

