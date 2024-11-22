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
           m_angular_scattering(TAKIZUKA),
           m_weight_method(PROBABILISTIC),
           m_energy_fraction(0.05),
           m_energy_fraction_max(0.5)
      {
         m_scatter_dt = DBL_MAX;
         parseSpecies( a_pp );
         parseParameters( a_pp );
      }
      
      Coulomb( const int   a_sp1,
               const int   a_sp2,
               const ParmParse&  a_pp_scatterC,
               const int   a_verbosity )
         : Scattering(a_verbosity),
           m_angular_scattering(TAKIZUKA),
           m_weight_method(PROBABILISTIC),
           m_energy_fraction(0.05),
           m_energy_fraction_max(0.5)
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
            if (m_angular_scattering==TAKIZUKA) {
              std::cout << " angular scattering:     TAKIZUKA " << std::endl;
            }
            else if (m_angular_scattering==NANBU) {
              std::cout << " angular scattering:     NANBU " << std::endl;
            }
            else if (m_angular_scattering==NANBU_FAS) {
              std::cout << " angular scattering:     NANBU_FAS " << std::endl;
              std::cout << "   limit_fas_theta_max:   " << (m_limit_fas_theta_max?"true":"false") << std::endl;
              std::cout << "   exclude electron fas:  " << (m_exclude_electron_fas?"true":"false") << std::endl;
            }
            else if (m_angular_scattering==NANBU_FAS_v2) {
              std::cout << " angular scattering:     NANBU_FAS_v2 " << std::endl;
              std::cout << "   limit_fas_theta_max:   " << (m_limit_fas_theta_max?"true":"false") << std::endl;
              std::cout << "   exclude electron fas:  " << (m_exclude_electron_fas?"true":"false") << std::endl;
            }
            else if (m_angular_scattering==ISOTROPIC) {
              std::cout << " angular scattering:     ISOTROPIC " << std::endl;
            }
            if (m_weight_method==PROBABILISTIC) {
              std::cout << " weight method:          PROBABILISTIC " << std::endl;
              std::cout << " NxN pairings:           " << (m_NxN?"true":"false") << std::endl;
              if (!m_NxN) {
                std::cout << " NxN Nthresh =           " << m_NxN_Nthresh << std::endl;
              }
              if (m_enforce_conservations) {
                std::cout << " enforce conservations:  true" << std::endl;
                std::cout << "    print_correction_count  " << (m_print_correction_count?"true":"false") << std::endl;
                std::cout << "    sort_weighted_particles " << (m_sort_weighted_particles?"true":"false") << std::endl;
                std::cout << "    energy_fraction =       " << m_energy_fraction << std::endl;
                std::cout << "    energy_fraction_max =   " << m_energy_fraction_max << std::endl;
                std::cout << "    beta_weight_exponent =  " << m_beta_weight_exponent << std::endl;
              }
              else std::cout << " enforce conservations:  false" << std::endl;
            }
            if (m_weight_method==CONSERVATIVE) {
              std::cout << " weight method:          CONSERVATIVE " << std::endl;
            }
            if (m_Clog>0.0) {
               std::cout << " fixed Clog =            " << m_Clog << std::endl;
            }
            std::cout << " num subcycles =         " << m_Nsubcycles << std::endl << std::endl;
         }
      }
      
      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const;
      
   private:

      enum ANGULAR_SCATTERING {
         TAKIZUKA,
         NANBU,
         NANBU_FAS,
         NANBU_FAS_v2,
         ISOTROPIC
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
                      const Real                 a_sigma_max,
                      const Real                 a_dt_sec ) const;
     
      void LorentzScatter( std::array<Real,3>&  a_up1,
                           std::array<Real,3>&  a_up2,
                     const bool                 a_scatter2,
                     const long double          a_mass1,
                     const long double          a_mass2,
                     const Real                 a_den12,
                     const Real                 a_bmax,
                     const Real                 a_sigma_max,
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
         if (a_pp.contains("coulomb_logarithm")) {
            a_pp.get( "coulomb_logarithm", m_Clog ); 
            CH_assert(m_Clog>=2.0);
         }
         a_pp.query( "num_subcycles", m_Nsubcycles );
         CH_assert(m_Nsubcycles>0);
         a_pp.query( "NxN", m_NxN );
         a_pp.query( "NxN_Nthresh", m_NxN_Nthresh );
        
         // look for specified scattering model
         std::string scatter_type_str;
         a_pp.query("angular_scattering",scatter_type_str);
         if (!scatter_type_str.empty()) {
            if (scatter_type_str=="TAKIZUKA" || scatter_type_str=="takizuka") {
               m_angular_scattering = TAKIZUKA;
            }
            else if (scatter_type_str=="NANBU" || scatter_type_str=="nanbu") {
               m_angular_scattering = NANBU;
            }
            else if (scatter_type_str=="NANBU_FAS" || scatter_type_str=="nanbu_fas") {
               m_angular_scattering = NANBU_FAS;
            }
            else if (scatter_type_str=="NANBU_FAS_v2" || scatter_type_str=="NANBU_FAS_JRA") {
               // NANBU_FAS_JRA is for legacy purposes
               m_angular_scattering = NANBU_FAS_v2;
            }
            else if (scatter_type_str=="isotropic") {
               m_angular_scattering = ISOTROPIC;
            }
            else {
               cout << "angular_scattering = " << scatter_type_str << endl;
               cout << " valid options: 1) takizuka, 2) nanbu " << endl;
               MayDay::Error( "Coulomb: invalid angular_scattering" );
            }
         }

         if (m_angular_scattering==NANBU_FAS || m_angular_scattering==NANBU_FAS_v2) {
            a_pp.query( "limit_fas_theta_max", m_limit_fas_theta_max );
            a_pp.query( "exclude_electron_fas", m_exclude_electron_fas );
         }
         
         // look for specified weighted particle method
         std::string weight_method_str;
         a_pp.query("weight_method",weight_method_str);
         if (!weight_method_str.empty()) {
            if (weight_method_str=="PROBABILISTIC" || weight_method_str=="probabilistic") {
               m_weight_method = PROBABILISTIC;
            }
            else if (weight_method_str=="CONSERVATIVE" || weight_method_str=="conservative") {
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
         if (m_enforce_conservations) {
            a_pp.query("print_correction_count",m_print_correction_count);
            a_pp.query("sort_weighted_particles",m_sort_weighted_particles);
            a_pp.query("energy_fraction",m_energy_fraction);
            a_pp.query("energy_fraction_max",m_energy_fraction_max);
            a_pp.query("beta_weight_exponent",m_beta_weight_exponent);
            CH_assert( m_beta_weight_exponent > 0 );
         }

      }

      std::string m_species1_name;
      std::string m_species2_name;
      mutable Real m_scatter_dt; // (mean free time)
      Real m_Clog = 0.0;       // Coulomb Logarithm
      int m_sp1;     // species 1 identifier
      int m_sp2;     // species 2 identifier
      signed int m_charge1; // species 1 charge / |electron charge|
      signed int m_charge2; // species 2 charge / |electron charge|
      long double m_mass1;  // species 1 mass / electron mass
      long double m_mass2;  // species 2 mass / electron mass
      long double m_mu;     // reduced mass / electron mass = mass1*mass2/(mass1+mass2)
      int m_Nsubcycles = 1;  // number of subcycling for scattering
      
      bool m_NxN = false;
      int m_NxN_Nthresh = 11;

      mutable LevelData<FArrayBox> m_count;

      ANGULAR_SCATTERING m_angular_scattering;
      WeightMethod m_weight_method;
   
      // precompute a conversion factor used for each scattering event
      Real m_b90_fact, m_bqm_fact;
      Real m_EF_fact = 0.0;
      mutable Real m_EF_norm;

      // member variables for energy-conserving method for weighted particle collisions
      bool m_enforce_conservations = false;
      bool m_print_correction_count = false;
      bool m_sort_weighted_particles = false;
      bool m_limit_fas_theta_max = false;
      bool m_exclude_electron_fas = true;
      Real m_energy_fraction;
      Real m_energy_fraction_max;
      int m_beta_weight_exponent = 1;
      int m_conservation_Nmin_save = 100000;

      // preallocate some constants to avoid redeclaring many times
      mutable Real m_s12, m_costh, m_sinth, m_phi, m_cosphi, m_sinphi;

      mutable Real m_A12, m_s12sq, m_s12cu;
      inline void setNANBUcosthsinth(const Real a_s12) const 
      {

         // costh = 1 + 1/A*ln[ 1 - rand()*(1-exp(-2*A)) ]
         // and A is found from coth(A) - 1/A = exp(-s12)
         if (a_s12<0.1466) {
            m_A12 = 1.0/(a_s12*(1.0 - a_s12/2.0 + a_s12*a_s12/6.0));
            m_costh = 1.0 + 1.0/m_A12*std::log(1.0 - MathUtils::rand()*(1.0-std::exp(-2.0*m_A12)));
         }
         else if (a_s12<3.0) {
            m_s12sq = a_s12*a_s12;
            m_s12cu = a_s12*m_s12sq;
            m_A12 = 1.0/(0.0056958 + 0.9560202*a_s12 - 0.508139*m_s12sq
                       + 0.47913906*m_s12cu - 0.12788975*m_s12sq*m_s12sq 
                       + 0.02389567*m_s12cu*m_s12sq);
            m_costh = 1.0 + 1.0/m_A12*std::log(1.0 - MathUtils::rand()*(1.0-std::exp(-2.0*m_A12)));
         }
         else if (a_s12<6.0) {
            m_A12 = 3.0*std::exp(-a_s12);
            m_costh = 1.0 + 1.0/m_A12*std::log(1.0 - MathUtils::rand()*(1.0-std::exp(-2.0*m_A12)));
         }
         else m_costh = 2.0*MathUtils::rand() - 1.0;
         m_sinth = std::sqrt(1.0 - m_costh*m_costh);

      }

      inline void setNANBUFAScosthsinth( const Real  s12,
                                         const Real  Clog,
                                         const Real  b0,
                                         const Real  bmin_qm,
                                         const Real  bmax,
                                         const Real  sigma_eff ) const
      {
         // set costh and sinth using NANBU method extended to include single large angle
         // scattering events (i.e., full angle scattering) via Higginson JCP 2017.

         // s12_min (where N = 1) is Clog*b0^2/(bmax^2-bmin_R^2).
         // However, when bmin in Clog is set by bmin_qm rather than b0/2,
         // then the FAS model is not solvable. Using bmin_R set to some fraction
         // of bmin_qm helps with solvability, but not gauranteed. Default is to use
         // the expression by Higginson 2017 for s12_min so the equations are solvable.
         // However, it should be mentioned that this means one will do pure Rutherford
         // even though N > 1 when bmin>b0/2.

         // set s12_min (where N ~= 1). Factor of 1.33 is for solvability at Clog = 2
         Real s12_min = 1.33*4.0*Clog/(std::exp(2.0*Clog)-1.0);

         // set minimum impact parameter for Rutherford scattering.
         // bmin_R = 0.14*bmin_qm matches transition for 27.5 MeV alphas
         // at 60 degrees observed experimentally. Using 0.5 for factor, which works
         // better for electrons. However, really shouldn't use this for electrons.
         Real bmin_R = 0.0;
         Real theta_max = Constants::PI;
         if (m_limit_fas_theta_max && bmin_qm>0.0) {
            bmin_R = 0.5*bmin_qm;
            theta_max = 2.0*std::atan(b0/2.0/bmin_R);
            s12_min = 1.33*Clog*b0*b0/(bmax*bmax-bmin_R*bmin_R);
         }
         const Real costhmax = std::cos(theta_max);
         const Real mu_max = (1.0 - costhmax)/2.0; // sin^2(theta_max/2)

         if (bmin_R>bmax/2.0 || b0>bmax) {
            setNANBUcosthsinth(s12);
            return;
         }

         if (s12 < s12_min) { // do single scatter event
            const Real theta_min = 2.0*std::atan(b0/2.0/bmax);  // minimum scatter angle
            const Real mu_min = (1.0 - std::cos(theta_min))/2.0;
            const Real sigma_R = Constants::PI*b0*b0/4.0*(1.0/mu_min - 1.0/mu_max);
            const Real s12_L = sigma_R/sigma_eff*s12;
            const Real PL = 1.0 - std::exp(-s12_L);
            if (s12_L > 2.0) { setNANBUcosthsinth(s12); }
            else if (MathUtils::rand()<PL) {
               const Real RL = MathUtils::rand();
               const Real costhmin = std::cos(theta_min);
               const Real numer = RL*(costhmin - costhmax) - costhmin*(1.0-costhmax);
               const Real denom = RL*(costhmin - costhmax) - (1.0 - costhmax);
               m_costh = numer/denom;
               m_sinth = std::sqrt(1.0-m_costh*m_costh);
            }
            else {
               m_costh = 1.0;
               m_sinth = 0.0;
            }
         }
         else if (s12 < 0.5) { // use nanbu full angle scattering model

            // set coefficients alpha_g, sA, and muc = sin^2(theta_c/2.0)
            Real alpha_g, sA, muc;
            int success = setFAScoefficients( alpha_g, sA, muc, Clog, s12, mu_max );
            if (success<0 || muc!=muc) {
               cout << "setFAScoefficients() did not have success" << endl;
               cout << "  Clog = " << Clog << endl;
               cout << "  s12_min = " << s12_min << endl;
               cout << "  s12 = " << s12 << endl << endl;
               setNANBUcosthsinth(s12);
            }
            else {
               const Real Uc = 1.0 - s12/(4.0*Clog)*(mu_max-muc)/(mu_max*muc);
               const Real R = MathUtils::rand();
               if (R<Uc) {
                  m_costh = 1.0 + sA*std::log(1.0 - R/Uc*(1.0-std::exp(-2.0*muc/sA)));
                  m_sinth = std::sqrt(1.0-m_costh*m_costh);
               }
               else {
                  const Real costhc = 1.0 - 2.0*muc;
                  const Real R2 = MathUtils::rand();
                  const Real numer = R2*(costhc - costhmax) - costhc*(1.0-costhmax);
                  const Real denom = R2*(costhc - costhmax) - (1.0 - costhmax);
                  m_costh = numer/denom;
                  m_sinth = std::sqrt(1.0-m_costh*m_costh);
               }
            }
         }
         else { // use standard cumulative small angle nanbu model
            setNANBUcosthsinth(s12);
         }
      }

      inline void setNANBUFAS_v2_costhsinth( const Real  s12,
                                             const Real  Clog,
                                             const Real  b0,
                                             const Real  bmin_qm,
                                             const Real  bmax,
                                             const Real  sigma_eff ) const
      {
         // set costh and sinth using probability distribution function that is a linear
         // combination of that used by NANBU (N97) for cumulative small angle events and
         // a modified Rutherford formula for single scattering events. This method is
         // similar to the full angle scattering (FAS) method by Higginson JCP 2017, but
         // 1) Must solve a nonlinear 2-variable system rather than a 3-variable one (easier and faster)
         // 2) P(theta) is continuous and smooth
         // 3) P(theta) merges identically to modified Rutherfod for N = 1
         // 4) No solvability issues when b_min in Clog is not set by b0/2

         // set minimum impact parameter for Rutherford scattering
         // bmin_R = 0.14*bmin_qm matches transition for 27.5 MeV alphas
         // at 60 degrees observed experimentally. 0.5 works better for original FAS.
         Real bmin_R = 0.0;
         Real theta_max = Constants::PI;
         if (m_limit_fas_theta_max && bmin_qm>0.0) {
            bmin_R = 0.5*bmin_qm;
            theta_max = 2.0*std::atan(b0/2.0/bmin_R);
         }
         const Real mu_max = (1.0 - std::cos(theta_max))/2.0; // sin^2(theta_max/2)

         if (bmin_R>bmax/2.0 || b0>bmax || s12>0.6) {
            setNANBUcosthsinth(s12);
            return;
         }

         // set target value of number of Rutherford collisions to start transition region
         // and the scale factor for mu_tr in the modified Rutherford formula.
         // Nmax = 2*Clog and mutr_factor = Clog works well for Clog >= 2.
         Real Nmax = 2.0*Clog;
         const Real mutr_factor = Clog;

         // compute s12_Nmin where N = sigma_R*n*u12*dt = 1
         // and s12_Nmax where N = Nmax
         const Real sig_ratio = b0*b0/4.0/(bmax*bmax - bmin_R*bmin_R); // sig90/sigR
         const Real s12_Nmin = 4.0*Clog*sig_ratio;
         Real s12_Nmax = Nmax*s12_Nmin; // may be overwritten below

         // there is a lower cutoff to Nmax set by variance of Rutherford
         // equal to variance of Nanbu 97:
         // 1 - exp(-s12) ~= s12 = s12/(2*Clog)*log((mu_max+mu_tr)/mu_tr)
         // This only occurs for small Clog and when bmin in Clog is set by
         // bmin_R and not by b0/2. Ensure Nmax is at least twice this value.
         const Real Nmax_min = 2.0*mutr_factor*mu_max/(std::exp(2.0*Clog)-1.0)/s12_Nmin;
         if (Nmax<Nmax_min) {
            Nmax = Nmax_min;
            s12_Nmax = Nmax*s12_Nmin;
         }

         if (s12 < s12_Nmax) { // tansition region to pure Rutherford

            // compute probability that a scattering occurs
            const Real Ntot = s12/s12_Nmin;
            const Real Pscatter = 1.0 - std::exp(-Ntot);
            if (MathUtils::rand()>Pscatter) { // no scattering occurs
               m_costh = 1.0;
               m_sinth = 0.0;
               return;
            }

            // solver for coefficients alpha_g, sA using s12_Nmax
            Real alpha_g, sA;
            Real mu_tr = s12_Nmax/mutr_factor;
            int success = setFAS_v2_coefficients( alpha_g, sA, Clog, s12_Nmax, mu_max, mu_tr );
            if (success<0 || sA!=sA) {
               cout << "setFAS_v2_coefficients failed:" << endl;
               cout << "  Clog     = " << Clog << endl;
               cout << "  s12_Nmax = " << s12_Nmax << endl;
               cout << "  mu_max   = " << mu_max << endl;
               cout << "  mu_tr    = " << mu_tr << endl;
               exit(EXIT_FAILURE);
            }

            // set alpha_g such that it goes to zero at s12 = s12_Nmin
            alpha_g = std::max(0.0,alpha_g*(s12-s12_Nmin)/(s12_Nmax-s12_Nmin));

            // set sA to scale with s12
            if (alpha_g==0.0) { sA = mu_max; }
            else { sA *= s12/s12_Nmax; }

            // compute mu0 = mu_tr where modified Rutherford cross section is
            // sigma_{R,m} = sig90*(1/mu_tr - 1/(mu_max+mu_tr)) = pi*(bmax^2 - bmin_R^2)
            Real mu0;
            const Real coefc = 4.0*sig_ratio/mu_max;
            if (coefc<1.0e-10) { mu0 = sig_ratio; }
            else { mu0 = mu_max*(-1.0 + std::sqrt(1.0 + coefc))/2.0; }

            // compute mu_tr for Rutherford such that S_N97 + S_L = 1
            const Real S_N97 = alpha_g*(1.0-std::exp(-2.0*mu_max/sA));
            if (s12<=s12_Nmin) { mu_tr = mu0; }
            else {
               const Real C0 = (1.0 - S_N97)*4.0*Clog/s12;
               const Real coef = 4.0/mu_max/C0;
               if (coef<1.0e-10) { mu_tr = 1.0/C0; }
               else { mu_tr = mu_max*(-1.0 + std::sqrt(1.0 + coef))/2.0; }
               mu_tr = std::max(mu_tr,mu0);
            }

            // Sample from cumulative P with probability S_N97, else Rutherford
            const Real R = MathUtils::rand();
            if (R<S_N97) { // cumulative scattering
               m_costh = 1.0 + sA*std::log(1.0 - R/S_N97*(1.0-std::exp(-2.0*mu_max/sA)));
            }
            else { // Rutherford scattering
               const Real RL = MathUtils::rand();
               m_costh = 1.0 - 2.0*RL*mu_tr*mu_max/(mu_max*(1.0-RL)+mu_tr);
            }
            m_sinth = std::sqrt(1.0-m_costh*m_costh);

         }
         else { // use full angle scattering model

            // set coefficients alpha_g, sA for Nanbu part of distribution
            Real alpha_g, sA;
            Real mu_tr = s12/mutr_factor;
            int success = setFAS_v2_coefficients( alpha_g, sA, Clog, s12, mu_max, mu_tr );
            if (success<0 || sA!=sA) {
               cout << "setFAS_v2_coefficients failed:" << endl;
               cout << "  Clog   = " << Clog << endl;
               cout << "  s12    = " << s12_Nmax << endl;
               cout << "  mu_max = " << mu_max << endl;
               cout << "  mu_tr  = " << mu_tr << endl;
               exit(EXIT_FAILURE);
            }

            // compute relative contributions of cumulative and Rutherford
            // probability distribution functions
            const Real S_L = s12/(4.0*Clog)*mu_max/mu_tr/(mu_max+mu_tr);
            const Real S_N97 = 1.0 - S_L;

            // determine if scattering is cumulative or Rutherford
            const Real R = MathUtils::rand();
            if (R<S_N97) { // cumulative scattering
               m_costh = 1.0 + sA*std::log(1.0 - R/S_N97*(1.0-std::exp(-2.0*mu_max/sA)));
            }
            else { // Rutherford scattering
               const Real RL = MathUtils::rand();
               m_costh = 1.0 - 2.0*RL*mu_tr*mu_max/(mu_max*(1.0-RL)+mu_tr);
            }
            m_sinth = std::sqrt(1.0-m_costh*m_costh);

         }

      }

      inline int setFAScoefficients( Real& alpha_g,
                                     Real& sA,
                                     Real& muc,
                               const Real  Clog,
                               const Real  s12,
                               const Real  mu_max ) const
      {

         Real Xc, sAold, f1, f2;

         // initial guess for variables
         alpha_g = 1.0;
         sA = s12/2.0;
         Xc = getTransitionX_NANBU( Clog, s12, alpha_g, sA );
         muc = sA*Xc; // sin^2(theta_c/2)

         int iter = 0;
         int success = 1;
         Real error = 1.0;
         while (error>1.0e-4) {

            sAold = sA;
            f1 = 1.0 - std::exp(-s12) + s12/Clog*0.5*std::log(Xc*sAold/mu_max);
            f2 = (1.0-std::exp(-2.0*Xc))/(1.0-(1.0+Xc)*std::exp(-2.0*Xc));
            sA = (4.0*mu_max*Clog/(4.0*mu_max*Clog+s12))*(s12/(4.0*Xc*Clog) + f1*f2);

            //  update muc = sin^2(theta_c/2)
            Xc = getTransitionX_NANBU( Clog, s12, alpha_g, sA );
            muc = sA*Xc; // sin^2(theta_c/2)

            // update alpha_g
            alpha_g = (1.0 - s12/(4.0*Clog)*(mu_max-muc)/mu_max/muc)/(1.0-std::exp(-2.0*muc/sA));

            // check error
            error = std::abs(1.0-sAold/sA);
            iter += 1;
            if (iter>20) {
               success = -1;
               break;
            }

         }

         return success;

      };

      inline int setFAS_v2_coefficients( Real& alpha_g,
                                         Real& sA,
                                   const Real  Clog,
                                   const Real  s12,
                                   const Real  mu_max,
                                   const Real  mu_tr ) const
      {
         // Newton Routine to solve for alpha_g and sA for the following equations:
         // 1 = S_C + S_L (total probability = 1)
         // mu_N97 = mu_C + mu_L (total variance of mu=sin^2(theta/2) matches Nanbu 97)
         // S_C = alpha_g*(1-exp(-2*mu_max/sA))           (integral of P_C)
         // S_L = s12/4/Clog*(1/mu_tr - 1/(mu_max+mu_tr)) (integral of P_L)
         // mu_N97 = 1/2*(1-exp(-s12))
         // mu_C = 0.5*alpha_g*sA*(1 - (1+mu_max/sA)*exp(-2*mu_max/sA))
         // mu_L = s12/4/Clog*log((mu_max+mu_tr)/mu_tr)

         Real sAold, f1, f2;
         int success = 1;

         const Real S_L = s12/(4.0*Clog)*mu_max/mu_tr/(mu_max+mu_tr);
         const Real mu_L = s12/(2.0*Clog)*std::log((mu_max+mu_tr)/mu_tr); // actually 2X mu_L
         const Real mu_N97 = 1.0-std::exp(-s12); // actually 2X mu_N97

         if (mu_L>mu_N97 || S_L>1.0) {
            alpha_g = 0.0;
            sA = mu_tr;
            success = -1;
            return success;
         }

         // initial guess for alpha_g and sA
         // P_N97 = alpha_g/sA*sin(theta)*std::exp(-2*sin(theta/2)^2/sA)
         sA = s12/2.0;
         alpha_g = std::max(0.0,(1.0-S_L)/(1.0-std::exp(-2.0*mu_max/sA)));

         int iter = 0;
         Real error = 1.0;
         while (error>1.0e-4) {

            sAold = sA;

            f1 = mu_N97 - mu_L + alpha_g*mu_max*std::exp(-2.0*mu_max/sAold);
            f2 = 1.0 - S_L;
            sA = f1/f2;

            // update alpha_g
            alpha_g = std::max(0.0,(1.0-S_L)/(1.0-std::exp(-2.0*mu_max/sA)));

            // check error
            error = std::abs(1.0-sAold/sA);
            iter += 1;
            if (iter>20) {
               success = -1;
               break;
            }

         }

         return success;

      };

      inline Real getTransitionX_NANBU( const Real  Clog,
                                        const Real  s12,
                                        const Real  alpha_g,
                                        const Real  sA ) const
      {
         // Newton method solve xc^2*exp(-2xc) = C0 with C0 = s12/(8*Clog*sA*alpha_g)
         // xc = sin^2(theta_c/2)/sA. This is the solution where the probability
         // distribution function for cumulative small angle scattering used in
         // the NANBU model and single scattering described by Rutherford intersect
         //
         // NOTE:
         // 1) The physical solution must be >= 1
         // 2) There is no solution for C0 < exp(-2)

         Real xc = 2.0;
         const Real C0 = s12/(8.0*Clog*alpha_g*sA);
         if (C0 > std::exp(-2.0)) {
            xc = 1.0;
            return xc;
         }

         Real error = 1.0;
         int iter = 0;

         Real xold, y0, dy0dx;
         while (error>1.0e-4) {
             xold = xc;
             y0 = xc*xc*std::exp(-2.0*xc) - C0;
             dy0dx = 2.0*xc*(1.0 - xc)*std::exp(-2.0*xc);
             xc = xc - y0/dy0dx; // update xc

             // check error
             error = std::abs(1.0-xold/xc);
             iter += 1;
             if (iter>20) { break; }
         }

         if (sA*xc > 1.0) { xc = 1.0/sA; }

         return xc;

      };

      inline Real getTransitionX_TA(const Real  Clog) const
      {
         // Newton method solve Clog*xc^3*exp(-xc^2/2) = sqrt(pi/2)
         // with xc = 2*tan(theta_c/2)/sqrt(<theta_sq>), this is
         // the solution where the probability distribution function for cumulative
         // small angle scattering used in the TA model and single scattering
         // described by Rutherford intersect

         Real xc = 2.57; // solution at a_Clog = 2 used for initial guess
         Real error = 1.0;
         int iter = 0;

         Real xold, xcsq, y0, dy0dx;
         while (error>1.0e-4) {
             xold = xc;
             xcsq = xc*xc;
             y0 = Clog*xc*xcsq*std::exp(-xcsq/2.0) - std::sqrt(Constants::PI/2.0);
             dy0dx = Clog*(3.0-xcsq)*xcsq*std::exp(-xcsq/2.0);
             xc = xc - y0/dy0dx; // update xc

             // check error
             error = std::abs(1.0-xold/xc);
             iter += 1;
             if (iter>20) { break; }
         }

         return xc;

      };

      inline void largeAngleScatter_method0( Real&  costh,
                                             Real&  sinth,
                                       const Real  Clog,
                                       const Real  sinsqthmaxo2,
                                       const Real  s12 ) const
      {
         // check for large angle Coulomb scattering event
         // update costh/sinth if necessary
         // passed in costh and sinth are those determined
         // from cumulative scattering model

         //const Real C0 = 2.5676 + std::pow(Clog-2.0,0.4)/2.8; // tan(thc/2) = xc*s12
         //const Real theta_c = std::min(C0*std::sqrt(2.0*s12),Constants::PI);
         const Real Xc = getTransitionX_TA(Clog); // Xc = 2*tan(thc/2)/sqrt(<theta_sq>)
         const Real multiplier = 4.0/3.0; // 4/3 gives smoother transition for small Clog
         const Real theta_c = 2.0*std::atan(Xc*multiplier*std::sqrt(s12/2.0));

         const Real sinsqthco2 = std::pow(std::sin(theta_c/2.0),2.0);
         const Real s12_L = s12/(4.0*Clog)*(2.0/sinsqthco2-1.0/sinsqthmaxo2);
         const Real P_L = 1.0 - std::exp(-s12_L);
         if (MathUtils::rand()<P_L) {
            const Real costhc = 1.0 - 2.0*sinsqthco2;
            const Real costhmax = 1.0 - 2.0*sinsqthmaxo2;
            const Real R_L = MathUtils::rand();
            Real costh_L = 1.0;
            if (MathUtils::rand()<1.0/(2.0-sinsqthco2/sinsqthmaxo2)) {
               costh_L = 1.0 - R_L*(1.0 - costhc);
            }
            else {
               const Real numer = R_L*(costhc - costhmax) - costhc*(1.0-costhmax);
               const Real denom = R_L*(costhc - costhmax) - (1.0 - costhmax);
               costh_L = numer/denom;
            }
            if (costh_L < costh) {
               costh = costh_L;
               sinth = std::sqrt(1.0 - costh*costh);
            }
         }

      }
      
      inline void modEnergyPairwise( std::array<Real,3>&  a_betap1,
                                     std::array<Real,3>&  a_betap2,
                               const long double          a_wpmp1,
                               const long double          a_wpmp2,
                                    const double          a_mult,
                                          double&         a_Erel_cumm,
                                     long double&         a_deltaE ) const
      { 
         // method to add/subtract a small percent of energy from a 
         // pair of particles such that momentum is not disturbed
         //
         // This is done by treating it as an inelastic scattering event
         // with potential U = a_deltaE and zero scattering angle
         // a_deltaE > 0 ==> need to take energy away
         // a_deltaE < 0 ==> need to add energy

         int sign = 1;
         if (a_deltaE<0.0) { sign = -1; }

         const Real ux = a_betap1[0] - a_betap2[0];
         const Real uy = a_betap1[1] - a_betap2[1];
         const Real uz = a_betap1[2] - a_betap2[2];
         const long double usq = ux*ux + uy*uy + uz*uz;
#ifdef RELATIVISTIC_PARTICLES
         long double gbsq1 = 0.0;
         for (int n=0; n<3; n++) { gbsq1 += a_betap1[n]*a_betap1[n]; }
         long double gbsq2 = 0.0;
         for (int n=0; n<3; n++) { gbsq2 += a_betap2[n]*a_betap2[n]; }
         const long double gamma1 = std::sqrt(1.0 + gbsq1);
         const long double gamma2 = std::sqrt(1.0 + gbsq2);
         const long double E1 = a_wpmp1*gamma1;
         const long double E2 = a_wpmp2*gamma2;
         const long double Etot = E1 + E2;
         const long double pxtot = a_wpmp1*a_betap1[0] + a_wpmp2*a_betap2[0];
         const long double pytot = a_wpmp1*a_betap1[1] + a_wpmp2*a_betap2[1];
         const long double pztot = a_wpmp1*a_betap1[2] + a_wpmp2*a_betap2[2];
         const long double Ecm = std::sqrt(Etot*Etot - pxtot*pxtot - pytot*pytot - pztot*pztot);
         const long double Erel = Ecm - a_wpmp1 - a_wpmp2;
#else
         const long double muR = a_wpmp1*a_wpmp2/(a_wpmp1 + a_wpmp2);
         const long double Erel = muR/2.0*usq;
#endif
         if (Erel<=0.0) { return; }

         long double deltaE = sign*m_energy_fraction*a_mult*Erel;
         if (std::abs(deltaE) > std::abs(a_deltaE)) {
            deltaE = a_deltaE;
            a_deltaE = 0.0;
         }
         else { a_deltaE -= deltaE; }
         a_Erel_cumm += Erel - deltaE;

#ifdef RELATIVISTIC_PARTICLES
         const long double A = Etot - deltaE;
         const long double D = A*A + E2*E2 - E1*E1;
         const long double p2dotu = a_wpmp2*(a_betap2[0]*ux + a_betap2[1]*uy + a_betap2[2]*uz);
         const long double ptdotu = pxtot*ux + pytot*uy + pztot*uz;

         // compute coefficients for quadratic equation for alpha
         const long double a = A*A*usq - ptdotu*ptdotu;
         const long double b = D*ptdotu - 2*A*A*p2dotu;
         const long double c = A*A*E2*E2 - D*D/4.0;

         const long double root = b*b - 4.0*a*c;
         if (root < 0.0 || a == 0.0) { return; }
         const long double alpha = (-b + std::sqrt(root))/(2.0*a);

         // update particle velocities
         const long double ratio1 = alpha/a_wpmp1;
         const long double ratio2 = alpha/a_wpmp2;
         a_betap1[0] += ratio1*ux;
         a_betap1[1] += ratio1*uy;
         a_betap1[2] += ratio1*uz;
         a_betap2[0] -= ratio2*ux;
         a_betap2[1] -= ratio2*uy;
         a_betap2[2] -= ratio2*uz;
#else
         const long double uprime_over_u = std::sqrt(1.0 - deltaE/Erel);
         std::array<Real,3> deltaU;
         deltaU[0] = uprime_over_u*ux - ux;
         deltaU[1] = uprime_over_u*uy - uy;
         deltaU[2] = uprime_over_u*uz - uz;

         for (int n=0; n<3; n++) {
            a_betap1[n] += muR/a_wpmp1*deltaU[n];
            a_betap2[n] -= muR/a_wpmp2*deltaU[n];
         }
#endif
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
         if (a_Eafter<Eafter2) return; // can happen for particles where weights are off to machine precision
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

      inline void printCorrectionCount() const
      {

         Real count_avg_local = 0.0;
         const DisjointBoxLayout& grids = m_count.disjointBoxLayout();
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            Box cell_box = grids[dit];
            IntVect range;
            for (int dir=0; dir<SpaceDim; dir++) {
               range[dir] = cell_box.bigEnd(dir) - cell_box.smallEnd(dir);
            }
            const Real num_cells = range.product();
            const FArrayBox& this_count = m_count[dit];
            for (BoxIterator bit(cell_box); bit.ok(); ++bit) {
               count_avg_local += this_count(bit(),0);
            }
            count_avg_local = count_avg_local/num_cells;
         }

         Real count_avg_global = 0.0;
#ifdef CH_MPI
         MPI_Reduce( &count_avg_local,
                     &count_avg_global,
                     1,
                     MPI_DOUBLE,
                     MPI_SUM,
                     0,
                     MPI_COMM_WORLD );
#else
         count_avg_global = avg_local;
#endif

         if (procID()==0) {
            double count_avg = count_avg_global / (double)numProc();
            cout << "JRA: numProc() = " << numProc() << endl;
            cout << "JRA: count_avg = " << count_avg << endl;
         }

      }

};

#include "NamespaceFooter.H"

#endif

