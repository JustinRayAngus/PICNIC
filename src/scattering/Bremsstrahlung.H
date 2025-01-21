#ifndef _BREMSSTRAHLUNG_H_
#define _BREMSSTRAHLUNG_H_

#include "Scattering.H"
#include "PicnicConstants.H"
#include "JustinsParticle.H"
#include "ParticleData.H"

#include <cmath>

#include "NamespaceHeader.H"

/**
 * Bremsstrahlung scattering class:
 * e + H  ==> e + H + photon
 * e + He ==> e + He + photon
*/

class Bremsstrahlung : public Scattering
{
   public:

      Bremsstrahlung( ParmParse&  a_pp,
        const int&        a_verbosity )
         : Scattering(a_verbosity)
      {
         m_scatter_dt = DBL_MAX;
         m_deltaE_bremsstrahlung = 0.0;
         parseParameters( a_pp );
      }

      virtual ~Bremsstrahlung() = default;

      virtual ScatteringType getScatteringType() const { return BREMSSTRAHLUNG; }

      virtual std::string getScatteringSubTypeName() const {
          if      (m_bremsstrahlung_type==eH)  { return "eH"; }
          else if (m_bremsstrahlung_type==eHe) { return "eHe"; }
          else if (m_bremsstrahlung_type==eB)  { return "eB"; }
          else if (m_bremsstrahlung_type==eC)  { return "eC"; }
          else { return "invalid"; }
      }

      virtual int species1() const { return m_sp1; }

      virtual int species2() const { return m_sp2; }

      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const;

      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh );

      virtual Real scatterDt() const { return m_scatter_dt; }

      virtual Real getDeltaEbremsstrahlung() const { return m_deltaE_bremsstrahlung; }
      virtual void zeroDeltaEbremsstrahlung() const { m_deltaE_bremsstrahlung=0.0; }

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Bremsstrahlung scattering parameters:" << std::endl;
            std::cout << "  species 1 name:   " << m_species1_name << std::endl;
            std::cout << "  species 2 name:   " << m_species2_name << std::endl;
            std::cout << "  species 3 name:   " << m_species3_name << std::endl;
            std::cout << "  create_photons:   " << (m_create_photons?"true":"false") << std::endl;
            std::cout << "  ion atomic number = " << m_atomic_number << std::endl;
            std::cout << "  fmulti = " << m_fmulti << std::endl;
            std::cout << std::endl;
         }
      }

      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const;

   private:

      void setInterMFT( const LevelData<FArrayBox>&  a_numberDensity1,
                        const LevelData<FArrayBox>&  a_energyDensity1,
                        const LevelData<FArrayBox>&  a_numberDensity2,
                        const LevelData<FArrayBox>&  a_energyDensity2 ) const;

      void readCrossSectionFile();

      void interSpeciesBremsstrahlung( PicChargedSpecies&  a_picSpecies1,
                               PicChargedSpecies&  a_picSpecies2,
                               PicPhotonSpecies&  a_picSpecies3,
                         const DomainGrid&  a_mesh,
                         const Real         a_dt_sec ) const;

      enum XSEC_TYPE {
         MARTINEZ,
         TEXT_FILE
      };

      enum BREMSSTRAHLUNG_TYPE {
         eH,  // e + H  ==> e + H + photon
         eHe, // e + He ==> e + He + photon
         eB,  // e + B  ==> e + B + photon
         eC   // e + B  ==> e + C + photon
      };

      Real getSigma( const Real  KEcut_eV,
                     const Real  betarel,
                     const Real  KErel_eV ) const;
      Real getSigma_SB( const Real  KEcut_eV,
                        const Real  betarel,
                        const Real  KErel_eV ) const;

      // Parse the input database for parameters.
      void parseParameters( ParmParse&  a_pp )
      {

          // get number associated with each species
          a_pp.get( "species1", m_sp1 );
          a_pp.get( "species2", m_sp2 );
          a_pp.get( "species3", m_sp3 );
          a_pp.get( "atomic_number", m_atomic_number );
          if (m_atomic_number==1)  { m_bremsstrahlung_type = eH; }
          else if (m_atomic_number==2)  { m_bremsstrahlung_type = eHe; }
          else if (m_atomic_number==5)  { m_bremsstrahlung_type = eB; }
          else if (m_atomic_number==6)  { m_bremsstrahlung_type = eC; }
          else {
              MayDay::Error( "Bremsstrahlung: invalid atomic_number." );
          }
          a_pp.query( "create_photons", m_create_photons );
          a_pp.query( "fmulti", m_fmulti );
#ifdef RELATIVISTIC_PARTICLES
#else
          MayDay::Error( "Bremsstrahlung scattering model requires RELATIVISTIC_PARTICLES" );
#endif

      }

      std::string m_species1_name;
      std::string m_species2_name;
      std::string m_species3_name;
      mutable Real m_scatter_dt;    // (mean free time)
      mutable Real m_deltaE_bremsstrahlung; // running sum of energy converted to brems
      int m_sp1;      // species 1 identifier (electron)
      int m_sp2;      // species 2 identifier (ion)
      int m_sp3;      // species 3 identifier (photon)
      Real m_mass1;   // species 1 mass / electron mass
      Real m_mass2;   // species 2 mass / electron mass
      Real m_mass3 = 0.0;   // photon mass / electron mass
      int m_atomic_number; // atomic number of ion species
      Real m_fmulti = 1.0; // production factor (for better statistics)
   
      // relative photon energy grid ( koT1 = Ephoton_eV/KErel_eV )
      const std::vector<Real> m_koT1_grid = {0., 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.99, 1.00};
      
      // energy grid for electrons in eV
      const std::vector<Real> m_KEgrid_eV = {1.0e3, 2.0e3, 5.0e3, 1.0e4, 2.0e4, 5.0e4, 1.0e5, 2.0e5, 5.0e5, 1.0e6, 2.0e6};
      
      // cumulative cross section on photon energy grid
      mutable std::vector<Real> m_k_eV, m_sigmaC;
      mutable std::vector<std::vector<Real>> m_kdsigdk;

      // scaled energy-weighted differential total cross section for Bremsstrahlung (kdsigdk_n + Z*kdsigdk_e)
      // for e + H from Table I of Seltzer and Berger, ATOMIC DATA AND NUCLEAR DATA TABLES 35,345-418 (1985).
      // The values here are actualy beta1^2/Z^2*k*dsig/dk in units of mBarns = 1e-31 m^2 where Z is the
      // atomic number and beta1 = sqrt(1.0 - 1/gamma1^2) with gamma1 = 1 + T1/(me*c^2). These vectors are
      // converted to k*dsig/dk in units of m^2 on initialization.

      // atomic number Z = 1
      std::vector<vector<Real>> m_kdsigdk_eH = {
          {7.853, 7.849, 7.833, 7.800, 7.746, 7.446, 7.040, 6.586, 6.124, 5.664, 5.230, 4.841, 4.521, 4.400, 4.372, 4.362, 4.360}, // 1.0 keV
          {8.805, 8.817, 8.801, 8.738, 8.638, 8.059, 7.377, 6.699, 6.052, 5.431, 4.839, 4.263, 3.682, 3.437, 3.374, 3.326, 3.320}, // 2.0 keV
          {10.32, 10.25, 10.15, 9.976, 9.753, 8.703, 7.661, 6.728, 5.899, 5.148, 4.424, 3.697, 2.897, 2.446, 2.286, 2.161, 2.104}, // 5.0 keV
          {11.55, 11.40, 11.19, 10.89, 10.55, 9.087, 7.795, 6.711, 5.776, 4.952, 4.172, 3.387, 2.509, 1.982, 1.762, 1.548, 1.439}, // 10 keV
          {13.07, 12.61, 12.13, 11.62, 11.11, 9.370, 7.915, 6.674, 5.617, 4.734, 3.933, 3.141, 2.239, 1.684, 1.424, 1.129, 0.967}, // 20 keV
          {14.90, 14.17, 13.42, 12.65, 11.92, 9.620, 7.858, 6.426, 5.295, 4.370, 3.575, 2.777, 1.940, 1.401, 1.112, 0.754, 0.551}, // 50 keV
          {16.94, 15.72, 14.55, 13.47, 12.51, 9.686, 7.689, 6.107, 4.937, 3.951, 3.179, 2.438, 1.670, 1.180, 0.903, 0.548, 0.343}, // 100 keV
          {20.31, 18.17, 16.25, 14.71, 13.46, 9.837, 7.509, 5.808, 4.595, 3.549, 2.694, 2.022, 1.362, 0.938, 0.697, 0.380, 0.197}, // 200 keV
          {27.90, 23.39, 19.68, 17.28, 15.68, 10.84, 7.831, 5.945, 4.555, 3.400, 2.386, 1.634, 1.045, 0.706, 0.505, 0.244, 0.094}, // 500 keV
          {35.57, 28.83, 23.45, 20.32, 18.41, 12.21, 8.960, 6.849, 5.241, 3.927, 2.763, 1.732, 1.038, 0.688, 0.483, 0.218, 0.065}, // 1.0 MeV
          {43.42, 34.55, 27.56, 23.66, 21.45, 14.61, 11.03, 8.605, 6.744, 5.195, 3.809, 2.472, 1.306, 0.846, 0.587, 0.250, 0.057}  // 2.0 MeV
      };

      // atomic number Z = 2
      std::vector<vector<Real>> m_kdsigdk_eHe = {
          {7.167, 7.192, 7.206, 7.201, 7.181, 7.001, 6.726, 6.409, 6.077, 5.740, 5.422, 5.150, 4.939, 4.858, 4.837, 4.825, 4.821}, // 1.0 keV
          {8.232, 8.239, 8.229, 8.187, 8.120, 7.713, 7.200, 6.666, 6.147, 5.646, 5.173, 4.729, 4.313, 4.145, 4.099, 4.064, 4.053}, // 2.0 keV
          {9.678, 9.640, 9.570, 9.444, 9.276, 8.439, 7.568, 6.770, 6.055, 5.397, 4.770, 4.135, 3.510, 3.180, 3.070, 2.986, 2.951}, // 5.0 keV
          {10.81, 10.71, 10.56, 10.33, 10.06, 8.834, 7.706, 6.747, 5.914, 5.166, 4.461, 3.762, 3.009, 2.590, 2.430, 2.283, 2.213}, // 10 keV
          {12.18, 11.80, 11.40, 10.98, 10.55, 9.062, 7.784, 6.678, 5.721, 4.897, 4.148, 3.420, 2.605, 2.131, 1.927, 1.711, 1.600}, // 20 keV
          {13.49, 12.92, 12.33, 11.71, 11.11, 9.139, 7.592, 6.336, 5.325, 4.473, 3.707, 2.950, 2.148, 1.657, 1.413, 1.130, 0.974}, // 50 keV
          {14.54, 13.71, 12.88, 12.05, 11.28, 8.928, 7.241, 5.917, 4.892, 4.017, 3.265, 2.544, 1.797, 1.334, 1.093, 0.792, 0.624}, // 100 keV
          {16.12, 14.79, 13.54, 12.44, 11.49, 8.747, 6.868, 5.470, 4.418, 3.513, 2.751, 2.092, 1.433, 1.023, 0.802, 0.525, 0.365}, // 200 keV
          {19.94, 17.44, 16.27, 13.69, 12.50, 9.013, 6.783, 5.288, 4.140, 3.184, 2.345, 1.667, 1.080, 0.745, 0.556, 0.315, 0.178}, // 500 keV
          {24.04, 20.39, 17.38, 15.42, 14.06, 9.865, 7.470, 5.820, 4.535, 3.479, 2.549, 1.714, 1.060, 0.713, 0.517, 0.266, 0.123}, // 1.0 MeV
          {28.11, 23.46, 19.71, 17.43, 15.98, 11.46, 8.864, 7.025, 5.582, 4.375, 3.303, 2.268, 1.320, 0.866, 0.613, 0.292, 0.108}  // 2.0 MeV
      };

      // atomic number Z = 5
      std::vector<vector<Real>> m_kdsigdk_eB = {
          {5.670, 5.717, 5.760, 5.793, 5.820, 5.871, 5.856, 5.797, 5.713, 5.610, 5.502, 5.411, 5.350, 5.323, 5.311, 5.299, 5.293}, // 1.0 keV
          {6.814, 6.863, 6.903, 6.925, 6.932, 6.863, 6.698, 6.484, 6.250, 6.010, 5.789, 5.596, 5.438, 5.366, 5.340, 5.318, 5.306}, // 2.0 keV
          {8.325, 8.350, 8.360, 8.324, 8.265, 7.870, 7.382, 6.895, 6.441, 6.020, 5.635, 5.285, 4.987, 4.880, 4.847, 4.814, 4.799}, // 5.0 keV
          {9.457, 9.432, 9.377, 9.272, 9.128, 8.384, 7.614, 6.920, 6.310, 5.759, 5.259, 4.795, 4.372, 4.201, 4.150, 4.111, 4.097}, // 10 keV
          {10.56, 10.44, 10.29, 10.07, 9.819, 8.684, 7.651, 6.784, 6.042, 5.381, 4.785, 4.222, 3.676, 3.414, 3.332, 3.275, 3.252}, // 20 keV
          {12.00, 11.56, 11.11, 10.67, 10.23, 8.718, 7.455, 6.396, 5.515, 4.765, 4.097, 3.463, 2.785, 2.431, 2.302, 2.185, 2.136}, // 50 keV
          {12.69, 12.07, 11.45, 10.86, 10.30, 8.452, 7.031, 5.897, 4.996, 4.223, 3.526, 2.867, 2.182, 1.807, 1.650, 1.490, 1.410}, // 100 keV
          {13.49, 12.57, 11.70, 10.91, 10.22, 8.039, 6.503, 5.341, 4.421, 3.625, 2.914, 2.286, 1.651, 1.287, 1.121, 0.940, 0.844}, // 200 keV
          {15.27, 13.74, 12.38, 11.33, 10.49, 7.918, 6.188, 4.945, 3.953, 3.121, 2.392, 1.765, 1.185, 0.865, 0.706, 0.521, 0.419}, // 500 keV
          {17.23, 15.27, 13.58, 12.36, 11.41, 8.400, 6.550, 5.235, 4.166, 3.260, 2.464, 1.749, 1.131, 0.788, 0.617, 0.409, 0.293}, // 1.0 MeV
          {19.34, 16.89, 14.86, 13.50, 12.53, 9.504, 7.545, 6.077, 4.899, 3.908, 3.027, 2.183, 1.368, 0.924, 0.693, 0.414, 0.256}  // 2.0 MeV
      };

      // atomic number Z = 6
      std::vector<vector<Real>> m_kdsigdk_eC = {
          {5.336, 5.384, 5.427, 5.464, 5.494, 5.570, 5.585, 5.557, 5.501, 5.425, 5.337, 5.259, 5.206, 5.185, 5.175, 5.164, 5.158}, // 1.0 keV
          {6.498, 6.553, 6.600, 6.630, 6.648, 6.628, 6.518, 6.359, 6.174, 5.976, 5.790, 5.619, 5.470, 5.400, 5.375, 5.355, 5.346}, // 2.0 keV
          {8.028, 8.065, 8.084, 8.070, 8.031, 7.721, 7.315, 6.897, 6.502, 6.135, 5.801, 5.503, 5.251, 5.160, 5.129, 5.096, 5.080}, // 5.0 keV
          {9.168, 9.159, 9.123, 9.041, 8.922, 8.281, 7.593, 6.964, 6.411, 5.913, 5.467, 5.064, 4.712, 4.577, 4.536, 4.500, 4.485}, // 10 keV
          {10.27, 10.18, 10.05, 9.865, 9.638, 8.606, 7.644, 6.832, 6.140, 5.522, 4.973, 4.465, 3.994, 3.779, 3.715, 3.671, 3.654}, // 20 keV
          {11.72, 11.31, 10.90, 10.49, 10.08, 8.651, 7.450, 6.436, 5.587, 4.862, 4.221, 3.622, 2.997, 2.684, 2.578, 2.489, 2.455}, // 50 keV
          {12.38, 11.81, 11.24, 10.68, 10.16, 8.400, 7.026, 5.924, 5.045, 4.293, 3.611, 2.971, 2.314, 1.967, 1.830, 1.701, 1.638}, // 100 keV
          {13.10, 12.26, 11.45, 10.72, 10.06, 7.971, 6.481, 5.348, 4.447, 3.668, 2.970, 2.351, 1.726, 1.378, 1.225, 1.068, 0.987}, // 200 keV
          {14.65, 13.27, 12.04, 11.06, 10.26, 7.806, 6.132, 4.921, 3.953, 3.133, 2.417, 1.797, 1.222, 0.904, 0.756, 0.586, 0.494}, // 500 keV
          {16.39, 14.64, 13.12, 11.99, 11.09, 8.243, 6.460, 5.179, 4.135, 3.248, 2.467, 1.769, 1.153, 0.814, 0.650, 0.453, 0.345}, // 1.0 MeV
          {18.14, 16.07, 14.31, 13.09, 12.18, 9.254, 7.367, 5.974, 4.837, 3.869, 3.000, 2.184, 1.385, 0.942, 0.720, 0.453, 0.301}  // 2.0 MeV
      };

      bool m_create_photons = false;

      BREMSSTRAHLUNG_TYPE m_bremsstrahlung_type;
      Real m_mcSq;

};

#include "NamespaceFooter.H"

#endif

