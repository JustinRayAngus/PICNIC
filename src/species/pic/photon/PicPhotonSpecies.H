/***
 *
 * photon particle-in-cell species class header file
 *
***/

#ifndef PIC_PHOTON_SPECIES_H
#define PIC_PHOTON_SPECIES_H

#include "ParmParse.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DomainGrid.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"
#include "Vector.H"
#include <vector>
#include <array>
#include "RefCountedPtr.H"

#include "BinFab.H"
#include "ParticleData.H"
#include "Particle.H"
#include "PhotonParticle.H"
#include "PhotonParticlePtr.H"
#include "PicPhotonSpeciesBC.H"
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "MeshInterp.H"
#include "CodeUnits.H"
#include "Species.H"

#include "NamespaceHeader.H"

class PicPhotonSpecies : public Species
{

public:

  PicPhotonSpecies( ParmParse&   a_ppspc,
              const int          a_species,
              const string&      a_name,
              const DomainGrid&  a_mesh );
  
  virtual ~PicPhotonSpecies();

  virtual void initialize( const CodeUnits&    a_units,
                           const Real          a_time,
                           const std::string&  a_restart_file_name );

  virtual SpeciesType getSpeciesType() const { return PIC_PHOTON; }

  void initializeFromRestartFile( const Real          a_time,
                                  const std::string&  a_restart_file_name );

  virtual int numParticles() const { return m_data.numParticles(); }

  const RealVect& getMassOut_lo() const { return m_MassOut_lo; }
  const RealVect& getMassOut_hi() const { return m_MassOut_hi; }
  const RealVect& getMomXOut_lo() const { return m_MomXOut_lo; }
  const RealVect& getMomXOut_hi() const { return m_MomXOut_hi; }
  const RealVect& getMomYOut_lo() const { return m_MomYOut_lo; }
  const RealVect& getMomYOut_hi() const { return m_MomYOut_hi; }
  const RealVect& getMomZOut_lo() const { return m_MomZOut_lo; }
  const RealVect& getMomZOut_hi() const { return m_MomZOut_hi; }
  const RealVect& getEnergyOut_lo() const { return m_EnergyOut_lo; }
  const RealVect& getEnergyOut_hi() const { return m_EnergyOut_hi; }

  void advancePhotons( const Real  a_full_dt,
                       const LevelData<FArrayBox>&  a_Ne );

  void applyBCs( const bool  a_intermediate_advance,
                 const Real  a_time );

  void removeOutflowParticles();

  // references to photon particle data containers
  ParticleData<PhotonParticle>& partData() {return m_data;}
  LevelData<BinFab<PhotonParticlePtr>>& partData_binfab() {return m_data_binfab_ptr;}
  
  // return const reference to particle data containers
  const ParticleData<PhotonParticle>& partData() const {return m_data;}
  const LevelData<BinFab<PhotonParticlePtr>>& partData_binfab() const {return m_data_binfab_ptr;}
  
  // ways to get/set number density member
  const LevelData<FArrayBox>& getNumberDensity( const bool  setData = false ) const
  { 
    if (setData) { setNumberDensity(); }
    return m_density;
  }
  
  const LevelData<FArrayBox>& getNumberDensityFromBinFab( const bool  setData = false ) const
  { 
    if (setData) { setNumberDensityFromBinFab(); }
    return m_density;
  }

  // ways to get/set momentum density member
  const LevelData<FArrayBox>& getMomentumDensity( const bool  setData = false ) const
  { 
    if(setData) setMomentumDensity();
    return m_momentumDensity;
  }
  

  // ways to get/set number density member
  const LevelData<FArrayBox>& getEnergyDensity( const bool  setData = false ) const
  { 
    if (setData) { setEnergyDensity(); }
    return m_energyDensity;
  }

  void setNppc() const;
  void setNumberDensity() const;
  void setNumberDensityFromBinFab() const;
  void setEnergyDensity() const;
  void setMomentumDensity() const;
  
  const LevelData<FArrayBox>& getNppc() const {return m_Nppc;}
  
  // bin the particle up for scattering
  void binTheParticles();
  
  void globalMoments(std::vector<Real>&  a_global_moments) const;

  void bdryMoments(std::vector<Real>&  a_bdry_moments);
  
private:

  PicPhotonSpeciesBC* m_species_bc;
  MeshInterp* m_meshInterp;

  Real m_cvac_norm;
  
  void createMeshInterp();
  
  ParticleData<PhotonParticle> m_data; // container for storing list of particles for each box
  LevelData<BinFab<PhotonParticlePtr>> m_data_binfab_ptr; // cell-based container for storing lists of ptrs to particles 

  mutable LevelData<FArrayBox> m_Nppc;             // number of macro particles per cell
  mutable LevelData<FArrayBox> m_density;          // density moment at cell center
  mutable LevelData<FArrayBox> m_energyDensity;    // energy density moment at cell center
  mutable LevelData<FArrayBox> m_momentumDensity;  // momentum density moment at cell center

  // boundary probes used for diagnostics
  RealVect m_MassOut_lo;         // total weight of outflow particles at lo side
  RealVect m_MassOut_hi;         // total weight of outflow particles at hi side
  RealVect m_MomXOut_lo;         // total X-momentum of outflow particles at lo side
  RealVect m_MomXOut_hi;         // total X-momentum of outflow particles at hi side
  RealVect m_MomYOut_lo;         // total Y-momentum of outflow particles at lo side
  RealVect m_MomYOut_hi;         // total Y-momentum of outflow particles at hi side
  RealVect m_MomZOut_lo;         // total Z-momentum of outflow particles at lo side
  RealVect m_MomZOut_hi;         // total Z-momentum of outflow particles at hi side
  RealVect m_EnergyOut_lo;       // total energy of outflow particles at lo side
  RealVect m_EnergyOut_hi;       // total energy of outflow particles at hi side
  
};

typedef RefCountedPtr<PicPhotonSpecies> PicPhotonSpeciesPtr;
typedef Vector<PicPhotonSpeciesPtr> PicPhotonSpeciesPtrVect;

#include "NamespaceFooter.H"

#endif
