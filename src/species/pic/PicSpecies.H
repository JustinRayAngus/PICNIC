/***
 *
 * base particle-in-cell species class header file
 *
***/

#ifndef PicSpecies_h
#define PicSpecies_h

#include "ParmParse.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "DomainGrid.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"
#include "Vector.H"
#include "RefCountedPtr.H"

#include "BinFab.H"
#include "ParticleData.H"
#include "Particle.H"
#include "JustinsParticle.H"
#include "JustinsParticlePtr.H"
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "MeshInterp.H"
#include "CodeUnits.H"

#include "NamespaceHeader.H"

class PicSpecies
{

public:

  PicSpecies( ParmParse&         a_ppspc,
              const string&      a_name,
              const MeshInterp&  a_meshInterp,
              //const RefCountedPtr<MeshInterp>&  a_meshInterp,
              //MeshInterp*  a_meshInterp,
              const DomainGrid&  a_mesh );
  
  virtual ~PicSpecies();

  void initialize( const CodeUnits&  a_units );

  virtual Real mass() const { return m_mass; }
  virtual int charge() const { return m_charge; }
  virtual Real Uint() const { return m_Uint; }
  virtual Real stableDt() const { return m_stable_dt; }

  virtual const std::string& name() const { return m_name; }

  virtual bool motion() const { return m_motion; }
  virtual bool forces() const { return m_forces; }
  virtual bool scatter() const { return m_scatter; }
  virtual bool writeAll() const { return m_write_all_part_comps; }
  
  virtual bool isSpecies( const string& name ) const;

  // return non-const reference to particle data containers
  virtual ParticleData<JustinsParticle>& partData() {return m_data;}
  virtual LevelData<BinFab<JustinsParticlePtr>>& partData_binfab() {return m_data_binfab_ptr;}
  
  // return const reference to particle data containers
  virtual const ParticleData<JustinsParticle>& partData() const {return m_data;}
  virtual const LevelData<BinFab<JustinsParticlePtr>>& partData_binfab() const {return m_data_binfab_ptr;}
  
  // ways to get/set number density member
  //
  virtual void numberDensity( LevelData<FArrayBox>& a_rho );
  
  virtual const LevelData<FArrayBox>& getNumberDensity() const {return m_density;}
  virtual const LevelData<FArrayBox>& getNumberDensity( const bool&  setData )
  { 
    if(setData) setNumberDensity();
    return m_density;
  }
  
  void setNumberDensity();
  
  // ways to get/set momentum density member
  //
  virtual void momentumDensity( LevelData<FArrayBox>& a_mom );
  
  virtual const LevelData<FArrayBox>& getMomentumDensity() const {return m_momentum;}
  virtual const LevelData<FArrayBox>& getMomentumDensity( const bool&  setData )
  { 
    if(setData) setMomentumDensity();
    return m_momentum;
  }
  
  void setMomentumDensity();
  
  // ways to get/set energy density member
  //
  virtual void energyDensity( LevelData<FArrayBox>& a_ene );
  
  virtual const LevelData<FArrayBox>& getEnergyDensity() const {return m_energy;}
  virtual const LevelData<FArrayBox>& getEnergyDensity( const bool setData = false )
  { 
    if(setData) setEnergyDensity();
    return m_energy;
  }
  
  void setEnergyDensity();
  
  // ways to get/set current density member
  //
  void currentDensity( LevelData<FArrayBox>&  a_currentDensity );
  
  void setCurrentDensity();
  const LevelData<EdgeDataBox>& getCurrentDensity() const {return m_currentDensity;}
  const LevelData<EdgeDataBox>& getCurrentDensity( const bool&  setData )
  { 
    if(setData) setCurrentDensity();
    return m_currentDensity;
  }
  
  const LevelData<NodeFArrayBox>& getCurrentDensity_virtual() const {return m_currentDensity_virtual;}
  
  void setChargeDensity();
  void setChargeDensityOnFaces();
  const LevelData<FArrayBox>& getChargeDensity() const {return m_chargeDensity;}
  const LevelData<FluxBox>& getChargeDensityOnFaces() const {return m_chargeDensity_FB;}
 
  // advance particle positions
  void advancePositions( const Real&  a_cnormDt );
  void advancePositions_2ndHalf();
  
  // compute stable time step based on particle CFL
  void setStableDt( const CodeUnits&  a_units );
  
  // bin the particle up for scattering
  void binTheParticles();
  
  // test particle shuffling
  void testParticleShuffling( const Real&  a_dt );

  //const PicSpecies& operator=( const PicSpecies& a_species );
  
  void inspectBinFab()
  {
     inspectBinFab(m_data_binfab_ptr);
  }
  
private:
  
  Real m_mass;     // species mass / electron mass
  int  m_charge;   // species charge / electron charge
  Real m_Uint;     // species internal energy [eV]

  Real m_volume_scale; 
  Real m_stable_dt = DBL_MAX; 

  std::string m_name;          // species name
  InterpType m_interpToGrid;  // method to interpolate charge/current density to grid

  bool m_motion;
  bool m_forces;
  bool m_scatter;
  bool m_write_all_part_comps;

  const DomainGrid& m_mesh;
  const MeshInterp& m_meshInterp;
  
  //MeshInterp* m_meshInterp;
  
  //ParticleData<Particle> m_data; // container for storing list of particle data
  ParticleData<JustinsParticle> m_data; // container for storing list of particles for each box
  LevelData<BinFab<JustinsParticlePtr>> m_data_binfab_ptr; // container for storing array of lists of ptrs to JustinsParticles 
                                         // each grid cell has a list of all particles living in that cell

  mutable LevelData<FArrayBox> m_density;     // density moment at cell center
  mutable LevelData<FArrayBox> m_momentum;    // mean momentum density at cell center
  mutable LevelData<FArrayBox> m_energy;      // mean energy density at cell center
  
  // define containers for current density
  mutable LevelData<EdgeDataBox> m_currentDensity;          
  mutable LevelData<NodeFArrayBox> m_currentDensity_virtual;  // needed for 1D/2D simulations
  //mutable BoxLayoutData<FArrayBox> m_testBLD;          
  mutable LevelData<FArrayBox> m_chargeDensity;          
  mutable LevelData<FluxBox> m_chargeDensity_FB;          

  // physical quantities
  //
  mutable LevelData<FArrayBox> m_temperature; // temperature [eV] at cell center
  mutable LevelData<FArrayBox> m_velocity;    // velocity [m/s] at cell center


  void createMeshInterp();

  void inspectBinFab(const LevelData<BinFab<JustinsParticlePtr>>&  a_binfab_ptr);

};

typedef RefCountedPtr<PicSpecies> PicSpeciesPtr;
typedef Vector<PicSpeciesPtr> PicSpeciesPtrVect;

#include "NamespaceFooter.H"

#endif
