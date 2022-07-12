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
#include <vector>
#include <array>
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

#include "ElectroMagneticFields.H"
#include "PicSpeciesBC.H"

#include "NamespaceHeader.H"

class PicSpecies
{

public:

  PicSpecies( ParmParse&   a_ppspc,
        const int          a_species,
        const string&      a_name,
        const MeshInterp&  a_meshInterp,
        //const RefCountedPtr<MeshInterp>&  a_meshInterp,
        //MeshInterp*  a_meshInterp,
        const DomainGrid&  a_mesh );
  
  virtual ~PicSpecies();

  void initialize( const CodeUnits&    a_units,
                   const Real          a_time,
                   const std::string&  a_restart_file_name );

  void initializeFromInputFile( const CodeUnits&    a_units );
  void initializeFromRestartFile( const Real          a_time,
                                  const std::string&  a_restart_file_name );

  virtual Real mass() const { return m_mass; }
  virtual int charge() const { return m_charge; }
  virtual Real Uint() const { return m_Uint; }
  virtual Real stableDt() const { return m_stable_dt; }

  virtual const std::string& name() const { return m_name; }
  
  virtual int numParticles() const { return m_data.numParticles(); }

  bool motion() const { return m_motion; }
  bool forces() const { return m_forces; }
  bool scatter() const { return m_scatter; }
  bool writeAll() const { return m_write_all_part_comps; }
  bool chargeConservingDeposit() const { return m_charge_conserving_deposit; }
  
  InterpType getInterpType() const {return m_interpJToGrid;}
  
  virtual bool borisInertia() const { return m_use_axisymmetric_boris; }
  virtual bool isSpecies( const string& name ) const;
  

  // return non-const reference to particle data containers
  virtual ParticleData<JustinsParticle>& partData() {return m_data;}
  virtual LevelData<BinFab<JustinsParticlePtr>>& partData_binfab() {return m_data_binfab_ptr;}
  
  // return const reference to particle data containers
  virtual const ParticleData<JustinsParticle>& partData() const {return m_data;}
  virtual const LevelData<BinFab<JustinsParticlePtr>>& partData_binfab() const {return m_data_binfab_ptr;}
  
  // ways to get/set number density member
  virtual void numberDensity( LevelData<FArrayBox>& a_rho );
  
  virtual const LevelData<FArrayBox>& getNumberDensity() const {return m_density;}
  virtual const LevelData<FArrayBox>& getNumberDensity( const bool&  setData )
  { 
    if(setData) setNumberDensity();
    return m_density;
  }
  
  virtual const LevelData<FArrayBox>& getNumberDensityFromBinFab()
  { 
    setNumberDensityFromBinFab();
    return m_density_binfab;
  }

  void setNumberDensity();
  void setNumberDensityFromBinFab();
  
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

  void accumulateMassMatrices( LevelData<EdgeDataBox>&    a_sigma_xx, 
                               LevelData<EdgeDataBox>&    a_sigma_xy,
                               LevelData<EdgeDataBox>&    a_sigma_xz,
                               LevelData<NodeFArrayBox>&  a_sigma_yx,
                               LevelData<NodeFArrayBox>&  a_sigma_yy,
                               LevelData<NodeFArrayBox>&  a_sigma_yz,
                               LevelData<NodeFArrayBox>&  a_sigma_zx,
                               LevelData<NodeFArrayBox>&  a_sigma_zy,
                               LevelData<NodeFArrayBox>&  a_sigma_zz,  
                               LevelData<EdgeDataBox>&    a_Jtilde,
                               LevelData<NodeFArrayBox>&  a_Jtildev,
                         const ElectroMagneticFields&     a_em_fields,
                         const Real                       a_dt ) const;
  
  const LevelData<NodeFArrayBox>& getCurrentDensity_virtual() const {return m_currentDensity_virtual;}
  
  void setChargeDensity();
  void setChargeDensityOnFaces();
  void setChargeDensityOnNodes();
  
  const LevelData<FArrayBox>& getChargeDensity() const {return m_chargeDensity;}
  const LevelData<FluxBox>& getChargeDensityOnFaces() const {return m_chargeDensity_faces;}
  const LevelData<NodeFArrayBox>& getChargeDensityOnNodes() const {return m_chargeDensity_nodes;}

  void interpolateFieldsToParticles( const ElectroMagneticFields&  a_em_fields );

  void interpolateFieldsToOutcasts( const ElectroMagneticFields&  a_em_fields );

  void addExternalFieldsToParticles( const ElectroMagneticFields&  a_em_fields );
  
  void repositionInflowParticles(); 

  void repositionOutcastsAndApplyForces( const ElectroMagneticFields&   a_em_fields,
                                         const Real  a_dt, 
                                         const bool  a_byHalfDt );

  void advancePositions( const Real  a_full_dt,
                         const bool  a_half_step = false );
  void advancePositions_2ndHalf();
  
  void applyBCs( const bool  a_intermediate_advance );

  void advanceVelocities( const Real  a_full_dt, 
                          const bool  a_half_step );
  void advanceVelocities_2ndHalf();
  
  void applyInertialForces( const Real  a_dt,
                            const bool  a_use_avg_velocity,
                            const bool  a_update_positions,
                            const bool  a_half_positions = false );

  void advanceParticles( const ElectroMagneticFields&  a_em_fields,
                         const Real                    a_dt,
                         const bool                    a_iter_order_swap );
  
  void advanceParticlesIteratively( const ElectroMagneticFields&  a_em_fields,
                                    const Real                    a_dt,
                                    const bool                    a_iter_order_swap,
                                    const int                     a_iter_max );
  
  void removeOutflowParticles();
  
  void createInflowParticles( const Real&  a_time,
                              const Real&  a_dt );

  void updateOldParticlePositions();
  void updateOldParticleVelocities();
  
  // compute stable time step based on particle CFL
  void setStableDt( const CodeUnits&  a_units );
  
  // bin the particle up for scattering
  void binTheParticles();
  
  //const PicSpecies& operator=( const PicSpecies& a_species );
  
  void inspectBinFab()
  {
     inspectBinFab(m_data_binfab_ptr);
  }
  
  Real max_wpdt( const CodeUnits&  a_units,
                 const Real&       a_dt );
  void globalMoments(std::vector<Real>&  a_global_moments) const;
  
private:

  int  m_species;  // species number identifier 
  Real m_mass;     // species mass / electron mass
  int  m_charge;   // species charge / electron charge
  Real m_Uint;     // species internal energy [eV]
  
  bool m_use_axisymmetric_boris; 
  int  m_axisymmetric_iter_max; 

  PicSpeciesBC* m_species_bc;  
  std::string m_BC_type_lower, m_BC_type_upper;

  Real m_cvac_norm; 
  Real m_volume_scale; 
  Real m_stable_dt = DBL_MAX; 
  Real m_fnorm_const; // normalization constant for particle force

  std::string m_name;          // species name
  InterpType m_interpRhoToGrid;  // method to interpolate charge density to grid
  InterpType m_interpJToGrid;    // method to interpolate current density to grid
  InterpType m_interpEToParts;   // method to interpolate electric field to particle

  bool m_motion;
  bool m_forces;
  bool m_scatter;
  bool m_write_all_part_comps;
  bool m_charge_conserving_deposit;

  const DomainGrid& m_mesh;
  const MeshInterp& m_meshInterp;
  
  //MeshInterp* m_meshInterp;
  
  ParticleData<JustinsParticle> m_data; // container for storing list of particles for each box
  LevelData<BinFab<JustinsParticlePtr>> m_data_binfab_ptr; // container for storing array of lists of ptrs to JustinsParticles 
                                         // each grid cell has a list of all particles living in that cell

  mutable LevelData<FArrayBox> m_density;     // density moment at cell center
  mutable LevelData<FArrayBox> m_density_binfab;   // density moment from binfab
  mutable LevelData<FArrayBox> m_momentum;    // mean momentum density at cell center
  mutable LevelData<FArrayBox> m_energy;      // mean energy density at cell center
  
  // define containers for current density
  mutable LevelData<EdgeDataBox> m_currentDensity;          
  mutable LevelData<NodeFArrayBox> m_currentDensity_virtual;  // needed for 1D/2D simulations
  //mutable BoxLayoutData<FArrayBox> m_testBLD;          
  mutable LevelData<FArrayBox> m_chargeDensity;          
  mutable LevelData<FluxBox> m_chargeDensity_faces;          
  mutable LevelData<NodeFArrayBox> m_chargeDensity_nodes;          

  // physical quantities
  mutable LevelData<FArrayBox> m_temperature; // temperature [eV] at cell center
  mutable LevelData<FArrayBox> m_velocity;    // velocity [m/s] at cell center

  void advancePositions( List<JustinsParticle>&  a_pList,
                   const Real   a_cnormDt );

  void stepNormParticleTransfer( List<JustinsParticle>&  a_pList,
                                 List<JustinsParticle>&  a_temp_pList, 
                           const Real                    a_cnormHalfDt, 
                           const bool                    a_reverse );

  void applyForces( List<JustinsParticle>&  a_pList,
              const Real   a_cnormDt, 
              const bool   a_byHalfDt );
  
  void interpolateFieldsToParticles( List<JustinsParticle>&  a_pList,
                               const EdgeDataBox&  a_Efield_inPlane,
                               const FluxBox&      a_Bfield_inPlane );
  
  void interpolateFieldsToParticles( List<JustinsParticle>&  a_pList,
                               const EdgeDataBox&  a_Efield_inPlane,
                               const FluxBox&      a_Bfield_inPlane,
                               const FArrayBox&    a_Efield_virtual,
                               const FArrayBox&    a_Bfield_virtual );

  void newtonParticleUpdate( List<JustinsParticle>&  a_pList,
                       const Real          a_cnormDt, 
                       const EdgeDataBox&  a_Efield_inPlane,
                       const FluxBox&      a_Bfield_inPlane,
                       const FArrayBox&    a_Efield_virtual,
                       const FArrayBox&    a_Bfield_virtual );
   
  void cellCrossingPositionUpdate( List<JustinsParticle>&  a_pList, const Real  a_cnormDt );
  
  void addExternalFieldsToParticles( List<JustinsParticle>&  a_pList, 
                               const ElectroMagneticFields&  a_em_fields );
 
  void createMeshInterp();

  void inspectBinFab(const LevelData<BinFab<JustinsParticlePtr>>&  a_binfab_ptr);
  
  Real minimumParticleWeight( const LevelData<FArrayBox>&  a_density,
                              const RealVect&              a_sXmin,
                              const RealVect&              a_sXmax,
                              const Real                   a_numDen_scale,
                              const Real                   a_cell_volume,
                              const int                    a_partsPerCell ) const;

};

typedef RefCountedPtr<PicSpecies> PicSpeciesPtr;
typedef Vector<PicSpeciesPtr> PicSpeciesPtrVect;

#include "NamespaceFooter.H"

#endif
