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

#include "EMFields.H"
#include "PicSpeciesBC.H"

#include "NamespaceHeader.H"

class PicSpecies
{

public:

  PicSpecies( ParmParse&   a_ppspc,
        const int          a_species,
        const string&      a_name,
        const DomainGrid&  a_mesh );
  
  virtual ~PicSpecies();

  void initialize( const CodeUnits&    a_units,
                   const Real          a_time,
                   const std::string&  a_restart_file_name );

  void initializeFromInputFile( const CodeUnits&    a_units );
  void perturbPositions();
  void initializeFromRestartFile( const Real          a_time,
                                  const std::string&  a_restart_file_name );

  virtual Real mass() const { return m_mass; }
  virtual int charge() const { return m_charge; }
  virtual Real Uint() const { return m_Uint; }
  virtual Real stableDt() const { return m_stable_dt; }
   
  bool suborbit_inflowJ() const { return m_suborbit_inflowJ; }
  bool use_suborbit_model() const { return m_use_suborbit_model; }
  bool suborbit_fast_particles() const { return m_suborbit_fast_particles; }
  
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
  
  const RealVect& getMassIn_lo() const { return m_MassIn_lo; }
  const RealVect& getMassIn_hi() const { return m_MassIn_hi; }
  const RealVect& getMomXIn_lo() const { return m_MomXIn_lo; }
  const RealVect& getMomXIn_hi() const { return m_MomXIn_hi; }
  const RealVect& getMomYIn_lo() const { return m_MomYIn_lo; }
  const RealVect& getMomYIn_hi() const { return m_MomYIn_hi; }
  const RealVect& getMomZIn_lo() const { return m_MomZIn_lo; }
  const RealVect& getMomZIn_hi() const { return m_MomZIn_hi; }
  const RealVect& getEnergyIn_lo() const { return m_EnergyIn_lo; }
  const RealVect& getEnergyIn_hi() const { return m_EnergyIn_hi; }

  virtual const std::string& name() const { return m_name; }
  
  virtual int numParticles() const { return m_data.numParticles(); }

  bool motion() const { return m_motion; }
  bool forces() const { return m_forces; }
  bool scatter() const { return m_scatter; }
  bool writeAll() const { return m_write_all_part_comps; }
  bool chargeConservingDeposit() const { return m_charge_conserving_deposit; }
  
  const InterpType& getInterpType() const {return m_interpJToGrid;}
  
  virtual bool planarPush() const { 
     if(m_push_type==PLANAR) return true;
     else return false; 
  }
  
  virtual bool borisInertia() const { 
     if(m_push_type==CYL_BORIS) return true;
     else return false; 
  }
  
  virtual bool hybridInertia() const { 
     if(m_push_type==CYL_HYB) return true;
     else return false; 
  }

  virtual bool isSpecies( const string& name ) const;  

  void inspectParticles( List<JustinsParticle>&  a_pList ) const;

  void inspectParticles( List<JustinsParticle>&  a_pList,
		   const bool          a_print_fields,
                   const Box&          a_cell_box,
                   const EdgeDataBox&  a_Efield_inPlane,
                   const FluxBox&      a_Bfield_inPlane,
                   const FArrayBox&    a_Efield_virtual,
                   const FArrayBox&    a_Bfield_virtual ) const;
  
  virtual void setParticleSolverParams( const bool  a_verbose_particles,
                                        const bool  a_iter_order_swap,
                                        const int   a_iter_max,
                                        const Real  a_rtol, 
                                        const int   a_newton_maxits, 
                                        const int   a_newton_num_guess) 
  { 
     m_verbose_particles = a_verbose_particles; 
     m_iter_order_swap = a_iter_order_swap; 
     m_iter_max = a_iter_max; 
     m_rtol = a_rtol; 
     m_newton_maxits = a_newton_maxits;
     m_newton_num_guess = a_newton_num_guess;
  }

  // return non-const reference to particle data containers
  ParticleData<JustinsParticle>& partData() {return m_data;}
  LevelData<BinFab<JustinsParticlePtr>>& partData_binfab() {return m_data_binfab_ptr;}
  
  // return const reference to particle data containers
  const ParticleData<JustinsParticle>& partData() const {return m_data;}
  const LevelData<BinFab<JustinsParticlePtr>>& partData_binfab() const {return m_data_binfab_ptr;}
  
  // ways to get/set number density member
  const LevelData<FArrayBox>& getNumberDensity( const bool  setData = false ) const
  { 
    if(setData) setNumberDensity();
    return m_density;
  }
  
  const LevelData<FArrayBox>& getNumberDensityFromBinFab( const bool  setData = false ) const
  { 
    if(setData) setNumberDensityFromBinFab();
    return m_density;
  }

  void setNppc() const;
  void setNumberDensity() const;
  void setNumberDensityFromBinFab() const;
  
  const LevelData<FArrayBox>& getNppc() const {return m_Nppc;}
  
  // ways to get/set momentum density member
  const LevelData<FArrayBox>& getMomentumDensity( const bool  setData = false ) const
  { 
    if(setData) setMomentumDensity();
    return m_momentum;
  }
  
  const LevelData<FArrayBox>& getMomentumDensityFromBinFab( const bool  setData = false ) const
  { 
    if(setData) setMomentumDensityFromBinFab();
    return m_momentum;
  }
  
  void setMomentumDensity() const;
  void setMomentumDensityFromBinFab() const;
  
  // ways to get/set energy density member
  const LevelData<FArrayBox>& getEnergyDensity( const bool  setData = false ) const
  { 
    if(setData) setEnergyDensity();
    return m_energy;
  }
  
  const LevelData<FArrayBox>& getEnergyDensityFromBinFab( const bool  setData = false ) const
  { 
    if(setData) setEnergyDensityFromBinFab();
    return m_energy;
  }
  
  void setEnergyDensity() const;
  void setEnergyDensityFromBinFab() const;
  
  // ways to get/set energyOffDiag member
  const LevelData<FArrayBox>& getEnergyOffDiag( const bool  setData = false ) const
  { 
    if(setData) setEnergyOffDiag();
    return m_energyOffDiag;
  }
  void setEnergyOffDiag() const;
  
  // ways to get/set energy density flux member
  const LevelData<FArrayBox>& getEnergyDensityFlux( const bool  setData = false ) const
  { 
    if(setData) setEnergyDensityFlux();
    return m_energyFlux;
  }
  void setEnergyDensityFlux() const;
  
  // ways to get/set current density member
  void currentDensity( LevelData<FArrayBox>&  a_currentDensity );
  
  void setCurrentDensity( const Real  a_dt,
                          const bool  a_from_explicit_solver = false );
  const LevelData<EdgeDataBox>& getCurrentDensity() const {return m_currentDensity;}
  const LevelData<NodeFArrayBox>& getCurrentDensity_virtual() const {return m_currentDensity_virtual;}
  LevelData<EdgeDataBox>& getCurrentDensity() {return m_currentDensity;}
  LevelData<NodeFArrayBox>& getCurrentDensity_virtual() {return m_currentDensity_virtual;}
  
  const LevelData<EdgeDataBox>& getInflowJ() const {return m_inflowJ;}
  const LevelData<NodeFArrayBox>& getInflowJ_virtual() const {return m_inflowJ_virtual;}
  
  void advanceSubOrbitParticlesAndSetJ( const EMFields&  a_emfields,
		                        const Real       a_dt,
                                        const bool       a_from_emjacobian );
  
  void advanceInflowParticlesAndSetJ( const EMFields&  a_emfields,
		                      const Real       a_dt,
                                      const bool       a_from_emjacobian );

  const LevelData<EdgeDataBox>& getSubOrbitJ() const {return m_suborbitJ;}
  const LevelData<NodeFArrayBox>& getSubOrbitJ_virtual() const {return m_suborbitJ_virtual;}

  void accumulateMassMatrices( LevelData<EdgeDataBox>&    a_sigma_xx, 
                               LevelData<EdgeDataBox>&    a_sigma_xy,
                               LevelData<EdgeDataBox>&    a_sigma_xz,
#if CH_SPACEDIM==1
                               LevelData<NodeFArrayBox>&  a_sigma_yx,
                               LevelData<NodeFArrayBox>&  a_sigma_yy,
                               LevelData<NodeFArrayBox>&  a_sigma_yz,
#endif
                               LevelData<NodeFArrayBox>&  a_sigma_zx,
                               LevelData<NodeFArrayBox>&  a_sigma_zy,
                               LevelData<NodeFArrayBox>&  a_sigma_zz,  
                               LevelData<EdgeDataBox>&    a_J0,
                               LevelData<NodeFArrayBox>&  a_J0v,
                         const EMFields&                  a_emfields,
                         const Real                       a_dt ) const;
  
  
  void setChargeDensity();
  void setChargeDensityOnFaces();
  void setChargeDensityOnNodes( const bool  a_use_filtering );

#if CH_SPACEDIM>1
  // these two functions are used when writing the surface charge for each species
  void preWriteSurfaceCharge() const {
    if(m_charge==0) { return; }
    LDaddNodeOp<NodeFArrayBox> addNodeOp;
    m_surfaceCharge_nodes.exchange( m_surfaceCharge_nodes.interval(), 
                                    m_mesh.reverseCopier(), addNodeOp );
    SpaceUtils::exchangeNodeFArrayBox( m_surfaceCharge_nodes );
  }
  
  // m_surfaceCharge_nodes is a cummulative diagnostic. After doing addOp exchange,
  // need to mask the shared locations and then zero the ghost cells
  void postWriteSurfaceCharge() const {
    if(m_charge==0) { return; }
    const LevelData<NodeFArrayBox>& node_mask = m_mesh.getNodeMask();
    const DisjointBoxLayout& grids(m_mesh.getDBL());
    for(DataIterator dit(grids); dit.ok(); ++dit) {
      const FArrayBox&  mask = node_mask[dit].getFab();
      FArrayBox& sigma = m_surfaceCharge_nodes[dit].getFab();
      sigma.mult(mask,mask.box(),0,0,1);

      const Box node_box = surroundingNodes(grids[dit]);
      FArrayBox sigma_temp(node_box,1);
      sigma_temp.copy(sigma,node_box);
      sigma.setVal(0.0);
      sigma.copy(sigma_temp,node_box);
    } 
  }
#endif
  
  const LevelData<FArrayBox>& getChargeDensity() const {return m_chargeDensity;}
  const LevelData<FluxBox>& getChargeDensityOnFaces() const {return m_chargeDensity_faces;}
  const LevelData<NodeFArrayBox>& getChargeDensityOnNodes() const {return m_chargeDensity_nodes;}
  const LevelData<NodeFArrayBox>& getSurfaceChargeOnNodes() const {return m_surfaceCharge_nodes;}
  
  void getSurfaceChargeOnNodes( LevelData<NodeFArrayBox>& a_sigma ) const
  {
   
     SpaceUtils::zero( a_sigma );
     if(m_charge == 0) return;

     const DisjointBoxLayout& grids(m_mesh.getDBL());
     for(DataIterator dit(grids); dit.ok(); ++dit) {
        const FArrayBox& member_sigma = m_surfaceCharge_nodes[dit].getFab();
        FArrayBox& passed_sigma = a_sigma[dit].getFab();
        passed_sigma.copy(member_sigma);
     }

     LDaddNodeOp<NodeFArrayBox> addNodeOp;
     a_sigma.exchange( a_sigma.interval(), 
                       m_mesh.reverseCopier(), addNodeOp );
     SpaceUtils::exchangeNodeFArrayBox( a_sigma);
     
  }
  
  void interpolateEfieldToParticles( const EMFields&  a_emfields );
  
  void interpolateFieldsToParticles( const EMFields&  a_emfields );

  void addExternalFieldsToParticles( const EMFields&  a_emfields );

  void advancePositionsExplicit( const Real  a_full_dt,
                                 const bool  a_half_step = false );
  
  void advancePositionsExplicit( List<JustinsParticle>&  a_pList,
                           const Real                    a_cnormDt );
  
  void advancePositionsImplicit( const Real  a_full_dt );
  
  void advancePositionsImplicit( List<JustinsParticle>&  a_pList,
                           const Real                    a_cnormDt );
 
  void advancePositions_CYL_CAR( List<JustinsParticle>&  a_pList,
                          const Real                    a_cnormDt );

  void stepNormTransfer( List<JustinsParticle>&  a_in_pList,
                         List<JustinsParticle>&  a_out_pList,
                   const Real                    a_cnormDt, 
                   const bool                    a_reverse );
  
  void stepNormTransfer_CYL_CAR( List<JustinsParticle>&  a_in_pList,
                                 List<JustinsParticle>&  a_out_pList,
                           const Real                    a_cnormDt, 
                           const bool                    a_reverse );
  
  void stepNormTransfer_SPH_CAR( List<JustinsParticle>&  a_in_pList,
                                 List<JustinsParticle>&  a_out_pList,
                           const Real                    a_cnormDt, 
                           const bool                    a_reverse );
  
  void transferFastParticles();

  void advancePositions_2ndHalf();
  
  void advancePositions_2ndHalf( List<JustinsParticle>&  a_pList );
  
  void applyBCs( const bool  a_intermediate_advance, 
	         const Real  a_time );

  void checkForNAN( List<JustinsParticle>&  a_pList,
	      const std::string&  a_string );
  
  void checkForNAN( const std::string&  a_string );

  void advanceVelocities( const Real  a_full_dt, 
                          const bool  a_half_step );
  
  void advanceVelocities_2ndHalf( const Real  a_dt = 0.0 );
  
  void advanceVelocities_2ndHalf( List<JustinsParticle>&  a_pList );

  void advanceVelocities_CYL_HYB_2ndHalf( const Real  a_dt,
		              List<JustinsParticle>&  a_pList );
  
  void advanceVelocities_SPH_HYB_2ndHalf( const Real  a_dt,
		              List<JustinsParticle>&  a_pList );

  void advanceParticles_CYL_CAR_2ndHalf( List<JustinsParticle>&  a_pList );
  
#if CH_SPACEDIM==1
  void advanceParticles_SPH_CAR_2ndHalf( List<JustinsParticle>&  a_pList );
  void advancePositions_SPH_CAR( List<JustinsParticle>&  a_pList,
                          const Real                    a_cnormDt );
#endif

#if CH_SPACEDIM<3
  void rebaseVirtualPositions( );
  void rebaseVirtualPositions( List<JustinsParticle>&  a_pList );
#endif

  void applyInertialForces( const Real  a_dt,
                            const bool  a_use_avg_velocity,
                            const bool  a_update_positions,
                            const bool  a_half_positions = false );

  void advanceParticles( const EMFields&  a_emfields,
                         const Real       a_dt );
  
  void advanceParticlesIteratively( const EMFields&  a_emfields,
                                    const Real       a_dt );
  
  void mergeSubOrbitParticles();
  
  void removeOutflowParticles();
  
  void createInflowParticles( const Real  a_time,
                              const Real  a_dt );
  void injectInflowParticles();

  void updateOldParticlePositions();
  void updateOldParticleVelocities();
  
  // compute stable time step based on particle CFL
  void setStableDt();
  
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
  
  void bdryMoments(std::vector<Real>&  a_bdry_moments);
  
  void picardParams(std::vector<Real>&  a_picard_params);
  Real avg_picard_its();
  Real max_avg_picard_its();
  
private:

  int  m_species;  // species number identifier 
  Real m_mass;     // species mass / electron mass
  int  m_charge;   // species charge / electron charge
  Real m_Uint;     // species internal energy [eV]
  
  enum PushType {
     PLANAR, 
     CYL_BORIS,  // uses Boris trick for cylindrical push
     CYL_CYL,    // uses Cyl for velocity and Cyl for position
     CYL_HYB,    // uses Car for velocity and Cyl for position
     CYL_CAR,    // uses Car for velocity and Car for position
     SPH_SPH,    // uses Sph for velocity and Sph for position
     SPH_HYB,    // uses Car for velocity and Sph for position
     SPH_CAR,    // uses Car for velocity and Car for position
  };
  
  PushType m_push_type;

  PicSpeciesBC* m_species_bc;
  std::string m_BC_type_lower, m_BC_type_upper;

  Real m_cvac_norm; 
  Real m_stable_dt = DBL_MAX; 
  Real m_fnorm_const; // normalization constant for particle force

  bool m_iter_order_swap = false; 
  int m_iter_max = 0; 
  bool m_verbose_particles = false; 
  Real m_rtol = 1.0e-12; 
  int m_newton_maxits = 20; 
  int m_newton_num_guess = 1; 
  uint64_t m_num_parts_its = 0;
  uint64_t m_num_apply_its = 0;

  // boundary probes used for diagnostics
  RealVect m_MassOut_lo;         // total mass of outflow particles at lo side 
  RealVect m_MassOut_hi;         // total mass of outflow particles at hi side
  RealVect m_MomXOut_lo;         // total X-momentum of outflow particles at lo side 
  RealVect m_MomXOut_hi;         // total X-momentum of outflow particles at hi side
  RealVect m_MomYOut_lo;         // total Y-momentum of outflow particles at lo side 
  RealVect m_MomYOut_hi;         // total Y-momentum of outflow particles at hi side
  RealVect m_MomZOut_lo;         // total Z-momentum of outflow particles at lo side 
  RealVect m_MomZOut_hi;         // total Z-momentum of outflow particles at hi side
  RealVect m_EnergyOut_lo;       // total energy of outflow particles at lo side 
  RealVect m_EnergyOut_hi;       // total energy of outflow particles at hi side
  //
  RealVect m_MassIn_lo;         // total mass of inflow particles at lo side 
  RealVect m_MassIn_hi;         // total mass of inflow particles at hi side
  RealVect m_MomXIn_lo;         // total X-momentum of inflow particles at lo side 
  RealVect m_MomXIn_hi;         // total X-momentum of inflow particles at hi side
  RealVect m_MomYIn_lo;         // total Y-momentum of inflow particles at lo side 
  RealVect m_MomYIn_hi;         // total Y-momentum of inflow particles at hi side
  RealVect m_MomZIn_lo;         // total Z-momentum of inflow particles at lo side 
  RealVect m_MomZIn_hi;         // total Z-momentum of inflow particles at hi side
  RealVect m_EnergyIn_lo;       // total energy of inflow particles at lo side 
  RealVect m_EnergyIn_hi;       // total energy of inflow particles at hi side

  std::string m_name;          // species name
  InterpType m_interpRhoToGrid;  // method to interpolate charge density to grid
  InterpType m_interpJToGrid;    // method to interpolate current density to grid
  InterpType m_interpEToParts;   // method to interpolate electric field to particle
  bool m_interp_bc_check;
  bool m_suborbit_inflowJ;
  bool m_use_suborbit_model;
  bool m_suborbit_fast_particles;
  bool m_suborbit_testing = false; 

  bool m_motion;
  bool m_forces;
  bool m_scatter;
  bool m_write_all_part_comps;
  bool m_charge_conserving_deposit;
  
  const DomainGrid& m_mesh;
  MeshInterp* m_meshInterp;
  
  ParticleData<JustinsParticle> m_data; // container for storing list of particles for each box
  LevelData<BinFab<JustinsParticlePtr>> m_data_binfab_ptr; // cell-based container for storing lists of ptrs to particles 

  ParticleData<JustinsParticle> m_data_suborbit; // container for temporarily storing particles that are
                                                 // A) difficult to converge with Picard (subcycling needed)
						 // B) CFL-violating (to avoid need for extended MM stencil)
						 // These particles will not contribute to the mass matrices.
						 // Instead, their contribution to J will be computed directly
						 // during the linear iterations of PS-JFNK. At end of time step, 
						 // the particles are transferred back to the main container (m_data).

  mutable LevelData<FArrayBox> m_Nppc;           // number of macro particles per cell
  mutable LevelData<FArrayBox> m_density;        // density moment at cell center
  mutable LevelData<FArrayBox> m_momentum;       // momentum density at cell center
  mutable LevelData<FArrayBox> m_energy;         // energy density (diagonal elements of energy tensory) at cell center
  mutable LevelData<FArrayBox> m_energyOffDiag;  // off-diagonal (XY, XZ, and YZ) elements of energy tensor
  mutable LevelData<FArrayBox> m_energyFlux;     // energy density flux at cell center
  
  // define containers for current density
  mutable LevelData<EdgeDataBox> m_currentDensity, m_inflowJ, m_suborbitJ;
  mutable LevelData<NodeFArrayBox> m_currentDensity_virtual, m_inflowJ_virtual, m_suborbitJ_virtual;
  //mutable BoxLayoutData<FArrayBox> m_testBLD;          
  mutable LevelData<FArrayBox> m_chargeDensity;          
  mutable LevelData<FluxBox> m_chargeDensity_faces;          
  mutable LevelData<NodeFArrayBox> m_chargeDensity_nodes;          
  mutable LevelData<NodeFArrayBox> m_surfaceCharge_nodes;          

  // physical quantities
  mutable LevelData<FArrayBox> m_temperature; // temperature [eV] at cell center
  mutable LevelData<FArrayBox> m_velocity;    // velocity [m/s] at cell center

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
   
  void interpolateFieldsToParticle( JustinsParticle&  a_particle,
                              const EdgeDataBox&      a_Efield_inPlane,
                              const FluxBox&          a_Bfield_inPlane,
                              const FArrayBox&        a_Efield_virtual,
                              const FArrayBox&        a_Bfield_virtual );
  
  void addExternalFieldsToParticles( List<JustinsParticle>&  a_pList, 
                               const EMFields&               a_emfields );
  
  void advanceSubOrbitParticlesAndSetJ( List<JustinsParticle>&  a_pList,
		                        EdgeDataBox&            a_local_J,
		                        FArrayBox&              a_local_Jv,
		                  const EMFields&     a_emfields,
		                  const EdgeDataBox&  a_Efield_inPlane,
		                  const FluxBox&      a_Bfield_inPlane,
		                  const FArrayBox&    a_Efield_virtual,
		                  const FArrayBox&    a_Bfield_virtual,
		                  const int           a_bdry_dir,
		                  const int           a_bdry_side,
		                  const Real          a_dt,
                                  const bool          a_from_emjacobian );

  void advanceInflowPartToBdry( RealVect&            a_xpold,
                                Real&                a_cnormDt_sub,
                          const Real                 a_cnormDt,
		          const std::array<Real,3>&  a_upold,
                          const int                  a_bdry_dir,
                          const int                  a_bdry_side );
 
  void createMeshInterp();

  void inspectBinFab(const LevelData<BinFab<JustinsParticlePtr>>&  a_binfab_ptr);
  
  Real minimumParticleWeight( const LevelData<FArrayBox>&  a_density,
                              const RealVect&              a_sXmin,
                              const RealVect&              a_sXmax,
                              const Real                   a_numDen_scale,
                              const Real                   a_cell_volume,
                              const int                    a_partsPerCell ) const;
  
  Real minimumJacobian( ) const;

  void zeroFields() {
    CH_TIME("PicSpecies::zeroFields()");
    const BoxLayout& BL = m_data.getBoxes();
    DataIterator dit(BL);
    for(dit.begin(); dit.ok(); ++dit) {
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      ListIterator<JustinsParticle> lit(pList);
      for(lit.begin(); lit.ok(); ++lit) {
        std::array<Real,3> zero_array = {0.0,0.0,0.0};
        lit().setElectricField(zero_array);
        lit().setMagneticField(zero_array);
      }
    }
  }

};

typedef RefCountedPtr<PicSpecies> PicSpeciesPtr;
typedef Vector<PicSpeciesPtr> PicSpeciesPtrVect;

#include "NamespaceFooter.H"

#endif
