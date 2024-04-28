/***
 *
 * PicSpeciesInterface class header file
 *
 * This class owns a pointer to all pic species. 
 * An instance of this class is owned by System.
 *
***/

#ifndef PicSpeciesInterface_h
#define PicSpeciesInterface_h

#include "PicSpecies.H"
#include "CodeUnits.H"
#include "DomainGrid.H"

#include "EMFields.H"

#include "NamespaceHeader.H"

class PicSpeciesInterface
{

public:

  PicSpeciesInterface( const CodeUnits&          a_units,
                       const DomainGrid&         a_mesh );
  
  virtual ~PicSpeciesInterface() {;}
     
  void initialize( const CodeUnits&    a_units,
                   const bool          a_implicit_advance,
                   const Real          a_cur_time,
                   const std::string&  a_restart_file_name );
 
  void initializeMassMatrices( EMFields&     a_emfields,
		         const Real&         a_dt,
		         const std::string&  a_restart_file_name );
  
  Real courantDt();
  Real cyclotronDt( const Real  a_wc0 );

  int numSpecies() const {return m_num_species;}
  int writeSpeciesRho() const {return m_writeSpeciesChargeDensity;}
  int writeSpeciesSigma() const {return m_writeSpeciesSurfaceCharge;}
  int writeSpeciesJ() const {return m_writeSpeciesCurrentDensity;}
  int writeSpeciesNppc() const {return m_writeSpeciesNppc;}
  int writeSpeciesEnergyOffDiag() const {return m_writeSpeciesEnergyOffDiag;}
  int writeSpeciesEnergyFlux() const {return m_writeSpeciesEnergyFlux;}

  PicSpeciesPtrVect& getPtrVect() {return m_pic_species_ptr_vect;}; 
  const PicSpeciesPtrVect& getPtrVect() const {return m_pic_species_ptr_vect;}; 
 
  //void setChargeDensity();
  //void setChargeDensityOnFaces();
  void setChargeDensityOnNodes( const bool  a_use_filtering );
  void setSurfaceChargeOnNodes( const bool  a_use_filtering );

  const LevelData<FArrayBox>& getChargeDensity() const {return m_chargeDensity;}
  const LevelData<FluxBox>& getChargeDensityOnFaces() const {return m_chargeDensity_faces;}
  const LevelData<NodeFArrayBox>& getChargeDensityOnNodes() const {return m_chargeDensity_nodes;}
  LevelData<NodeFArrayBox>& getChargeDensityOnNodes() {return m_chargeDensity_nodes;}
  
  const LevelData<NodeFArrayBox>& getSurfaceChargeOnNodes() const {return m_surfaceCharge_nodes;}
  LevelData<NodeFArrayBox>& getSurfaceChargeOnNodes() {return m_surfaceCharge_nodes;}
  
  void setCurrentDensity( const EMFields&  a_emfields,
                          const Real       a_dt,
		          const bool       a_finalizeJ,
		          const bool       a_from_explicit_solver = false );

  void filterJ( const EMFields&  a_emfields,
                const Real       a_time );

  const LevelData<EdgeDataBox>& getCurrentDensity() const {return m_currentDensity;}
  const LevelData<NodeFArrayBox>& getVirtualCurrentDensity() const {return m_currentDensity_virtual;}
  
  const LevelData<FArrayBox>& getDebyeLength() const {return m_DebyeLength;}

  void preRHSOp( const bool       a_from_emjacobian, 
                 const EMFields&  a_emfields,
                 const Real       a_dt,
                 const int        a_nonlinear_iter );

  void postNewtonUpdate( const EMFields&  a_emfields,
                         const Real       a_time,
                         const Real       a_dt );
 
  void setMassMatrices( const EMFields&  a_emfields,
                        const Real       a_time );
  
  void computeJfromMassMatrices( const EMFields&  a_emfields );

  void finalizeSettingJ( LevelData<EdgeDataBox>&    a_currentDensity,
                         LevelData<NodeFArrayBox>&  a_currentDensity_virtual,
                   const EMFields&                  a_emfields );
  
  void setMassMatricesForPC( EMFields&  a_emfields );
  
  void addInflowJ( LevelData<EdgeDataBox>&    a_J,
                   LevelData<NodeFArrayBox>&  a_Jv,
             const EMFields&                  a_emfields,
             const Real                       a_dt,
             const bool                       a_from_emjacobian ); 
  
  void addSubOrbitJ( LevelData<EdgeDataBox>&   a_J,
                    LevelData<NodeFArrayBox>&  a_Jv,
              const EMFields&                  a_emfields,
              const Real                       a_dt,
              const bool                       a_from_emjacobian ); 
 
  void prepForScatter( const int   a_num_coulomb,
		       const bool  a_from_scatterDt );
  void setDebyeLength();
  
  void postTimeStep() {
     for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
        PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
        species->removeOutflowParticles();
#if CH_SPACEDIM<3
	// this is done here in order correctly set the virtual
	// position of particles created after the pic advance,
	// such as in scatterParticles();
        species->rebaseVirtualPositions();
#endif
     }
  }
  
  bool forces() const {
     bool forces = false; 
     for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
        PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
        if( (species->charge()!=0) && (species->forces()) ) {
           forces = true;
           break;
        }
     }
     return forces;
  }
  
  InterpType getInterpForMassMatrix() const { 
     InterpType sp_interp = UNKNOWN;
     for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
        PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
        if(species->charge() == 0) continue;
        const InterpType& this_interp = species->getInterpType();
        if(sp_interp==UNKNOWN) sp_interp = this_interp;
        else CH_assert(sp_interp==this_interp);
     }
     return sp_interp;
  }
  
  bool useMassMatrices() const { return m_use_mass_matrices; }

  // const references to mass matrix containers
  const LevelData<EdgeDataBox>& getJ0() const {return m_J0;}
  const LevelData<EdgeDataBox>& getSigmaxx() const {return m_sigma_xx;}
  const LevelData<EdgeDataBox>& getSigmaxy() const {return m_sigma_xy;}
  const LevelData<EdgeDataBox>& getSigmaxz() const {return m_sigma_xz;}
  const LevelData<NodeFArrayBox>& getVirtualJ0() const {return m_J0_virtual;}
#if CH_SPACEDIM==1
  const LevelData<NodeFArrayBox>& getSigmayx() const {return m_sigma_yx;}
  const LevelData<NodeFArrayBox>& getSigmayy() const {return m_sigma_yy;}
  const LevelData<NodeFArrayBox>& getSigmayz() const {return m_sigma_yz;}
#endif
  const LevelData<NodeFArrayBox>& getSigmazx() const {return m_sigma_zx;}
  const LevelData<NodeFArrayBox>& getSigmazy() const {return m_sigma_zy;}
  const LevelData<NodeFArrayBox>& getSigmazz() const {return m_sigma_zz;}

#ifdef MASS_MATRIX_TEST
  const LevelData<EdgeDataBox>& getJtest() const {return m_currentDensity_TEST;}
  const LevelData<NodeFArrayBox>& getVirtualJtest() const {return m_currentDensity_virtual_TEST;}
#endif

  // non-const references to mass matrix containers
  LevelData<EdgeDataBox>& getJ0() {return m_J0;}
  LevelData<EdgeDataBox>& getSigmaxx() {return m_sigma_xx;}
  LevelData<EdgeDataBox>& getSigmaxy() {return m_sigma_xy;}
  LevelData<EdgeDataBox>& getSigmaxz() {return m_sigma_xz;}
  LevelData<NodeFArrayBox>& getVirtualJ0() {return m_J0_virtual;}
#if CH_SPACEDIM==1
  LevelData<NodeFArrayBox>& getSigmayx() {return m_sigma_yx;}
  LevelData<NodeFArrayBox>& getSigmayy() {return m_sigma_yy;}
  LevelData<NodeFArrayBox>& getSigmayz() {return m_sigma_yz;}
#endif
  LevelData<NodeFArrayBox>& getSigmazx() {return m_sigma_zx;}
  LevelData<NodeFArrayBox>& getSigmazy() {return m_sigma_zy;}
  LevelData<NodeFArrayBox>& getSigmazz() {return m_sigma_zz;}
  
  const IntVect& getNcompxx() const {return m_ncomp_xx;}
  const IntVect& getNcompxy() const {return m_ncomp_xy;}
  const IntVect& getNcompxz() const {return m_ncomp_xz;}
  const IntVect& getNcompyx() const {return m_ncomp_yx;}
  const IntVect& getNcompyy() const {return m_ncomp_yy;}
  const IntVect& getNcompyz() const {return m_ncomp_yz;}
  const IntVect& getNcompzx() const {return m_ncomp_zx;}
  const IntVect& getNcompzy() const {return m_ncomp_zy;}
  const IntVect& getNcompzz() const {return m_ncomp_zz;}
  
private:

  void createAllPicSpecies( const DomainGrid& a_mesh );
  void speciesPairingInit();
  
  bool m_verbosity;    
  
  bool m_writeSpeciesChargeDensity;
  bool m_writeSpeciesSurfaceCharge;
  bool m_writeSpeciesCurrentDensity;
  bool m_writeSpeciesNppc;
  bool m_writeSpeciesEnergyOffDiag;
  bool m_writeSpeciesEnergyFlux;
      
  bool m_verbose_particles;
  bool m_part_order_swap;
  int m_iter_max_particles;
  Real m_rtol_particles;
  int m_newton_maxits;
  int m_newton_num_guess;
      
  bool m_freeze_particles_jacobian;
  bool m_quasi_freeze_particles_jacobian;
  bool m_use_mass_matrices;
  bool m_mod_init_advance;
  
  int m_num_species;

  const DomainGrid& m_mesh;
  Real m_courant_dt;

  Real m_cvac_norm;

  PicSpeciesPtrVect m_pic_species_ptr_vect;

  mutable LevelData<FArrayBox> m_chargeDensity;
  mutable LevelData<FluxBox> m_chargeDensity_faces;
  mutable LevelData<NodeFArrayBox> m_chargeDensity_nodes;
  mutable LevelData<NodeFArrayBox> m_surfaceCharge_nodes;
  
  LevelData<EdgeDataBox> m_currentDensity;
  LevelData<NodeFArrayBox> m_currentDensity_virtual;

#ifdef MASS_MATRIX_TEST
  LevelData<EdgeDataBox> m_currentDensity_TEST;          
  LevelData<NodeFArrayBox> m_currentDensity_virtual_TEST;
#endif
  
  LevelData<FArrayBox> m_DebyeLength;
  
  // containers for semi-implicit mass matrix
  LevelData<EdgeDataBox> m_J0, m_E0;
  LevelData<EdgeDataBox> m_sigma_xx, m_sigma_xy, m_sigma_xz;
  LevelData<NodeFArrayBox> m_J0_virtual, m_E0_virtual;
#if CH_SPACEDIM==1
  LevelData<NodeFArrayBox> m_sigma_yx, m_sigma_yy, m_sigma_yz;
#endif
  LevelData<NodeFArrayBox> m_sigma_zx, m_sigma_zy, m_sigma_zz;
  IntVect m_ncomp_xx, m_ncomp_xy, m_ncomp_xz;
  IntVect m_ncomp_yx, m_ncomp_yy, m_ncomp_yz;
  IntVect m_ncomp_zx, m_ncomp_zy, m_ncomp_zz;

//#define PRINT_COMPS
#ifdef PRINT_COMPS
  bool m_print_comps;
#endif

};

#include "NamespaceFooter.H"

#endif
