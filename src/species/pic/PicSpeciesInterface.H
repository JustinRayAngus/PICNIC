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

#include "ElectroMagneticFields.H"

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
 
  void initializeMassMatrices( const std::string&  a_restart_file_name );
  
  Real courantDt();

  int numSpecies() const {return m_num_species;}
  int writeSpeciesRho() const {return m_writeSpeciesChargeDensity;}
  int writeSpeciesJ() const {return m_writeSpeciesCurrentDensity;}

  PicSpeciesPtrVect& getPtrVect() {return m_pic_species_ptr_vect;}; 
  const PicSpeciesPtrVect& getPtrVect() const {return m_pic_species_ptr_vect;}; 
 
  //void setChargeDensity();
  //void setChargeDensityOnFaces();
  void setChargeDensityOnNodes();
  
  const LevelData<FArrayBox>& getChargeDensity() const {return m_chargeDensity;}
  const LevelData<FluxBox>& getChargeDensityOnFaces() const {return m_chargeDensity_faces;}
  const LevelData<NodeFArrayBox>& getChargeDensityOnNodes() const {return m_chargeDensity_nodes;}
  
  void setCurrentDensity( const bool  a_from_explicit_solver = false );
  const LevelData<EdgeDataBox>& getCurrentDensity() const {return m_currentDensity;}
  const LevelData<NodeFArrayBox>& getVirtualCurrentDensity() const {return m_currentDensity_virtual;}

  void preRHSOp( const bool                    a_from_emjacobian, 
                 const ElectroMagneticFields&  a_emfields,
                 const Real                    a_dt,
                 const int                     a_nonlinear_iter );

  void postNewtonUpdate( const ElectroMagneticFields&  a_emfields,
                         const Real                    a_time,
                         const Real                    a_dt );
 
  void setMassMatrices( const ElectroMagneticFields&  a_emfields,
                        const Real                    a_time );
  
  void computeJfromMassMatrices( const ElectroMagneticFields&  a_emfields );
  void redefineJtilde( const ElectroMagneticFields&  a_emfields );
  
  void addInflowJ( LevelData<EdgeDataBox>&    a_J,
                   LevelData<NodeFArrayBox>&  a_Jv,
             const ElectroMagneticFields&     a_emfields,
             const Real                       a_dt );
 
  void prepForScatter();
  
  void postTimeStep() {
     for (int sp=0; sp<m_pic_species_ptr_vect.size(); sp++) {
        PicSpeciesPtr species(m_pic_species_ptr_vect[sp]);
        species->removeOutflowParticles();
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
  const LevelData<EdgeDataBox>& getJtilde() const {return m_Jtilde;}
  const LevelData<EdgeDataBox>& getSigmaxx() const {return m_sigma_xx;}
  const LevelData<EdgeDataBox>& getSigmaxy() const {return m_sigma_xy;}
  const LevelData<EdgeDataBox>& getSigmaxz() const {return m_sigma_xz;}
  const LevelData<NodeFArrayBox>& getVirtualJtilde() const {return m_Jtilde_virtual;}
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
  LevelData<EdgeDataBox>& getJtilde() {return m_Jtilde;}
  LevelData<EdgeDataBox>& getSigmaxx() {return m_sigma_xx;}
  LevelData<EdgeDataBox>& getSigmaxy() {return m_sigma_xy;}
  LevelData<EdgeDataBox>& getSigmaxz() {return m_sigma_xz;}
  LevelData<NodeFArrayBox>& getVirtualJtilde() {return m_Jtilde_virtual;}
#if CH_SPACEDIM==1
  LevelData<NodeFArrayBox>& getSigmayx() {return m_sigma_yx;}
  LevelData<NodeFArrayBox>& getSigmayy() {return m_sigma_yy;}
  LevelData<NodeFArrayBox>& getSigmayz() {return m_sigma_yz;}
#endif
  LevelData<NodeFArrayBox>& getSigmazx() {return m_sigma_zx;}
  LevelData<NodeFArrayBox>& getSigmazy() {return m_sigma_zy;}
  LevelData<NodeFArrayBox>& getSigmazz() {return m_sigma_zz;}
  
private:

  void createAllPicSpecies( const DomainGrid& a_mesh );
  void speciesPairingInit();
  
  bool m_verbosity;    
  bool m_writeSpeciesChargeDensity;
  bool m_writeSpeciesCurrentDensity;
      
  bool m_part_order_swap;
  bool m_iter_min_two;
  int m_iter_max_particles;
  Real m_rtol_particles;
      
  bool m_freeze_particles_jacobian;
  bool m_quasi_freeze_particles_jacobian;
  bool m_use_mass_matrices;
  bool m_mod_init_advance;
  
  int m_num_species;
 
  const DomainGrid& m_mesh;
  Real m_courant_dt;

  Real m_cvac_norm;
  Real m_volume_scale;

  PicSpeciesPtrVect m_pic_species_ptr_vect;

  LevelData<FArrayBox> m_chargeDensity;
  LevelData<FluxBox> m_chargeDensity_faces;
  LevelData<NodeFArrayBox> m_chargeDensity_nodes;
  
  LevelData<EdgeDataBox> m_currentDensity;
  LevelData<NodeFArrayBox> m_currentDensity_virtual;

#ifdef MASS_MATRIX_TEST
  LevelData<EdgeDataBox> m_currentDensity_TEST;          
  LevelData<NodeFArrayBox> m_currentDensity_virtual_TEST;
#endif
  
  // containers for semi-implicit mass matrix
  LevelData<EdgeDataBox> m_Jtilde;
  LevelData<EdgeDataBox> m_sigma_xx, m_sigma_xy, m_sigma_xz;
  LevelData<NodeFArrayBox> m_Jtilde_virtual;
#if CH_SPACEDIM==1
  LevelData<NodeFArrayBox> m_sigma_yx, m_sigma_yy, m_sigma_yz;
#endif
  LevelData<NodeFArrayBox> m_sigma_zx, m_sigma_zy, m_sigma_zz;
  IntVect m_num_xx_comps, m_num_xy_comps, m_num_xz_comps;
  IntVect m_num_yx_comps, m_num_yy_comps, m_num_yz_comps;
  IntVect m_num_zx_comps, m_num_zy_comps, m_num_zz_comps;

};

#include "NamespaceFooter.H"

#endif
