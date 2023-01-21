#ifndef ElectroMagneticFields_h
#define ElectroMagneticFields_h

#include <array>

#include "ParmParse.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "DomainGrid.H"
#include "CodeUnits.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"
#include "Vector.H"
#include "RefCountedPtr.H"
#include "FieldBC.H"
#include "GridFunctionFactory.H"

#include "MeshInterp.H"
//#include "PicSpecies.H" // JRA, can't do this. Circular

#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "NodePoissonUtilities.H"

#include "BandedMatrix.H"
#include "GlobalDOF.H"

#include "SpaceUtils.H"

#define _EM_PICARD_B_FIELD_BLOCK_ 0
#define _EM_PICARD_E_FIELD_BLOCK_ 1

#include "NamespaceHeader.H"

typedef enum{b_only, e_only, e_and_b} EMVecType;

class ElectroMagneticFields
{

public:

  ElectroMagneticFields( ParmParse&,
                   const DomainGrid&,
                   const CodeUnits&,
                   const bool&,
                   const EMVecType& a = e_and_b);
  
  virtual ~ElectroMagneticFields();

  void initialize( const Real          a_cur_time,
                   const std::string&  a_restart_file_name );

  void initExternalFields();
  
  bool advance() const { return m_advance; }
  
  bool externalFields() const { return m_external_fields; }
  
  std::array<bool,3> advanceE_comp() const { return m_advanceE_comp; }
  std::array<bool,3> advanceB_comp() const { return m_advanceB_comp; }
  
  bool usePoisson() const { return m_use_poisson; }
  
  bool writeDivs() const { return m_writeDivs; }
  bool writeCurls() const { return m_writeCurls; }
  bool writeRho() const { return m_writeRho; }
  bool writeExB() const { return m_writeExB; }
  
  Real stableDt() const { return m_stable_dt; }
  
  const FieldBC& getFieldBC() const { return *m_field_bc; }  

  const LevelData<FluxBox>& getMagneticField() const {return m_magneticField;}
  LevelData<FluxBox>& getMagneticField() {return m_magneticField;}
  const LevelData<FArrayBox>& getVirtualMagneticField() const {return m_magneticField_virtual;}
  LevelData<FArrayBox>& getVirtualMagneticField() {return m_magneticField_virtual;}
 
  const LevelData<FluxBox>& getMagneticField_old() const {return m_magneticField_old;}
  LevelData<FluxBox>& getMagneticField_old() {return m_magneticField_old;}
  const LevelData<FArrayBox>& getVirtualMagneticField_old() const {return m_magneticField_virtual_old;}
  LevelData<FArrayBox>& getVirtualMagneticField_old() {return m_magneticField_virtual_old;}

  const LevelData<EdgeDataBox>& getElectricField() const {return m_electricField;}
  LevelData<EdgeDataBox>& getElectricField() {return m_electricField;}
  const LevelData<NodeFArrayBox>& getVirtualElectricField() const {return m_electricField_virtual;}
  LevelData<NodeFArrayBox>& getVirtualElectricField() {return m_electricField_virtual;}
  
  const LevelData<EdgeDataBox>& getElectricField_old() const {return m_electricField_old;}
  LevelData<EdgeDataBox>& getElectricField_old() {return m_electricField_old;}
  const LevelData<NodeFArrayBox>& getVirtualElectricField_old() const {return m_electricField_virtual_old;}
  LevelData<NodeFArrayBox>& getVirtualElectricField_old() {return m_electricField_virtual_old;}

  const LevelData<FluxBox>& getMagneticFieldRHS() const {return m_magneticField_rhs;}
  LevelData<FluxBox>& getMagneticFieldRHS() {return m_magneticField_rhs;}
  const LevelData<FArrayBox>& getVirtualMagneticFieldRHS() const {return m_magneticField_virtual_rhs;}
  LevelData<FArrayBox>& getVirtualMagneticFieldRHS() {return m_magneticField_virtual_rhs;}

  const LevelData<EdgeDataBox>& getElectricFieldRHS() const {return m_electricField_rhs;}
  LevelData<EdgeDataBox>& getElectricFieldRHS() {return m_electricField_rhs;}
  const LevelData<NodeFArrayBox>& getVirtualElectricFieldRHS() const {return m_electricField_virtual_rhs;}
  LevelData<NodeFArrayBox>& getVirtualElectricFieldRHS() {return m_electricField_virtual_rhs;}

  // const references to charge/current density
  const LevelData<EdgeDataBox>& getCurrentDensity() const {return m_currentDensity;}
  const LevelData<NodeFArrayBox>& getVirtualCurrentDensity() const {return m_currentDensity_virtual;}
  const LevelData<NodeFArrayBox>& getChargeDensity() const {return m_chargeDensity;}
  
  // non-const references to charge/current density
  LevelData<EdgeDataBox>& getCurrentDensity() {return m_currentDensity;}
  LevelData<NodeFArrayBox>& getVirtualCurrentDensity() {return m_currentDensity_virtual;}
  LevelData<NodeFArrayBox>& getChargeDensity() {return m_chargeDensity;}
  
  // const reference to m_rhs_vector and m_phi_vector zero comp
  const LevelData<NodeFArrayBox>& getPotential() const {return *m_phi_vector[0];}
  const LevelData<NodeFArrayBox>& getRHS() const {return *m_rhs_vector[0];}
  const LevelData<EdgeDataBox>& getElectricFieldCorrection() const {return m_electricField_correction;}

  // const references to curl of electric field  
  const LevelData<FluxBox>& getCurlE() const {return m_curlE;}
  const LevelData<FArrayBox>& getCurlE_virtual() const {return m_curlE_virtual;}
  
  // const references to curl of magnetic field  
  const LevelData<EdgeDataBox>& getCurlB() const {return m_curlB;}
  const LevelData<NodeFArrayBox>& getCurlB_virtual() const {return m_curlB_virtual;}
  
  // const references to divE and divB
  const LevelData<NodeFArrayBox>& getDivE() const {return m_divE;}
  const LevelData<FArrayBox>& getDivB() const {return m_divB;}

  // get external electric field values
  std::array<Real,3> getExternalE( const RealVect&  a_X ) const
  {  
     std::array<Real,3> extE = {0.0,0.0,0.0};
     if(m_external_fields) {
        extE[0] = m_gridFunction_extE0->getValue(a_X);
        extE[1] = m_gridFunction_extE1->getValue(a_X);
        extE[2] = m_gridFunction_extE2->getValue(a_X);
     }
     return extE;
  }
  
  // get external magnetic field values
  std::array<Real,3> getExternalB( const RealVect&  a_X ) const
  {  
     std::array<Real,3> extB = {0.0,0.0,0.0};
     if(m_external_fields) {
        extB[0] = m_gridFunction_extB0->getValue(a_X);
        extB[1] = m_gridFunction_extB1->getValue(a_X);
        extB[2] = m_gridFunction_extB2->getValue(a_X);
     }
     return extB;
  }
 
  // advance the fields in time
  void computeRHSElectricField( const Real, const Real );
  void computeRHSMagneticField( const Real, const Real );
  
  // compute the poynting flux 
  void setPoyntingFlux();
  
  // const references to components of ExB  
  const LevelData<EdgeDataBox>& getExB() const {return m_ExB;}
  const LevelData<NodeFArrayBox>& getEvxB() const {return m_EvxB;}
  
  // compute field energy
  Real electricFieldEnergy() const;
  Real magneticFieldEnergy() const;

  Real max_wc0dt( const CodeUnits&  a_units,
                  const Real&       a_dt ) const;
 
  const RealVect& getIntSdA_lo() const { return m_intSdA_lo; }  
  const RealVect& getIntSdA_hi() const { return m_intSdA_hi; } 
  const RealVect& getIntSdAdt_lo() const { return m_intSdAdt_lo; }  
  const RealVect& getIntSdAdt_hi() const { return m_intSdAdt_hi; }  

  // compute the divegence of the fields 
  void setDivB();
  void setDivE();

  void enforceGaussLaw( const Real  a_time );
  void solvePoisson( const Real  a_time );
  void correctElectricField();
  
  bool advanceE() const { return m_advanceE; }
  bool advanceB() const { return m_advanceB; }

  inline void zeroRHS()
  {
    SpaceUtils::zero( m_magneticField_rhs );
    if (SpaceDim<3) SpaceUtils::zero( m_magneticField_virtual_rhs );
    SpaceUtils::zero( m_electricField_rhs );
    if (SpaceDim<3) SpaceUtils::zero( m_electricField_virtual_rhs );
  }

  inline int getVectorSize(const EMVecType& a_type = e_and_b) const
  {
    int retval(0);
    if (a_type == b_only) {
      retval += m_vec_size_bfield;
      if (SpaceDim < 3) retval += m_vec_size_bfield_virtual;
    } else if (a_type == e_only) {
      retval += m_vec_size_efield;
      if (SpaceDim < 3) retval += m_vec_size_efield_virtual;
    } else {
      retval += m_vec_size_bfield;
      if (SpaceDim < 3) retval += m_vec_size_bfield_virtual;
      retval += m_vec_size_efield;
      if (SpaceDim < 3) retval += m_vec_size_efield_virtual;
    }
    return retval;
  }

  inline int numPicardBlocks() const
  {
    return 2;
  }

  inline int vecOffsetBField() const { return m_vec_offset_bfield; }
  inline int vecOffsetBFieldVirtual() const { return m_vec_offset_bfield_virtual; }
  inline int vecOffsetEField() const { return m_vec_offset_efield; }
  inline int vecOffsetEFieldVirtual() const { return m_vec_offset_efield_virtual; }

  void postTimeStep( const Real  a_dt );
  
  // applyBCs to fields
  void applyBCs_electricField( const Real  a_time );
  void applyBCs_magneticField( const Real  a_time );

  void assemblePrecondMatrix( BandedMatrix&              a_banded_matrix,
                        const bool                       a_use_mass_matrices,
                        const LevelData<EdgeDataBox>&    a_sigma_xx,
#if CH_SPACEDIM==1
                        const LevelData<NodeFArrayBox>&  a_sigma_yy,
#endif
                        const LevelData<NodeFArrayBox>&  a_sigma_zz,
                        const Real                       a_cnormDt );

  inline int numPCMatBands() const 
  {
    return m_pcmat_nbands;
  }
  
private:
  
  // compute the curl of the fields 
  void setCurlB();
  void setCurlE();
  
  const bool m_verbosity; 
  bool m_use_poisson; 
  bool m_enforce_gauss_startup; 
  bool m_external_fields; 
  bool m_advance; 
  bool m_advanceE, m_advanceE_inPlane, m_advanceE_virtual; 
  bool m_advanceB, m_advanceB_inPlane, m_advanceB_virtual; 
  std::array<bool,3> m_advanceE_comp;
  std::array<bool,3> m_advanceB_comp;
  bool m_writeDivs, m_writeCurls, m_writeRho, m_writeExB; 
  Real m_stable_dt;
  Real m_Jnorm_factor; // mu0*qe*cvac*Xscale/Bscale
  Real m_rhoCnorm_factor; // qe/ep0*Xscale/Escale
  Real m_energyE_factor; // Escale*Escale*ep0
  Real m_energyB_factor; // Bscale*Bscale/mu0
  Real m_SdA_factor;   // Escale*Bscale/mu0*area_scale
  Real m_time_scale;   // time_scale
  Real m_volume_scale; 
   
  // surface integrals of Poynting flux used for diagnostics
  RealVect m_intSdA_lo;
  RealVect m_intSdA_hi;
  RealVect m_intSdAdt_lo;
  RealVect m_intSdAdt_hi;

  // grid functions for external fields
  RefCountedPtr<GridFunction> m_gridFunction_extE0;
  RefCountedPtr<GridFunction> m_gridFunction_extE1;
  RefCountedPtr<GridFunction> m_gridFunction_extE2;
  //
  RefCountedPtr<GridFunction> m_gridFunction_extB0;
  RefCountedPtr<GridFunction> m_gridFunction_extB1;
  RefCountedPtr<GridFunction> m_gridFunction_extB2;
  
  const DomainGrid& m_mesh;
  
  FieldBC* m_field_bc;  

  // containers for magnetic field
  LevelData<FluxBox> m_magneticField;          
  LevelData<FArrayBox> m_magneticField_virtual;       // needed for 1D/2D simulations          
  LevelData<FluxBox> m_magneticField_old;          
  LevelData<FArrayBox> m_magneticField_virtual_old;
  
  // containers for electric field
  LevelData<EdgeDataBox> m_electricField;          
  LevelData<NodeFArrayBox> m_electricField_virtual;   // needed for 1D/2D simulations
  LevelData<EdgeDataBox> m_electricField_old;
  LevelData<NodeFArrayBox> m_electricField_virtual_old;
  
  // containers for charge/current density
  LevelData<EdgeDataBox> m_currentDensity;          
  LevelData<NodeFArrayBox> m_currentDensity_virtual;  // needed for 1D/2D simulations
  LevelData<NodeFArrayBox> m_chargeDensity;

  // containers for curlE
  LevelData<FluxBox> m_curlE;          
  LevelData<FArrayBox> m_curlE_virtual;       // needed for 1D/2D simulations
  
  // containers for curlB
  LevelData<EdgeDataBox> m_curlB;          
  LevelData<NodeFArrayBox> m_curlB_virtual;       // needed for 1D/2D simulations          

  // containers for divE and divB
  LevelData<NodeFArrayBox> m_divE;
  LevelData<FArrayBox> m_divB;

  // containers for Poynting_flux: S = ExB
  LevelData<EdgeDataBox> m_ExB;
  LevelData<NodeFArrayBox> m_EvxB;
  
  // poisson solver
  PoissonParameters m_poisson_params;
  Vector<LevelData<NodeFArrayBox>* > m_phi_vector;
  Vector<LevelData<NodeFArrayBox>* > m_rhs_vector;
  AMRMultiGrid<LevelData<NodeFArrayBox> > m_poisson_solver;
  BiCGStabSolver<LevelData<NodeFArrayBox> > m_bottomSolver;
  LevelData<NodeFArrayBox> m_phi0;
  LevelData<NodeFArrayBox> m_rhs0;
  LevelData<EdgeDataBox> m_electricField_correction;

  // RHS arrays for EField
  LevelData<EdgeDataBox> m_electricField_rhs;
  LevelData<NodeFArrayBox> m_electricField_virtual_rhs;

  // RHS arrays for BField
  LevelData<FluxBox> m_magneticField_rhs;
  LevelData<FArrayBox> m_magneticField_virtual_rhs;

  // Vector offsets and sizes
  int m_vec_size_bfield, m_vec_size_bfield_virtual;
  int m_vec_size_efield, m_vec_size_efield_virtual;
  int m_vec_offset_bfield, m_vec_offset_bfield_virtual;
  int m_vec_offset_efield, m_vec_offset_efield_virtual;

  // global DOF
  GlobalDOFEMFields m_gdofs;

  // number of bands in Jacobian matrix
  int m_pcmat_nbands;
  
};


#include "NamespaceFooter.H"

#endif
