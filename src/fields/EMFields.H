#ifndef EMFields_h
#define EMFields_h

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
#include "ODEVector.H"

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

class EMFields
{

public:

  EMFields( ParmParse&,
                   const DomainGrid&,
                   const CodeUnits&,
                   const bool&,
                   const EMVecType& a = e_and_b);
  
  virtual ~EMFields();

  void initialize( const Real          a_cur_time,
                   const std::string&  a_restart_file_name );
  
  void initializeMassMatricesForPC( const IntVect&  a_ncomp_xx,
		                    const IntVect&  a_ncomp_xy,
		                    const IntVect&  a_ncomp_xz,
		                    const IntVect&  a_ncomp_yx,
		                    const IntVect&  a_ncomp_yy,
		                    const IntVect&  a_ncomp_yz,
		                    const IntVect&  a_ncomp_zx,
		                    const IntVect&  a_ncomp_zy,
		                    const IntVect&  a_ncomp_zz );

  void initExternalFields();
  
  bool advance() const { return m_advance; }
  
  bool externalFields() const { return m_external_fields; }
  
  std::array<bool,3> advanceE_comp() const { return m_advanceE_comp; }
  std::array<bool,3> advanceB_comp() const { return m_advanceB_comp; }
  
  bool useFiltering() const { return m_use_filtering; }
  bool usePoisson() const { return m_use_poisson; }
  
  bool writeDivs() const { return m_writeDivs; }
  bool writeCurls() const { return m_writeCurls; }
  bool writeRho() const { return m_writeRho; }
  bool writeSigma() const { return m_writeSigma; }
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
  const LevelData<NodeFArrayBox>& getSurfaceCharge() const {return m_surfaceCharge;}
  
  // non-const references to charge/current density
  LevelData<EdgeDataBox>& getCurrentDensity() {return m_currentDensity;}
  LevelData<NodeFArrayBox>& getVirtualCurrentDensity() {return m_currentDensity_virtual;}
  LevelData<NodeFArrayBox>& getChargeDensity() {return m_chargeDensity;}
  LevelData<NodeFArrayBox>& getSurfaceCharge() {return m_surfaceCharge;}
  
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
  
  // non-const references to sigmas for PC 
  LevelData<EdgeDataBox>& getSigmaxxPC() {return m_sigma_xx_pc;}
  LevelData<EdgeDataBox>& getSigmaxyPC() {return m_sigma_xy_pc;}
  LevelData<EdgeDataBox>& getSigmaxzPC() {return m_sigma_xz_pc;}
#if CH_SPACEDIM==1
  LevelData<NodeFArrayBox>& getSigmayxPC() {return m_sigma_yx_pc;}
  LevelData<NodeFArrayBox>& getSigmayyPC() {return m_sigma_yy_pc;}
  LevelData<NodeFArrayBox>& getSigmayzPC() {return m_sigma_yz_pc;}
#endif
  LevelData<NodeFArrayBox>& getSigmazxPC() {return m_sigma_zx_pc;}
  LevelData<NodeFArrayBox>& getSigmazyPC() {return m_sigma_zy_pc;}
  LevelData<NodeFArrayBox>& getSigmazzPC() {return m_sigma_zz_pc;}
  
  const int getMassMatrixPCwidth() const {return m_pc_mass_matrix_width;}
  const bool includeMassMatrixij() const {return m_pc_mass_matrix_include_ij;}
  
  const IntVect& getNcompxxPC() const {return m_ncomp_xx_pc;}
  const IntVect& getNcompxyPC() const {return m_ncomp_xy_pc;}
  const IntVect& getNcompxzPC() const {return m_ncomp_xz_pc;}
  const IntVect& getNcompyxPC() const {return m_ncomp_yx_pc;}
  const IntVect& getNcompyyPC() const {return m_ncomp_yy_pc;}
  const IntVect& getNcompyzPC() const {return m_ncomp_yz_pc;}
  const IntVect& getNcompzxPC() const {return m_ncomp_zx_pc;}
  const IntVect& getNcompzyPC() const {return m_ncomp_zy_pc;}
  const IntVect& getNcompzzPC() const {return m_ncomp_zz_pc;}

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
 
  // get dir derivative of external electric field
  std::array<Real,3> getExternaldEdX( const RealVect&  a_X,
                                      const int        a_dir ) const
  {  
     std::array<Real,3> extdEdX = {0.0,0.0,0.0};
     if(m_external_fields) {
        extdEdX[0] = m_gridFunction_extE0->getDerivative(a_X,a_dir);
        extdEdX[1] = m_gridFunction_extE1->getDerivative(a_X,a_dir);
        extdEdX[2] = m_gridFunction_extE2->getDerivative(a_X,a_dir);
     }
     return extdEdX;
  }
  
  // get dir derivative of external magnetic field
  std::array<Real,3> getExternaldBdX( const RealVect&  a_X,
                                      const int        a_dir ) const
  {  
     std::array<Real,3> extdBdX = {0.0,0.0,0.0};
     if(m_external_fields) {
        extdBdX[0] = m_gridFunction_extB0->getDerivative(a_X,a_dir);
        extdBdX[1] = m_gridFunction_extB1->getDerivative(a_X,a_dir);
        extdBdX[2] = m_gridFunction_extB2->getDerivative(a_X,a_dir);
     }
     return extdBdX;
  }

  // advance the fields in time
  void computeRHSElectricField( const Real, const Real );
  void computeRHSMagneticField( const Real, const Real );
  
  // compute the poynting flux 
  void setPoyntingFlux();
  
  // set the filtered fields
  void setFilteredFields() const {
    const DisjointBoxLayout& grids(m_mesh.getDBL());
    for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         m_electricField_filtered[dit][dir].copy(m_electricField[dit][dir]);
      }
    }
    SpaceUtils::applyBinomialFilter(m_electricField_filtered);
    const Real dummy_time = 0.0;
    SpaceUtils::exchangeEdgeDataBox(m_electricField_filtered);
    m_field_bc->applyEdgeBC( m_electricField_filtered, dummy_time );
  }
  const LevelData<EdgeDataBox>& getFilteredElectricField() const {
    if(m_use_filtering) {return m_electricField_filtered;}
    else {return m_electricField;}
  }
  
  void setChargeDensity( const LevelData<NodeFArrayBox>&  a_pic_rho ) {
    const DisjointBoxLayout& grids(m_mesh.getDBL());
    for(DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& rho = m_chargeDensity[dit].getFab();
      const FArrayBox& this_pic_rho = a_pic_rho[dit].getFab();
      rho.copy(this_pic_rho,0,0,rho.nComp());
    }
  }
  
  void setSurfaceCharge( const LevelData<NodeFArrayBox>&  a_pic_sigma ) {
    const DisjointBoxLayout& grids(m_mesh.getDBL());
    for(DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& sigma = m_surfaceCharge[dit].getFab();
      const FArrayBox& this_pic_sigma = a_pic_sigma[dit].getFab();
      sigma.copy(this_pic_sigma,0,0,sigma.nComp());
    }
  }
  
  void setCurrentDensity( const LevelData<EdgeDataBox>&    a_pic_J,
		          const LevelData<NodeFArrayBox>&  a_pic_Jv ) {
    const DisjointBoxLayout& grids(m_mesh.getDBL());
    for(DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
        FArrayBox& this_J = m_currentDensity[dit][dir];
        const FArrayBox& this_pic_J = a_pic_J[dit][dir];
        this_J.copy(this_pic_J);
      }
      FArrayBox& this_Jv = m_currentDensity_virtual[dit].getFab();
      const FArrayBox& this_pic_Jv = a_pic_Jv[dit].getFab();
      this_Jv.copy(this_pic_Jv,0,0,this_Jv.nComp());
    }
  }
  
  
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
  
  // compute the curl of the fields 
  void setCurlB();
  void setCurlE();

  void enforceGaussLaw( const LevelData<NodeFArrayBox>&  a_rho,
		        const Real                       a_time );
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
      
  void copyEFromVec( const ODEVector<EMFields>&  a_vec )
  {
    int offset(0);
    copyEFromVec( a_vec, offset );
  }
  void copyBFromVec( const ODEVector<EMFields>& a_vec)
  {
    int offset(0);
    copyBFromVec( a_vec, offset );
  }
  void copyBToVec( ODEVector<EMFields>& a_vec ) const
  {
    int offset(0);
    copyBToVec( a_vec, offset );
  }
  void copyEToVec( ODEVector<EMFields>& a_vec ) const
  {
    int offset(0);
    copyEToVec( a_vec, offset );
  }

  void copyEoldFromVec( const ODEVector<EMFields>&  a_vec )
  {
    int offset(0);
    copyEoldFromVec( a_vec, offset );
  }
  void copyBoldFromVec( const ODEVector<EMFields>& a_vec)
  {
    int offset(0);
    copyBoldFromVec( a_vec, offset );
  }
  void copyBoldToVec( ODEVector<EMFields>& a_vec ) const
  {
    int offset(0);
    copyBoldToVec( a_vec, offset );
  }
  void copyEoldToVec( ODEVector<EMFields>& a_vec ) const
  {
    int offset(0);
    copyEoldToVec( a_vec, offset );
  }
  
  void copyERHSFromVec( const ODEVector<EMFields>&  a_vec )
  {
    int offset(0);
    copyERHSFromVec( a_vec, offset );
  }
  void copyBRHSFromVec( const ODEVector<EMFields>& a_vec)
  {
    int offset(0);
    copyBRHSFromVec( a_vec, offset );
  }
  void copyBRHSToVec( ODEVector<EMFields>& a_vec ) const
  {
    int offset(0);
    copyBRHSToVec( a_vec, offset );
  }
  void copyERHSToVec( ODEVector<EMFields>& a_vec ) const
  {
    int offset(0);
    copyERHSToVec( a_vec, offset );
  }

  void copyEFromVec( const ODEVector<EMFields>&  a_vec,
                           int&                  a_offset )
  {
    a_offset += SpaceUtils::copyToLevelData( m_electricField, a_vec.dataAt(a_offset) );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyToLevelData( m_electricField_virtual, a_vec.dataAt(a_offset) );
#endif
    return;
  }
  void copyBFromVec( const ODEVector<EMFields>&  a_vec,
                           int&                  a_offset )
  {
    a_offset += SpaceUtils::copyToLevelData( m_magneticField, a_vec.dataAt(a_offset) );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyToLevelData( m_magneticField_virtual, a_vec.dataAt(a_offset) );
#endif
    return;
  }
  void copyEToVec( ODEVector<EMFields>&  a_vec,
                   int&                  a_offset ) const
  {
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_electricField );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_electricField_virtual );
#endif
    return;
  }
  void copyBToVec( ODEVector<EMFields>&  a_vec,
                   int&                  a_offset ) const
  {
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_magneticField );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_magneticField_virtual );
#endif
    return;
  }

  void copyEoldFromVec( const ODEVector<EMFields>&  a_vec,
                              int&                  a_offset )
  {
    a_offset += SpaceUtils::copyToLevelData( m_electricField_old, a_vec.dataAt(a_offset) );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyToLevelData( m_electricField_virtual_old, a_vec.dataAt(a_offset) );
#endif
    return;
  }
  void copyBoldFromVec( const ODEVector<EMFields>&  a_vec,
                              int&                  a_offset )
  {
    a_offset += SpaceUtils::copyToLevelData( m_magneticField_old, a_vec.dataAt(a_offset) );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyToLevelData( m_magneticField_virtual_old, a_vec.dataAt(a_offset) );
#endif
    return;
  }
  void copyEoldToVec( ODEVector<EMFields>&  a_vec,
                      int&                  a_offset ) const
  {
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_electricField_old );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_electricField_virtual_old );
#endif
    return;
  }
  void copyBoldToVec( ODEVector<EMFields>&  a_vec,
                      int&                  a_offset ) const
  {
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_magneticField_old );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_magneticField_virtual_old );
#endif
    return;
  }

  void copyERHSFromVec( const ODEVector<EMFields>&  a_vec,
                              int&                  a_offset )
  {
    a_offset += SpaceUtils::copyToLevelData( m_electricField_rhs, a_vec.dataAt(a_offset) );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyToLevelData( m_electricField_virtual_rhs, a_vec.dataAt(a_offset) );
#endif
    return;
  }
  void copyBRHSFromVec( const ODEVector<EMFields>&  a_vec,
                              int&                  a_offset )
  {
    a_offset += SpaceUtils::copyToLevelData( m_magneticField_rhs, a_vec.dataAt(a_offset) );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyToLevelData( m_magneticField_virtual_rhs, a_vec.dataAt(a_offset) );
#endif
    return;
  }
  void copyERHSToVec( ODEVector<EMFields>&  a_vec,
                      int&                  a_offset ) const
  {
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_electricField_rhs );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_electricField_virtual_rhs );
#endif
    return;
  }
  void copyBRHSToVec( ODEVector<EMFields>&  a_vec,
                      int&                  a_offset ) const
  {
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_magneticField_rhs );
#if CH_SPACEDIM<3
    a_offset += SpaceUtils::copyFromLevelData( a_vec.dataAt(a_offset), m_magneticField_virtual_rhs );
#endif
    return;
  }

  void updatePhysicalState( ODEVector<EMFields>&  a_U,
                      const Real                  a_time,
                      const EMVecType&            a_vec_type )
  {
    CH_TIME("EMFields::updatePhysicalState()");
    if (a_vec_type == b_only) {
      copyBFromVec( a_U );
      applyBCs_magneticField( a_time );
      copyBToVec( a_U );
    } else if (a_vec_type == e_only ) {
      copyEFromVec( a_U );
      applyBCs_electricField( a_time );
      copyEToVec( a_U );
    } else {
      int offset(0);
      copyBFromVec( a_U, offset );
      copyEFromVec( a_U, offset );
      applyBCs_magneticField( a_time );
      applyBCs_electricField( a_time );
      offset = 0;
      copyBToVec( a_U, offset );
      copyEToVec( a_U, offset );
    }
  }

  inline int vecOffsetBField() const { return m_vec_offset_bfield; }
  inline int vecOffsetBFieldVirtual() const { return m_vec_offset_bfield_virtual; }
  inline int vecOffsetEField() const { return m_vec_offset_efield; }
  inline int vecOffsetEFieldVirtual() const { return m_vec_offset_efield_virtual; }

  void updateBoundaryProbes( const Real  a_dt );
  
  // applyBCs to fields
  void applyBCs_electricField( const Real  a_time );
  void applyBCs_magneticField( const Real  a_time );

  void assemblePrecondMatrix( BandedMatrix&,
                              const bool,
                              const Real,
                              const EMVecType& );

  void computePrecondMatrixNNZ( const EMVecType& );

  inline int numPCMatBands() const 
  {
    return m_pcmat_nnz;
  }
  
private:
  
  const bool m_verbosity; 
  bool m_use_filtering; 
  bool m_use_poisson; 
  bool m_enforce_gauss_startup; 
  bool m_external_fields; 
  bool m_advance; 
  bool m_advanceE, m_advanceE_inPlane, m_advanceE_virtual; 
  bool m_advanceB, m_advanceB_inPlane, m_advanceB_virtual; 
  std::array<bool,3> m_advanceE_comp;
  std::array<bool,3> m_advanceB_comp;
  bool m_writeDivs, m_writeCurls, m_writeRho, m_writeSigma, m_writeExB; 
  Real m_stable_dt;
  Real m_cvacNorm;
  Real m_Jnorm_factor; // mu0*qe*cvac*Xscale/Bscale
  Real m_rhoCnorm_factor; // qe/ep0*Xscale/Escale
  Real m_energyE_factor; // Escale*Escale*ep0
  Real m_energyB_factor; // Bscale*Bscale/mu0
  Real m_SdA_factor;   // Escale*Bscale/mu0*area_scale
  Real m_time_scale;   // time_scale
   
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
  
  // containers for BC masking in preconditioner matrix
  LevelData<FluxBox> m_PC_mask_B;
  LevelData<FArrayBox> m_PC_mask_Bv;
  LevelData<EdgeDataBox> m_PC_mask_E;
  LevelData<NodeFArrayBox> m_PC_mask_Ev;

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
  
  // containers for filtered fields
  mutable LevelData<EdgeDataBox> m_electricField_filtered;
  
  // containers for charge/current density
  LevelData<EdgeDataBox> m_currentDensity;          
  LevelData<NodeFArrayBox> m_currentDensity_virtual;  // needed for 1D/2D simulations
  LevelData<NodeFArrayBox> m_chargeDensity;
  LevelData<NodeFArrayBox> m_surfaceCharge;

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

  // number of nonzeros per row in Jacobian matrix
  int m_pcmat_nnz;
  // number of nearest neighbors to include in
  // mass matrix contributions to preconditioner
  int m_pc_mass_matrix_width;
  // include sigma_ij terms in mass matrix contribution
  // to preconditioner
  bool m_pc_mass_matrix_include_ij;

  LevelData<EdgeDataBox> m_sigma_xx_pc;
  LevelData<EdgeDataBox> m_sigma_xy_pc;
  LevelData<EdgeDataBox> m_sigma_xz_pc;
#if CH_SPACEDIM==1
  LevelData<NodeFArrayBox> m_sigma_yx_pc;
  LevelData<NodeFArrayBox> m_sigma_yy_pc;
  LevelData<NodeFArrayBox> m_sigma_yz_pc;
#endif
  LevelData<NodeFArrayBox> m_sigma_zx_pc;
  LevelData<NodeFArrayBox> m_sigma_zy_pc;
  LevelData<NodeFArrayBox> m_sigma_zz_pc;
  IntVect m_ncomp_xx_pc, m_ncomp_xy_pc, m_ncomp_xz_pc;
  IntVect m_ncomp_yx_pc, m_ncomp_yy_pc, m_ncomp_yz_pc;
  IntVect m_ncomp_zx_pc, m_ncomp_zy_pc, m_ncomp_zz_pc;

  void assemblePrecondMatrixB(  BandedMatrix&,
                                const bool,
                                const Real );

  void assemblePrecondMatrixE(  BandedMatrix&,
                                const bool,
                                const bool,
                                const Real );

};


#include "NamespaceFooter.H"

#endif
