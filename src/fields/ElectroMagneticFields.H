#ifndef ElectroMagneticFields_h
#define ElectroMagneticFields_h

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

#include "PicSpecies.H"

#include "NamespaceHeader.H"

class ElectroMagneticFields
{

public:

  ElectroMagneticFields( ParmParse&   a_ppfields,
                   const DomainGrid&  a_mesh,
                   const CodeUnits&   a_units,
                   const bool&        a_verbosity );
  
  virtual ~ElectroMagneticFields();

  void initialize();
  
  bool advance() const { return m_advance; }
  
  Real stableDt() const { return m_stable_dt; }

  // const references to magnetic field  
  const LevelData<FluxBox>& getMagneticField() const {return m_magneticField;}
  const LevelData<FArrayBox>& getVirtualMagneticField() const {return m_magneticField_virtual;}
  
  // const references to electric field  
  const LevelData<EdgeDataBox>& getElectricField() const {return m_electricField;}
  const LevelData<NodeFArrayBox>& getVirtualElectricField() const {return m_electricField_virtual;}
  
  // const references to current density
  const LevelData<EdgeDataBox>& getCurrentDensity() const {return m_currentDensity;}
  const LevelData<NodeFArrayBox>& getVirtualCurrentDensity() const {return m_currentDensity_virtual;}
  
  // advance the fields in time
  void advanceElectricField( const Real& a_cnormDt );
  void advanceMagneticField( const Real& a_cnormDt );
  
  // convert from half time advance to full advance
  void advanceMagneticField_2ndHalf();
  void advanceElectricField_2ndHalf();
  
  // update the saved old field values
  void updateOldElectricField();
  void updateOldMagneticField();

  // compute the curl of the fields 
  void setCurlB();
  void setCurlE();
 
  // compute the current density for Ampere's law
  void setCurrentDensity( const PicSpeciesPtrVect&  a_pic_species_ptr_vect );
  
private:
  
  const DomainGrid& m_mesh;
  const bool m_verbosity; 
  bool m_advance; 
  Real m_stable_dt;
  Real m_Jnorm_factor; // mu0*qe*cvac*Xscale/Bscale
 
  // containers for magnetic field
  LevelData<FluxBox> m_magneticField;          
  LevelData<FArrayBox> m_magneticField_virtual;       // needed for 1D/2D simulations          
  LevelData<FluxBox> m_magneticField_old;          
  LevelData<FArrayBox> m_magneticField_virtual_old;       // needed for 1D/2D simulations          
  
  // containers for electric field
  LevelData<EdgeDataBox> m_electricField;          
  LevelData<NodeFArrayBox> m_electricField_virtual;   // needed for 1D/2D simulations
  LevelData<EdgeDataBox> m_electricField_old;          
  LevelData<NodeFArrayBox> m_electricField_virtual_old;   // needed for 1D/2D simulations
  
  // containers for current density
  LevelData<EdgeDataBox> m_currentDensity;          
  LevelData<NodeFArrayBox> m_currentDensity_virtual;  // needed for 1D/2D simulations
  
  // containers for curlE
  LevelData<FluxBox> m_curlE;          
  LevelData<FArrayBox> m_curlE_virtual;       // needed for 1D/2D simulations          
  
  // containers for curlB
  LevelData<EdgeDataBox> m_curlB;          
  LevelData<NodeFArrayBox> m_curlB_virtual;       // needed for 1D/2D simulations          
  
};


#include "NamespaceFooter.H"

#endif
