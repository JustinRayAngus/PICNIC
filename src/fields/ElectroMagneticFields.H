#ifndef ElectroMagneticFields_h
#define ElectroMagneticFields_h

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

#include "NamespaceHeader.H"

class ElectroMagneticFields
{

public:

  ElectroMagneticFields( ParmParse&   a_ppfields,
                   const DomainGrid&  a_mesh,
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
  
  void advanceMagneticField( const Real& a_dt );
  void advanceElectricField( const Real& a_dt );
  
private:
  
  const DomainGrid& m_mesh;
  const bool m_verbosity; 
  bool m_advance; 
  Real m_stable_dt;
 
  // containers for magnetic field
  LevelData<FluxBox> m_magneticField;          
  LevelData<FArrayBox> m_magneticField_virtual;       // needed for 1D/2D simulations          
  
  // containers for electric field
  LevelData<EdgeDataBox> m_electricField;          
  LevelData<NodeFArrayBox> m_electricField_virtual;   // needed for 1D/2D simulations
  
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
