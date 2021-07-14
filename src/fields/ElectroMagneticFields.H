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

#include "SpaceUtils.H"
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
  bool advanceE() const { return m_advanceE; }
  bool advanceB() const { return m_advanceB; }
  
  bool writeCurls() const { return m_writeCurls; }
  
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
  
  // const references to curl of electric field  
  const LevelData<FluxBox>& getCurlE() const {return m_curlE;}
  const LevelData<FArrayBox>& getCurlE_virtual() const {return m_curlE_virtual;}
  
  // const references to curl of magnetic field  
  const LevelData<EdgeDataBox>& getCurlB() const {return m_curlB;}
  const LevelData<NodeFArrayBox>& getCurlB_virtual() const {return m_curlB_virtual;}
  
  // advance the fields in time
  void advanceElectricField( const Real& a_cnormDt );
  void advanceMagneticField( const Real& a_cnormDt );
  
  // convert from half time advance to full advance
  void advanceMagneticField_2ndHalf( const Real&  a_theta=0.5 );
  void advanceElectricField_2ndHalf( const Real&  a_theta=0.5 );
 
  // save fields for error calc
  void saveElectricField();
  void saveMagneticField();
  
  // applyBCs to fields
  void applyBCs_electricField( const Real  a_time );
  void applyBCs_magneticField( const Real  a_time );
  
  // diff the fields for error calc
  Real diffElectricField();
  Real diffMagneticField();
 
  // update the saved old field values
  void updateOldElectricField();
  void updateOldMagneticField();

  // compute the curl of the fields 
  void setCurlB();
  void setCurlE();
 
  // compute the current density for Ampere's law
  void setCurrentDensity( const PicSpeciesPtrVect&  a_pic_species_ptr_vect );
  
private:
  
  const bool m_verbosity; 
  bool m_advance; 
  bool m_advanceE, m_advanceE_inPlane, m_advanceE_virtual; 
  bool m_advanceB, m_advanceB_inPlane, m_advanceB_virtual; 
  std::array<bool,3> m_advanceE_comp;
  std::array<bool,3> m_advanceB_comp;
  bool m_writeCurls; 
  Real m_stable_dt;
  Real m_Jnorm_factor; // mu0*qe*cvac*Xscale/Bscale
  
  const DomainGrid& m_mesh;
  
  FieldBC* m_electricField_bc;  
  FieldBC* m_magneticField_bc;  
  FieldBC* m_electricField_virtual_bc;  
  FieldBC* m_magneticField_virtual_bc;  

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

  // Diff arrays for EField
  LevelData<EdgeDataBox> m_electricField_diff;
  LevelData<NodeFArrayBox> m_electricField_virtual_diff;

  // Diff arrays for BField
  LevelData<FluxBox> m_magneticField_diff;
  LevelData<FArrayBox> m_magneticField_virtual_diff;

  inline Real computeEFieldDiffNorm() const
  {
    Real retval = SpaceUtils::sumOfSquares(m_electricField_diff);
    if (SpaceDim<3) {
      retval += SpaceUtils::sumOfSquares(m_electricField_virtual_diff);
    }
    return sqrt(retval);
  }

  inline Real computeBFieldDiffNorm() const
  {
    Real retval = SpaceUtils::sumOfSquares(m_magneticField_diff);
    if (SpaceDim<3) {
      retval += SpaceUtils::sumOfSquares(m_magneticField_virtual_diff);
    }
    return sqrt(retval);
  }
  
};


#include "NamespaceFooter.H"

#endif
