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
#include "LevelData.H"
#include "DomainGrid.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"

#include "ParticleData.H"
#include "Particle.H"
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "MeshInterp.H"

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

  void initialize();

  virtual Real mass() const { return m_mass; }
  virtual Real Uint() const { return m_Uint; }
  virtual Real charge() const { return m_charge; }
  virtual Real stableDt() const { return m_stable_dt; }

  virtual const std::string& name() const { return m_name; }

  virtual bool motion() const { return m_motion; }
  virtual bool forces() const { return m_forces; }
  
  virtual bool isSpecies( const string& name ) const;

  virtual ParticleData<Particle>& partData() {return m_data;}
  
  virtual const ParticleData<Particle>& partData() const {return m_data;}
  
  // ways to get/set number density member
  //
  virtual void numberDensity( LevelData<FArrayBox>& a_rho );
  
  virtual const LevelData<FArrayBox>& getNumberDensity( const bool setData )
  { 
    if(setData) setNumberDensity();
    return m_density;
  }
  
  void setNumberDensity();
 
  void advancePositions( const Real&  a_dt );

  //const PicSpecies& operator=( const PicSpecies& a_species );
  
private:
  
  Real m_mass;     // species mass [9.1094e-31 kg] 
  Real m_Uint;     // species internal energy [eV]
  Real m_charge;   // species charge [1.6022e-13 muC]
 
  Real rpiston = 1.0;
  Real vpiston = -1.0; 
  Real m_stable_dt = DBL_MAX; 

  std::string m_name;  // species name

  bool m_motion;
  bool m_forces;

  const DomainGrid& m_mesh;
  const MeshInterp& m_meshInterp;
  
  //MeshInterp* m_meshInterp;
  
  ParticleData<Particle> m_data; // container for storing list of particle data

  mutable LevelData<FArrayBox> m_density;     // density moment at cell center
  LevelData<FArrayBox> m_velocity;    // mean velocity at cell center
  LevelData<FArrayBox> m_temperature; // temperature at cell center


  void createMeshInterp();

};

#include "NamespaceFooter.H"

#endif
