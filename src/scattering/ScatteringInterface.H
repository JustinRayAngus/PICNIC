/***
 *
 * ScatteringInterface class header file
 *
 * An instance of this class is owned by System and all
 * scattering operations are handled here
 *
***/

#ifndef ScatteringInterface_h
#define ScatteringInterface_h

#include "PicSpecies.H"
#include "PicSpeciesInterface.H"

#include "Scattering.H"
#include "ScatteringFactory.H"

#include "NamespaceHeader.H"

class ScatteringInterface
{

public:

  ScatteringInterface( const PicSpeciesInterface&  a_pic_species_intf );

  virtual ~ScatteringInterface() {;}

  void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                   const DomainGrid&           a_mesh,
                   const std::string&          a_restart_file_name );

  void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                  const DomainGrid&           a_mesh,
                  const Real                  a_dt_sec );

  void setScatteringProbes();

  void setFusionProbes();

  const ScatteringPtrVect& getScatterPtrVect() { return m_scattering_ptr_vect; }

  int numScatter() const { return m_scattering_ptr_vect.size(); }

  int numCoulomb() const { return m_num_coulomb; }

  int numFusion() const { return m_num_fusion; }

  Real scatterDt( const PicSpeciesInterface&  a_pic_species_intf );

  Real ionizationEnergy() const { return m_izn_energy_joules; }

  Real excitationEnergy() const { return m_exc_energy_joules; }

  Real fusionEnergy() const { return m_fus_energy_joules; }

  void writeCheckpoint( HDF5Handle&  a_handle );

  void readCheckpoint( HDF5Handle&  a_handle,
                 const DomainGrid&  a_mesh );

private:

  void createAllCoulomb( ScatteringPtrVect&  a_coulomb_ptr_vect,
                   const PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                   const ParmParse&          a_pp_scatterC ) const;

  bool m_verbosity;    
  Real m_scatterDt = DBL_MAX;
  WeightMethod m_weight_method;
  int m_model_search_count_max;
  //
  int m_num_coulomb = 0;
  int m_num_elastic = 0;
  int m_num_fusion = 0;
  int m_num_inelastic = 0;
  ScatteringPtrVect m_scattering_ptr_vect;
  //
  Real m_izn_energy_joules = 0.0;
  Real m_exc_energy_joules = 0.0;
  Real m_fus_energy_joules = 0.0;

};

#include "NamespaceFooter.H"

#endif
