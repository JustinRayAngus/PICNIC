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
#include "CodeUnits.H"

#include "Scattering.H"
#include "ScatteringFactory.H"

#include "NamespaceHeader.H"

class ScatteringInterface
{

public:

  ScatteringInterface( const CodeUnits&          a_units,
                       const PicSpeciesPtrVect&  a_pic_species_ptr_vect );
  
  virtual ~ScatteringInterface() {;}
     
  void initialize( const PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                   const DomainGrid&         a_mesh,
                   const std::string&        a_restart_file_name );
   
  void applyScattering( PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                  const DomainGrid&         a_mesh,
                  const Real                a_dt_sec );

  void setScatteringProbes();

  int numScatter() const { return m_scattering_ptr_vect.size(); }

  Real scatterDt( const PicSpeciesPtrVect&  a_pic_species_ptr_vect );
  
  Real ionizationEnergy() const { return m_izn_energy_joules; }
  
  Real excitationEnergy() const { return m_exc_energy_joules; }

  void writeCheckpoint( HDF5Handle&  a_handle );
  
  void readCheckpoint( HDF5Handle&  a_handle );
  
private:

  void createAllCoulomb( ScatteringPtrVect&  a_coulomb_ptr_vect,
                   const PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                   const ParmParse&          a_pp_scatterC ) const;
  
  bool m_verbosity;    
  Real m_scatterDt = DBL_MAX;
  WeightMethod m_weight_method;
  int m_model_search_count_max;
  //
  int m_num_coulomb;
  int m_num_elastic;
  int m_num_inelastic;
  ScatteringPtrVect m_scattering_ptr_vect;
  //
  Real m_izn_energy_joules = 0.0;
  Real m_exc_energy_joules = 0.0;
 
};

#include "NamespaceFooter.H"

#endif
