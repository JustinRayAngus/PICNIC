#ifndef  _dataFileIO_H_
#define  _dataFileIO_H_

#include <string>
#include <vector>
#include <time.h>

#include "ParmParse.H"
#include "CH_HDF5.H"

#include "Particle.H"
#include "JustinsParticle.H"
#include "ParticleData.H"
#include "PicSpeciesInterface.H"
#include "BinFab.H"
#include "LevelData.H"
#include "FArrayBox.H"

#include "DomainGrid.H"
#include "CodeUnits.H"
#include "ElectroMagneticFields.H"

#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif

#include "NamespaceHeader.H"

//
//  class for reading/writing to hdf5 files
//

class dataFileIO
{
   public:

      dataFileIO( ParmParse&   a_ppsys,  
            const DomainGrid&  a_mesh,
            const CodeUnits&   a_units );
      ~dataFileIO() {;}
      
      void writeElectroMagneticFieldsDataFile( const ElectroMagneticFields&  a_emfield,
#ifdef MASS_MATRIX_TEST
                                               const PicSpeciesInterface&    a_pic_species,
#endif
                                               const int                     a_cur_step, 
                                               const double                  a_cur_time );
      
      void writeEMFields( HDF5Handle&             a_handle,
                    const ElectroMagneticFields&  a_emfield );
      
      void writeEMFields_old( HDF5Handle&             a_handle,
                        const ElectroMagneticFields&  a_emfield );
      
      void writeMassMatrices( HDF5Handle&           a_handle,
                        const PicSpeciesInterface&  a_pic_species );
      
#ifdef MASS_MATRIX_TEST
      void writeMassMatricesTest( HDF5Handle&           a_handle,
                                  HDF5HeaderData&       a_header,
                            const PicSpeciesInterface&  a_pic_species );
#endif 
     
      void writeEMFieldRho( const ElectroMagneticFields&  a_emfield,
                            HDF5Handle&                   a_handle,
                            HDF5HeaderData&               a_header );
      
      void writeEMFieldExB( const ElectroMagneticFields&  a_emfield,
                            HDF5Handle&                   a_handle,
                            HDF5HeaderData&               a_header );
      
      void writeEMFieldDivs( const ElectroMagneticFields&  a_emfield,
                             HDF5Handle&                   a_handle,
                             HDF5HeaderData&               a_header );
      
      void writeEMFieldCurls( const ElectroMagneticFields&  a_emfield,
                              HDF5Handle&                   a_handle,
                              HDF5HeaderData&               a_header );
      
      void writeEMFieldPotential( const ElectroMagneticFields&  a_emfield,
                                  HDF5Handle&                   a_handle,
                                  HDF5HeaderData&               a_header );
      
      void writePicSpecies( const PicSpeciesInterface&  a_pic_species,
                            const int                   a_cur_step,
                            const double                a_cur_time );
      
      void writeNeutralSpeciesDataFile( const PicSpecies&  a_picSpecies,
                                        const int          a_species, 
                                        const int          a_cur_step, 
                                        const double       a_cur_time,
                                        const bool         a_write_energy_off_diag,
                                        const bool         a_write_energy_flux,
                                        const bool         a_write_nppc );
      
      void writeChargedSpeciesDataFile( const PicSpecies&  a_picSpecies,
                                        const int          a_species, 
                                        const int          a_cur_step, 
                                        const double       a_cur_time,
                                        const bool         a_write_charge_density,
                                        const bool         a_write_current_density,
                                        const bool         a_write_energy_off_diag,
                                        const bool         a_write_energy_flux,
                                        const bool         a_write_nppc );
      
      void writeSpeciesParticleFile( const PicSpecies&  a_picSpecies,
                                     const int          a_species, 
                                     const int          a_cur_step, 
                                     const double       a_cur_time );
      
      void writeSpeciesMomentsFile( const PicSpecies&  a_picSpecies,
                                    const int          a_species, 
                                    const int          a_cur_step, 
                                    const double       a_cur_time,
                                    const bool         a_write_charge_density,
                                    const bool         a_write_current_density,
                                    const bool         a_write_energy_off_diag,
                                    const bool         a_write_energy_flux,
                                    const bool         a_write_nppc );
      
      void writeBinFabDataFile( const LevelData<BinFab<JustinsParticle>>&  a_Pdata,
                                const int                      a_cur_step, 
                                const double                   a_cur_time );
      
      void writeCheckpointEMFields( HDF5Handle&             a_handle,
                              const int                     a_write_old_data,
                              const ElectroMagneticFields&  a_emfield );
      
      void writeCheckpointPicSpecies( HDF5Handle&           a_handle,
                                const PicSpeciesInterface&  a_pic_species );
      
      void writeCheckpointParticles( HDF5Handle&  a_handle,
                               const PicSpecies&  a_picSpecies,
                               const int          a_species );
      
   private:

      const DomainGrid& m_mesh;
      const CodeUnits& m_units;
      bool m_verbosity;
      void writeMeshDataFile();

      bool m_nodeDataTest;

};

#include "NamespaceFooter.H"

#endif
