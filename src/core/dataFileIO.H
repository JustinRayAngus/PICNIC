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
#include "PicSpecies.H"
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
                                               const int                     a_cur_step, 
                                               const double                  a_cur_time );
      
      void writeEMFieldRho( const ElectroMagneticFields&  a_emfield,
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
      
      void writeNeutralSpeciesDataFile( const PicSpecies&  a_picSpecies,
                                        const int          a_species, 
                                        const int          a_cur_step, 
                                        const double       a_cur_time );
      
      void writeChargedSpeciesDataFile( const PicSpecies&  a_picSpecies,
                                        const int          a_species, 
                                        const int          a_cur_step, 
                                        const double       a_cur_time,
                                        const bool         a_write_charge_density,
                                        const bool         a_write_current_density );
      
      void writeSpeciesParticleFile( const PicSpecies&  a_picSpecies,
                                     const int          a_species, 
                                     const int          a_cur_step, 
                                     const double       a_cur_time );
      
      void writeSpeciesMomentsFile( const PicSpecies&  a_picSpecies,
                                    const int          a_species, 
                                    const int          a_cur_step, 
                                    const double       a_cur_time,
                                    const bool         a_write_charge_density,
                                    const bool         a_write_current_density );
      
      void writeBinFabDataFile( const LevelData<BinFab<JustinsParticle>>&  a_Pdata,
                                const int                      a_cur_step, 
                                const double                   a_cur_time );

   private:

      const DomainGrid& m_mesh;
      const CodeUnits& m_units;
      bool m_verbosity;
      void writeMeshDataFile();

      bool m_nodeDataTest;

};

#include "NamespaceFooter.H"

#endif
