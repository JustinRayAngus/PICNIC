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
#include "BinFab.H"
#include "LevelData.H"
#include "FArrayBox.H"

#include "DomainGrid.H"

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

      dataFileIO( ParmParse&         pp,  
                  const DomainGrid&  mesh );
      ~dataFileIO() {;}
      
      void writeParticleDataFile( const ParticleData<JustinsParticle>&  a_Pdata,
                                  const LevelData<FArrayBox>&    a_density, 
                                  const LevelData<FArrayBox>&    a_momentum, 
                                  const LevelData<FArrayBox>&    a_energy, 
                                  const int                      a_cur_step, 
                                  const double                   a_cur_time );
      
      void writeBinFabDataFile( const LevelData<BinFab<JustinsParticle>>&  a_Pdata,
                                const int                      a_cur_step, 
                                const double                   a_cur_time );

   private:

   private:

      bool m_plot_particle_moments;
      const DomainGrid& m_mesh;
      bool m_verbosity;
      
      void writeMeshDataFile();

};

#include "NamespaceFooter.H"

#endif
