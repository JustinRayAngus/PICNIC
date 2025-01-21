#ifndef  _SYSTEM_H_
#define  _SYSTEM_H_

#include <string>
#include <vector>
#include <time.h>

#include "ParmParse.H"
#include "CH_HDF5.H"

#include "dataFileIO.H"

#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "DomainGrid.H"
#include "CodeUnits.H"

#include "PicSpeciesInterface.H"
#include "ScatteringInterface.H"
#include "EMFields.H"
#include "SpecialOperator.H"

#include "PICTimeIntegratorLibrary.H"

#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif

#include "NamespaceHeader.H"

/**
 * System class.
 *
 * This class owns pretty much everything and orchestrates the 
 * user specified time advancement.
 *
 */

class System 
{
   public:

      System( ParmParse&  a_pp );
      virtual ~System();

      virtual void initialize( const int           a_cur_step, 
                               const Real          a_cur_time,
                               const Real          a_cur_dt,
                               const std::string&  a_restart_file_name );
      
      void writePlotFile( const int   a_cur_step, 
                          const Real  a_cur_time, 
                          const Real  a_cur_dt, 
                          const bool  a_plot_parts );
      
      void writeHistFile( const int   a_cur_step, 
                          const Real  a_cur_time, 
                          const Real  a_cur_dt, 
                          const bool  a_startup_flag);
      
      void writeCheckpointFile( HDF5Handle&  a_handle,
                          const int          a_cur_step, 
                          const Real         a_cur_time, 
                          const Real         a_cur_dt );
      
      void readCheckpointFile( const std::string&  a_chkpt_fname );
      
      bool useParts() {
         if (m_pic_species->numSpecies()>0) { return true; }
         else { return false; }
      }

      bool useScattering() {
         if (m_scattering->numScatter()>0) { return true; }
         else { return false; }
      }
      
      bool useSpecialOps() { return m_use_specialOps; }
      
      bool useFields() { return !m_emfields.isNull(); }

      Real fieldsDt( const int a_step_number );
      
      Real partsDt( const int a_step_number );
      
      Real scatterDt( const int a_step_number );
      
      Real cyclotronDt( const int a_step_number );
      
      Real specialOpsDt( const int a_step_number );
      
      void adaptDt( bool&  a_adapt_dt );

      void preTimeStep( const Real  a_time, 
                        const Real  a_dt,
                        const int   a_step_number );
      
      void timeStep( const Real, const Real, const int );
      
      void postTimeStep( Real&, const Real&, int& );

      void setChargeDensity();
      
      void setSurfaceCharge();

      void scatterParticles( const Real&  a_dt );
     
      inline void printWalltimes() const
      {
        if (!procID()) {
          printf("Time integration walltimes (seconds):\n  %1.3e (pretimestep), %1.3e (timestep), %1.3e (posttimestep), %1.3e (total)\n",
                 m_wt_pretimestep, m_wt_timestep, m_wt_posttimestep,
                 (m_wt_pretimestep+m_wt_timestep+m_wt_posttimestep) );
        }
      }

   private:

      void parseParameters( ParmParse& ppsys );
      
      void createProblemDomain();
      
      void getDisjointBoxLayout( DisjointBoxLayout& grids );
      
      void setupHistFile( const int  a_cur_step );
      
      void createState();
      
      void createEMfields();
      
      void createSpecialOperators();

      PICTIType m_advance_method;
      bool m_implicit_advance;

      DomainGrid* m_mesh;
      CodeUnits*  m_units;
      dataFileIO* m_dataFile;
      PicSpeciesInterface* m_pic_species;
      ScatteringInterface* m_scattering;
     
      RefCountedPtr<EMFields> m_emfields;
      RefCountedPtr<SpecialOperator> m_specialOps;

      bool m_solver_probes;
      bool m_write_wall_time;
      bool m_field_probes;
      bool m_scattering_probes;
      bool m_fusion_probes;
      bool m_bremsstrahlung_probes;
      std::vector<bool> m_species_probes;
      bool m_field_bdry_probes;
      bool m_species_bdry_probes;
      bool m_species_solver_probes;
      int m_history_precision;      

      bool m_verbosity;
      bool m_use_specialOps;

      // problem domain 
      ProblemDomain     m_domain;
      int               m_num_ghosts;
      std::vector<int>  m_num_cells;
      std::vector<bool> m_is_periodic;
      std::vector<int>  m_config_decomp;

      PICTimeIntegrator* m_time_integrator;

      double m_wt_pretimestep, m_wt_timestep, m_wt_posttimestep;
};

#include "NamespaceFooter.H"

#endif
