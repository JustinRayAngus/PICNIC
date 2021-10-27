#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <string>
#include <vector>
#include <time.h>

#include "RefCountedPtr.H"
#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "REAL.H"
#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif // CH_USE_TIMER

#include "System.H"

#ifdef USE_ARRAYVIEW
#include "UsingNamespace.H"
extern "C" {
#include <fpu_control.h>
}
/* IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask
  ---(pm is kinda stupid)
*/
static void __attribute__ ((constructor)) trapfpe(void)
{
   fpu_control_t cw =
      _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
   _FPU_SETCW(cw);
}
#endif // USE_ARRAYVIEW

#include "NamespaceHeader.H"

/**
 * Simulation driver class.
 *
 * Responsible for generic initialization, advancement, and finalization of
 * simulation.  This class specifically controls the monitoring of the number
 * of steps and of the simulation time in order to execute plotfile and
 * checkpoint file creation as well as simulation termination.
*/
class Simulation
{
   public:

      Simulation( ParmParse& pp );
      ~Simulation();

      bool notDone();
      void advance();
      void finalize();

   private:

      void initializeTimers();
      void parseParameters( ParmParse& ppsim );
      void printParameters();
      void loadRestartFile( ParmParse& ppsim );
      void writePlotFile();
      inline void writeHistFile(bool startup_flag);
      void writeCheckpointFile();
      void setFixedTimeStep(const Real& a_dt);
      void enforceTimeStep(const Real  a_local_dt);
      void preTimeStep();
      void postTimeStep();
      void writeFiles();
      
      int m_verbosity;

      int m_cur_step;
      int m_max_step;

      Real m_cur_time;
      Real m_max_time;
      Real m_max_wall_time_hrs;
      
      Real m_new_epsilon;

      Real m_cur_dt;
      Real m_fixed_dt;
      Real m_cfl;
      Real m_dt_light;
      Real m_dt_parts;
      Real m_cfl_scatter;
      Real m_dt_scatter;
      int m_dt_scatter_interval;
      bool m_adapt_dt;
      static const Real s_DT_EPS;

      int         m_plot_interval;
      Real        m_plot_time_interval;
      Real        m_plot_time;
      int         m_last_plot;
      std::string m_plot_prefix;

      bool        m_history;
      int         m_history_interval;
      int         m_last_history;
      
      bool        m_checkpoint;
      int         m_checkpoint_interval;
      Real        m_checkpoint_time_interval;
      Real        m_checkpoint_time;
      int         m_last_checkpoint;
      std::string m_checkpoint_prefix;

      std::string m_restart_file_name;

      int m_fixed_random_seed;
      
      int m_command_check_step_interval;
      enum COMMAND {
         NONE,
         STOP,
         ABORT
      };
      COMMAND m_command;
      
      System* m_system;

#ifdef CH_USE_TIMER
      Timer* m_all_timer;
      Timer* m_setup_timer;
      Timer* m_solve_timer;
      Timer* m_shutdown_timer;
#endif

      /* Walltime measurements */
      clock_t m_main_start, m_solve_start, m_main_end, m_solve_end;
};

#include "NamespaceFooter.H"

#endif // _SIMULATION_H_
