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
 *
 * Required input keys and data types:
 *
 *    -\b max_step
 *      positiveinteger value maximum number of simulation steps
 *
 *    -\b max_time
 *      positive real value final time at which simulation terminates
 *
 *    -\b fixed_dt
 *      positive real value, if set, for the time step to be used
 *
 *    -\b max_dt_grow
 *      positive real value of the multiple of the current time step to use for
 *      the next time step.  Must be grater than 1.  Multually exclusive with
 *      fixed_dt
 *
 *    -\b initial_dt_fraction
 *      positive real value fraction of the initial stable time step to use.
 *      Must be less than 1.  Multually exclusive with fixed_dt
 *
 *    -\b checkpoint_interval
 *      integer value specifying the number of steps between checkpoint dumps
 *
 *    -\b checkpoint_prefix
 *      string used as prefix for checkpoint file names ["chk"]
 *
 *    -\b plot_interval
 *      integer value specifying the number of steps between plot dumps
 *
 *    -\b plot_time_interval
 *      Real value specifying the time interval between plot dumps
 *
 *    -\b plot_prefix
 *      string used as prefix for plot file names ["plt"]
 *
 *    -\b verbosity
 *      integer flag specifying the verbosity of logging output; zero turns
 *      of output, and increasing values produce more detailed logs
 *
 * The following represents a sample input entry:
 *
 * \verbatim
 *  simulation.max_step = 1000
 *  simulation.max_time = 1.0
 *  simulation.max_dt_grow = 1.1
 *  simultaion.initial_dt_fraction = 0.1
 *  simulation.checkpoint_interval = 50
 *  simulation.checkpoint_prefix = "chk"
 *  simulation.plot_interval = 10
 *  simulation.plot_prefix = "plt"
 *  simulation.verbosity = 1
 * \endverbatim
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
      void setFixedTimeStep(const Real& a_dt_stable);
      void preTimeStep();
      inline void postTimeStep();
      
      void createProblemDomain( );
      void getDisjointBoxLayout( DisjointBoxLayout& grids );

      int m_verbosity;

      int m_cur_step;
      int m_max_step;

      Real m_cur_time;
      Real m_max_time;

      Real m_cur_dt;
      Real m_fixed_dt;
      Real m_max_dt_grow;
      Real m_cfl;
      bool m_adapt_dt;
      static const Real s_DT_EPS;

      int         m_checkpoint_interval;
      int         m_last_checkpoint;
      std::string m_checkpoint_prefix;

      int         m_plot_interval;
      Real        m_plot_time_interval;
      Real        m_plot_time;
      int         m_last_plot;
      std::string m_plot_prefix;

      bool        m_history;
      int         m_history_interval;
      int         m_last_history;

      ProblemDomain m_domain;
      std::string m_geom_type; 
      int               m_num_ghosts;
      std::vector<int>  m_num_cells;
      std::vector<bool> m_is_periodic;
      std::vector<int>  m_config_decomp;

      //Units m_units;
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
