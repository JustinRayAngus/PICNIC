#include "Simulation.H"
#include <limits>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <sys/stat.h>

#include "MathUtils.H"
#include "LoadBalance.H"

#include "NamespaceHeader.H"

Simulation::Simulation( ParmParse& a_pp, const std::time_t a_main_start )
   :   m_verbosity(0),
       m_cur_step(0),
       m_max_step(0),
       m_cur_time(0.0),
       m_max_time(0.0),
       m_max_wall_time_hrs(24.0),
       m_new_epsilon(s_DT_EPS),
       m_cur_dt(DBL_MAX),
       m_fixed_dt(-1.0),
       m_dt_light(DBL_MAX),
       m_dt_parts(DBL_MAX),
       m_dt_scatter(DBL_MAX),
       m_dt_specialOps(DBL_MAX),
       m_cfl(1.0),
       m_cfl_scatter(0.1),
       m_dt_scatter_interval(10),
       m_dt_check_interval(10),
       m_adapt_dt(true),
       m_plot_interval(0),
       m_plot_time_interval(0.0),
       m_plot_time(0.0),
       m_last_plot(0),
       m_plot_prefix( "plt" ),
       m_history(false),
       m_history_interval(1),
       m_last_history(0),
       m_checkpoint(false),
       m_checkpoint_interval(0),
       m_checkpoint_time_interval(0.0),
       m_checkpoint_time(0.0),
       m_last_checkpoint(0),
       m_checkpoint_prefix( "chk" ),
       m_fixed_random_seed(-1),
       m_command_check_step_interval(100),
       m_command(NONE),
       m_system( NULL )
#ifdef CH_USE_TIMER
,
       m_all_timer( NULL ),
       m_setup_timer( NULL ),
       m_solve_timer( NULL ),
       m_shutdown_timer( NULL )
#endif
{
   m_main_start = a_main_start;
#ifdef CH_USE_TIMER
   initializeTimers();
   m_all_timer->start();
   m_setup_timer->start();
#endif

#ifdef USE_ARRAYVIEW
   trapfpe();
#endif

   // parse simulation specific parameters from input file
   ParmParse ppsim( "simulation" );
   parseParameters( ppsim );

   // check if doing restart
   if (ppsim.contains( "restart_file" ) ) {
      loadRestartFile( ppsim );
   }
   
   if (m_verbosity) printParameters();

   // initialize the entire system 
   // (domain, mesh, species, fields, operators, ibcs...)
   m_system = new System( a_pp );
   m_system->initialize(m_cur_step, m_cur_time, m_restart_file_name);
   m_restart_step = m_cur_step;

   // write t=0 (or at restart time) values to plot files
   if(m_plot_interval>=0 || m_plot_time_interval>=0.0) {
      writePlotFile();
      m_last_plot = m_cur_step;
      if ( m_plot_time_interval>=0.0 ) { 
         m_plot_time = m_cur_time + m_plot_time_interval;
         if(!m_restart_file_name.empty()) {
            double delta_time;
            delta_time = fmod(m_plot_time, m_plot_time_interval);
            m_plot_time = m_plot_time - delta_time;
         }
      }
   }
   if(m_history) {
      writeHistFile(true);
      m_last_history = m_cur_step;
   }
   if(m_checkpoint_time_interval>=0.0) {
      m_checkpoint_time = m_cur_time + m_checkpoint_time_interval;
      if(!m_restart_file_name.empty()) {
         double delta_time;
         delta_time = fmod(m_checkpoint_time, m_checkpoint_time_interval);
         m_checkpoint_time = m_checkpoint_time - delta_time;
      }
   }
      
   //bool useSpecialOps() { return m_use_specialOps; }
 
   // initialize some time step
   m_dt_light = m_system->fieldsDt( m_cur_step );
   if(m_fixed_dt>0.0) {
      if(m_cur_dt>m_dt_light*m_cfl && !procID()) {
         std::cout << "WARNING: m_cur_dt > m_dt_light*m_cfl" << std::endl;
      }
      //CH_assert(m_cur_dt<=m_dt_light*m_cfl);
   }
   else {
      m_cur_dt = std::min( m_cur_dt, m_dt_light*m_cfl );
      m_system->adaptDt(m_adapt_dt);
   }
   if(m_adapt_dt) m_dt_check_interval = 1;

   enforceTimeStep( m_cur_dt );

#ifdef CH_USE_TIMER
   setup_timer->stop();
#endif

#ifdef CH_USE_TIMER
   solve_timer->start();
#endif
   auto now = std::chrono::system_clock::now();
   m_solve_start = std::chrono::system_clock::to_time_t(now);

}

Simulation::~Simulation()
{
   delete m_system;
#ifdef CH_USE_TIMER
   delete m_all_timer;
   delete m_setup_timer;
   delete m_solve_timer;
   delete m_shutdown_timer;
#endif
}

const Real
Simulation::s_DT_EPS = std::numeric_limits<Real>::epsilon();


void Simulation::printParameters()
{
   if(!procID()) {
      std::cout << "maximum step = " << m_max_step << endl;
      std::cout << "maximum time = " << m_max_time << endl;
      std::cout << "maximum wall time (hrs) = " << m_max_wall_time_hrs << endl;
      std::cout << "s_DT_EPS = " << s_DT_EPS << endl;
      if(!m_adapt_dt) std::cout << "fixed dt = " << m_fixed_dt << endl;
      if(m_plot_time_interval>0.0) {
         std::cout << "plot time interval = " << m_plot_time_interval << endl;
      }
      else {
         std::cout << "plot step interval = " << m_plot_interval << endl;
      }
      if(m_checkpoint) {
         if(m_checkpoint_time_interval>0.0) {
            std::cout << "checkpoint time interval = " << m_checkpoint_time_interval << endl;
         }
         else {
            std::cout << "checkpoint interval = " << m_checkpoint_interval << endl;
         }
      }
      if(m_cfl_scatter) std::cout << "scatter dt multiplier = " << m_cfl_scatter << endl;
      if(m_cfl_scatter) std::cout << "scatter dt check interval = " << m_dt_scatter_interval << endl;
      if(m_history)     std::cout << "history step interval = " << m_history_interval << endl;
      if(m_fixed_random_seed>=0)     std::cout << "fixed random seed = " << m_fixed_random_seed<< endl;
      std::cout << endl;
   }
}


void Simulation::loadRestartFile( ParmParse& a_ppsim )
{
   a_ppsim.query( "restart_file", m_restart_file_name );

   if(!procID()) cout << "Reading restart file " << m_restart_file_name << endl;
   
   ifstream f(m_restart_file_name.c_str());
   if(!f.good()) {
      if(!procID()) { 
         cout << "EXIT FAILURE!!!" << endl;
         cout << "restart file = " << m_restart_file_name;
         cout << " not found" << endl << endl;
      }
      exit(EXIT_FAILURE);
   }

#ifdef CH_USE_HDF5
   HDF5Handle handle( m_restart_file_name, HDF5Handle::OPEN_RDONLY );

   HDF5HeaderData header;
   header.readFromFile( handle );

   m_cur_step = header.m_int["cur_step"];
   m_cur_time = header.m_real["cur_time"];
   //m_cur_dt   = header.m_real["cur_dt"];

   handle.close();

#else
   MayDay::Error("restart only defined with hdf5");
#endif

   if(!procID()) {
      cout << "restart at step = " << m_cur_step << endl;
      cout << "restart at time = " << m_cur_time << endl;
      cout << "restart dt = " << m_cur_dt << endl << endl;
   }

}

//static int plottedonce = 0;

void Simulation::writePlotFile()
{
   if (!procID() && m_verbosity) {
      cout << "writing plot file at step " << m_cur_step << endl << endl;
   }

#ifdef CH_USE_HDF5
   m_system->writePlotFile(m_cur_step, m_cur_time, m_cur_dt);
#else
   MayDay::Error( "plotting only defined with hdf5" );
#endif

}

inline void Simulation::writeHistFile(bool startup_flag)
{
   // the startup_flag is used to force a write on start-up
   m_system->writeHistFile(m_cur_step, m_cur_time, m_cur_dt, startup_flag);

}

void Simulation::writeCheckpointFile()
{

   // need to write history at checkpoint time in order
   // for certain cummulative probes to be correct on restart
   if(m_history && m_last_history!=m_cur_step) {
      writeHistFile(false);
      m_last_history = m_cur_step;
   }

   std::string chk_dir("checkpoint_data");
#ifdef CH_MPI
   if (procID() == 0) {
#endif
      // only works the first time, subsequent failure is normal and expected
      mkdir( chk_dir.c_str(), 0777 );
#ifdef CH_MPI
   }
#endif
   std::string dir_prefix = std::string( chk_dir + "/" );

#ifdef CH_USE_HDF5
   char buffer[100];
   sprintf( buffer, "%s%06d.%dd.hdf5",
            m_checkpoint_prefix.c_str(), m_cur_step, SpaceDim );
            
   std::string filename(dir_prefix + buffer);

   HDF5Handle handle( filename, HDF5Handle::CREATE);
   m_system->writeCheckpointFile( handle, m_cur_step, m_cur_time, m_cur_dt );
   handle.close();
#else
   MayDay::Error( "file writing only defined with hdf5" );
#endif

}

void Simulation::initializeTimers()
{
#ifdef CH_USE_TIMER
   m_all_timer      = new Timer("All",0);
   m_setup_timer    = new Timer("Setup",*m_all_timer);
   m_solve_timer    = new Timer("Solve",*m_all_timer,1);
   m_shutdown_timer = new Timer("Shutdown",*n_all_timer);
   Timer::TimerInit( 0 );
#endif
}

bool Simulation::notDone()
{
   CH_TIMERS("Simulation::notDone()");
   
   //  look for command file with "stop" or "abort" command
   if((m_cur_step % m_command_check_step_interval)==0 ) {
      const std::string fileName("command");
      struct stat buffer;
      if(procID()==0 && (stat (fileName.c_str(),&buffer)==0)) {
         ifstream commandFile;
         commandFile.open(fileName,ios::in);
         std::string line;
         while(getline(commandFile,line)) {
            if ( line.compare("stop") == 0 ) {
               cout << " terminating simulation: found stop command " << endl;
               m_command = STOP;
               break;
            }
            if ( line.compare("abort") == 0 ) {
               m_command = ABORT;
               cout << " terminating simulation: found abort command " << endl;
               break;
            }
         }
         commandFile.close();
         std::remove(fileName.c_str());
      }
      
      // check run time vs max wall time
      auto now = std::chrono::system_clock::now();
      std::time_t now_time = std::chrono::system_clock::to_time_t(now);
      double wall_time_sec = ((double) (now_time - m_main_start));
      double wall_time_hrs = wall_time_sec/3600.0;
      double wall_time_hrs_g;

#ifdef CH_MPI
      MPI_Allreduce(&wall_time_hrs,&wall_time_hrs_g,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
      wall_time_hrs_g = wall_time_hrs;
#endif
      if(wall_time_hrs_g > m_max_wall_time_hrs && !procID()){
         cout << endl;
         cout << "Maximum wall time reached. Terminating simulation..." << endl;
         cout << endl;
         m_command = STOP;
      }

      // broadcast command to all other processors
#ifdef CH_MPI
      MPI_Bcast( &m_command, 1, MPI_INT, 0, MPI_COMM_WORLD ); 
#endif
   }

   if(m_command==STOP || m_command==ABORT) {
      if(m_command==ABORT) {
         m_last_plot = m_cur_step;
         m_last_checkpoint = m_cur_step;
      }
      return 0;
   }
   else {
      //Real delta = m_cur_time - m_max_time;
      //Real delta = m_cur_time - m_max_time*(1.0-s_DT_EPS);
      Real delta = m_cur_time - m_max_time + m_new_epsilon;
      return ( (m_cur_step<m_max_step) && (delta<0.0) );
   }
   
}

void Simulation::advance()
{
   CH_TIMERS("Simulation::advance()");
   CH_TIMER("print_diagnostics",t_print_diagnostics);
   int cout_step_interval = 1;
   if(m_cur_step<10) cout_step_interval = 1;
   else if(m_cur_step<100) cout_step_interval = 10;
   else cout_step_interval = 100;
   
   if (!procID() && (m_cur_step==0 || m_cur_step==m_restart_step) ) {
      cout << endl << "====================================================================" << endl;
   }

   preTimeStep();
   if(!procID() && ((m_cur_step+1) % cout_step_interval)==0 ) {
      cout << endl << "Step " << m_cur_step+1 << ", dt = " << m_cur_dt << endl;
   }
   m_system->timeStep( m_cur_time, m_cur_dt, m_cur_step );
   postTimeStep();
   
   auto now = std::chrono::system_clock::now();
   std::time_t now_time = std::chrono::system_clock::to_time_t(now);

   double walltime, walltime_g;
   walltime = ((double) (now_time - m_solve_start));
#ifdef CH_MPI
   MPI_Allreduce(&walltime,&walltime_g,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
   walltime_g = walltime;
#endif

   CH_START(t_print_diagnostics);
   //m_system->printDiagnostics();
   if(!procID() && (m_cur_step % cout_step_interval)==0 ) {
      cout << "Step " << m_cur_step 
           << " completed, simulation time is " << m_cur_time 
           << ", solver wall time is " << walltime_g << " seconds" << endl;
      if(m_system->useScattering()) cout << "mean free scattering time  = " << m_dt_scatter << endl;
      if(m_system->useParts()) cout << "particle courant time step = " << m_dt_parts << endl;
      if(m_system->useSpecialOps()) cout << "special ops time step = " << m_dt_specialOps << endl;
      cout << "----" << endl;
   }
   CH_STOP(t_print_diagnostics);
   
   writeFiles();
}

void Simulation::finalize()
{
   auto now = std::chrono::system_clock::now();
   std::time_t solve_end = std::chrono::system_clock::to_time_t(now);
#ifdef CH_USE_TIMER
   solve_timer->stop();
#endif
#ifdef CH_USE_TIMER
   shutdown_timer->start() ;
#endif

   if ( (m_plot_interval >= 0 || m_plot_time_interval >= 0.0) && (m_last_plot!=m_cur_step) ) {
      writePlotFile();
   }

   if ( m_checkpoint && (m_last_checkpoint!=m_cur_step) ) {
      writeCheckpointFile();
   }

   now = std::chrono::system_clock::now();
   std::time_t now_time = std::chrono::system_clock::to_time_t(now);
   double wall_time_sec = ((double) (now_time - m_main_start));
   double solve_walltime = ((double) (solve_end - m_solve_start));
   double wall_time_sec_g, solve_walltime_g;
#ifdef CH_MPI
   MPI_Allreduce(&wall_time_sec,&wall_time_sec_g,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   MPI_Allreduce(&solve_walltime,&solve_walltime_g,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
   wall_time_sec_g = wall_time_sec;
   solve_walltime_g = solve_walltime;
#endif
   if(!procID()) {
     cout << "Solve wall time (in seconds): " << solve_walltime_g << "\n";
     cout << "Total wall time (in seconds): " << wall_time_sec_g << "\n";
   }

#ifdef CH_USE_TIMER
   m_shutdown_timer->stop() ;
   m_all_timer->stop() ;
   Timer::TimerSummary();
#endif
}


void Simulation::parseParameters( ParmParse& a_ppsim )
{
   // This determines the amount of diagnositic output generated
   a_ppsim.query( "verbosity", m_verbosity );
   CH_assert( m_verbosity >= 0 );

   // Stop after this number of steps
   a_ppsim.query( "max_step", m_max_step );
   CH_assert( m_max_step >= 0 );

   // Stop when the simulation time get here
   a_ppsim.query( "max_time", m_max_time );
   CH_assert( m_max_time >= 0.0 );
   
   // get the wall time
   a_ppsim.query( "wall_time_hrs", m_max_wall_time_hrs );

   // If set, use as the fixed time step size
   if ( a_ppsim.query( "fixed_dt", m_fixed_dt ) ) {
      CH_assert( m_fixed_dt>0.0 );
      setFixedTimeStep( m_fixed_dt );
      m_adapt_dt = false;
   }

   // set cfl number for the case of dynamic timestep selection
   if ( a_ppsim.query( "cfl_number", m_cfl ) ) {
      CH_assert( m_cfl>0.0 && m_cfl<=2.0 );
      if (!m_adapt_dt) MayDay::Error( "fixed_dt and cfl are mutually exclusive!" );
   }
   a_ppsim.query("dt_check_interval", m_dt_check_interval);
   
   // set cfl number for scattering
   if ( a_ppsim.query( "cfl_scatter", m_cfl_scatter ) ) {
      CH_assert( m_cfl_scatter>0.0 && m_cfl_scatter<=2.0 );
      if (!m_adapt_dt) MayDay::Error( "fixed_dt and cfl_scatter are mutually exclusive!" );
   }
   a_ppsim.query("dt_scatter_interval", m_dt_scatter_interval);
 
   // Set up plot file writing
   a_ppsim.query( "plot_interval", m_plot_interval );
   a_ppsim.query( "plot_time_interval", m_plot_time_interval );
   if( m_plot_time_interval>0.0 ) m_plot_interval = 0;
   a_ppsim.query( "plot_prefix", m_plot_prefix );

   // get History parameter
   a_ppsim.query("history", m_history);
   a_ppsim.query("history_interval", m_history_interval );
   if(m_history) CH_assert( m_history_interval>0 );

   // Set up checkpoint file writing
   a_ppsim.query( "checkpoint_interval", m_checkpoint_interval );
   a_ppsim.query( "checkpoint_time_interval", m_checkpoint_time_interval );
   if( m_checkpoint_time_interval>0.0 ) m_checkpoint_interval = 0;
   if( m_checkpoint_time_interval>0.0 || m_checkpoint_interval > 0) m_checkpoint = true;
   a_ppsim.query( "checkpoint_prefix", m_checkpoint_prefix );
   
   // seed for the global random number generator
   a_ppsim.query( "fixed_random_seed", m_fixed_random_seed );
   if(m_fixed_random_seed>=0) m_fixed_random_seed = m_fixed_random_seed*(procID() + 1);
   MathUtils::seedRNG(m_fixed_random_seed);
   
}


void Simulation::setFixedTimeStep( const Real& a_dt )
{
   m_cur_dt = a_dt; 
}

void Simulation::enforceTimeStep(const Real  a_local_dt)
{
   // ensure time step is the same on all processors to machine precision
#ifdef CH_MPI
   MPI_Allreduce( &a_local_dt, &m_cur_dt, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD ); 
#endif

}


void Simulation::writeFiles()
{
   CH_TIMERS("Simulation::writeFiles()");
   
   if ( m_plot_time_interval>0.0 ) {
      if ( m_cur_time > m_plot_time - m_new_epsilon ) {
         writePlotFile();
         m_last_plot = m_cur_step;
         m_plot_time = m_plot_time + m_plot_time_interval;
      } 
   } 
   else {
      if ( (m_cur_step % m_plot_interval)==0 ) {
         writePlotFile();
         m_last_plot = m_cur_step;
      }
   }

   if(m_history) {
      if ( (m_cur_step % m_history_interval)==0 ) {
         writeHistFile(false);
         m_last_history = m_cur_step;
      }
   }
   
   if(m_checkpoint) {
      if ( m_checkpoint_time_interval>0.0 ) { 
         if ( m_cur_time > m_checkpoint_time - m_new_epsilon ) {
            writeCheckpointFile();
            m_last_checkpoint = m_cur_step;
            m_checkpoint_time = m_checkpoint_time + m_checkpoint_time_interval;
         }
      }
      else {
         if ( (m_cur_step % m_checkpoint_interval)==0 ) {
            writeCheckpointFile();
            m_last_checkpoint = m_cur_step;
         }
      }
   }

}

void Simulation::preTimeStep()
{
   CH_TIMERS("Simulation::preTimeStep()");

   // set the stable time step
   if ( m_cur_step==m_restart_step || ((m_cur_step+1) % m_dt_scatter_interval)==0 ) {
      if(m_system->useScattering()) m_dt_scatter = m_system->scatterDt( m_cur_step );
   }
   if ( m_cur_step==m_restart_step || ((m_cur_step+1) % m_dt_check_interval)==0 ) {
      if(m_system->useParts()) m_dt_parts = m_system->partsDt( m_cur_step );
      if(m_system->useSpecialOps()) m_dt_specialOps = m_system->specialOpsDt( m_cur_step );
   }
   if ( m_adapt_dt ) { 
      if(m_system->useParts()) m_cur_dt = std::min( m_dt_light*m_cfl, m_dt_parts*m_cfl );
      if(m_system->useScattering()) m_cur_dt = std::min( m_cur_dt, m_dt_scatter*m_cfl_scatter );
      if(m_system->useSpecialOps()) m_cur_dt = std::min( m_cur_dt, m_dt_specialOps );
      enforceTimeStep( m_cur_dt );
   } 
   
   // If less than a time step from the final time, adjust time step
   // to end just over the final time.
   Real timeRemaining = m_max_time - m_cur_time;
   if ( m_adapt_dt && m_cur_dt > timeRemaining ) {
      //m_cur_dt = timeRemaining + m_max_time * s_DT_EPS;
      m_cur_dt = timeRemaining;
   }
   
   m_new_epsilon = m_cur_dt/1000.0;

   m_system->preTimeStep( m_cur_time, m_cur_dt, m_cur_step );

}

void Simulation::postTimeStep()
{
  CH_TIMERS("Simulation::postTimeStep()");
  m_system->postTimeStep( m_cur_time, m_cur_dt, m_cur_step );
}

#include "NamespaceFooter.H"
