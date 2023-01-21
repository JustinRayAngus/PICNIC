#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <random>

#include "CH_HDF5.H"
#include "ParmParse.H"
#include "parstream.H"
#if 1  // warning, OS dependencies, will not work on all platforms
#include <sys/stat.h>
#endif

#include "DebugDump.H"
#include "memtrack.H"
#include "FABView.H"
#include "parstream.H"

#include "Simulation.H"
#include "MathUtils.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "UsingNamespace.H"

inline int checkCommandLineArgs( int a_argc, char* a_argv[] )
{  
   // Check for an input file
   if (a_argc<=1) {
      cout << "Usage: PICNIC...ex <inputfile>" << endl;
      cout << "No input file specified" << endl;
      return -1;
   }
   return 0;
}

#ifdef with_petsc
static const char help[] = "PICNIC";
#endif

// set the mt19937 random generator as global variable 
std::mt19937 global_rand_gen;

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
   // Start MPI
   MPI_Init( &a_argc, &a_argv );
   setChomboMPIErrorHandler();
#endif

#ifdef with_petsc
  PetscInitialize(&a_argc,&a_argv,(char*)0,help);
#endif

   int status = checkCommandLineArgs( a_argc, a_argv );

   if (status==0) {

      if(!procID()) {
         cout << endl;
         cout << "PICNIC: number of procs = " << numProc() << endl;
         cout << "PICNIC: SpaceDim = " << SpaceDim << endl;
#ifdef MASS_MATRIX_TEST
         cout << "MASS_MATRIX_TEST flag is defined" << endl;
#endif
         cout << "input file = " << a_argv[1] << endl << endl;
      }

      ParmParse pp( a_argc-2, a_argv+2, NULL, a_argv[1] );
      Simulation sim( pp ); 
      while ( sim.notDone() ) sim.advance();
      sim.finalize();
   }

#ifdef with_petsc
  PetscFinalize();
#endif

#ifdef CH_MPI
   CH_TIMER_REPORT();
   char* timerEnv = getenv("CH_TIMER");
   if(!procID() && timerEnv!=NULL) {
     mkdir( "time_tables", 0777 );
     std::system("mv time.* time_tables/");
   }
   MPI_Finalize();
#endif

   return status;
}
