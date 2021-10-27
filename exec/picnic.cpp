#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ParmParse.H"
//#define CH_SPACE_DIM DIM

#include "CH_HDF5.H"
#include "parstream.H"

#include "DebugDump.H"
#include "memtrack.H"
#include "CH_Attach.H"
#include "FABView.H"
#include "parstream.H"

#include "Simulation.H"
#include "MathUtils.H"

#include <random>

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

/***************/

// set the mt19937 random generator as global variable 
std::mt19937 global_rand_gen;

int main(int a_argc, char* a_argv[])
{

#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  setChomboMPIErrorHandler();
#endif

   int rank, num_procs;
#ifdef CH_MPI
   MPI_Comm_rank(Chombo_MPI::comm, &rank);
   MPI_Comm_size(Chombo_MPI::comm, &num_procs);
#else
   rank = 0;
   num_procs = 1;
#endif

   if(rank==0) {
      cout << endl;
      cout << "PICNIC: number of procs = " << num_procs << endl;
      cout << "PICNIC: SpaceDim = " << SpaceDim << endl;
   }

   // Check for an input file
   char* inFile = NULL;

   if(a_argc > 1)
   {
      inFile = a_argv[1];
      if(rank==0) cout << "input file = " << inFile << endl << endl;
   }
   else
   {
      if(rank==0) cout << "Usage: <executable name> <inputfile>" << endl;
      if(rank==0) cout << "No input file specified" << endl;
      return -1;
   }

   // Parse the command line and the input file (if any)
   ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
  
   // initialize the simulation
   Simulation sim( pp ); // initialization done in constructor

   // run/finalize the simulation
   while ( sim.notDone() ) sim.advance();
   sim.finalize();

#ifdef CH_MPI
   CH_TIMER_REPORT();
   //dumpmemoryatexit();
   MPI_Finalize();
#endif

   return 0;

}


