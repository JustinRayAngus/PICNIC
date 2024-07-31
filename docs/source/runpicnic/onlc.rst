.. _runpicnic-onlc:

.. raw::html

Run PICNIC on LC
=====

PICNIC can be run the CPU-based LC clusters Ruby and Dane. To see information about these machines and limits for jobs on their debug and batch queues, one can type

.. code-block:: bash

   news job.lim.dane
   news job.lim.ruby

Directly from the command line
----

The simplest way to run PICNIC on the LC is directly from the command line. This works for small jobs that use a small number of processors and run for a short amount of time. In a directory with a valid PICNIC input file named physics.in, run using

.. code-block:: bash

   srun-n4 -ppdebug path/to/picnic/exec/picnic.ex physics.in

where path/to/picnic.ex is the path to the appropriate picnic executable. The verbose information from the simulation will print directly to the screen during run time. To run in the backgroud and have the simulation information write to a text file, run using

.. code-block:: bash

   srun -n4 -ppdebug path/to/picnic.ex physics.in >> log.txt &

Running using a job script
----

For larger jobs, it is best to submit them to the appropriate queue using a job script. An example job script for running with 8 MPI ranks on the debug queue is

.. code-block:: bash

   #!/bin/csh 
   #SBATCH -t 0:10:00
   #SBATCH -J physics
   #SBATCH -N 1
   #SBATCH -n 8
   #SBATCH -A fipic
   #SBATCH -p pdebug

   setenv DECK   physics.in
   setenv PICNIC path/to/picnic/exec/picnic1d.Linux.64.mpicxx.mpif90.OPT.MPI.ex
   setenv NTASKS 8

   touch log.txt
   srun -n $NTASKS $PICNIC $DECK >>& log.txt

Assuming the above information is in a file called debug.job, submit the job using

.. code-block:: bash

   sbatch debug.job

One can see where the job is in the queue using

.. code-block:: bash

   squeue -ppdebug

