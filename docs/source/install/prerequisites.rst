.. _install-prerequisites:

.. only::html

Prerequisites
=====

Building PICNIC requires pre-installing various compilers/programs and setting corresponding environment variables. PICNIC is written in C++/fortran and uses MPI. HDF5 is used for writting. The environment variables presented below can be placed in your .bashrc file. 

Compilers
----

If installing on the LC, mpi and the various compilers are pre-installed. However, the default versions may not be the correct ones. The correct versions can be loaded using:

.. code-block:: bash
   module load gcc/10.3.1-magic
   module load mvapich2/2.3.7

If installing on a Mac, fortran compilers and mpi can be installed using brew:

.. code-block:: bash

   brew install gcc
   brew install open-mpi

After installation, set the following environment variables for compilers:

.. code-block:: bash

   export CC=$(which gcc)
   export CXX=$(which g++)
   export CPP=$(which cpp)
   export F77=$(which gfortran)
   export FC=$(which gfortran)
   export MPICC=$(which mpicc)
   export MPICXX=$(which mpicxx)
   export MPIF77=$(which mpif77)
   export MPIF90=$(which mpif90)

When installing gcc on Mac OS using brew, I have found that the directory containing the fortran libraries (e.g., libgfortran.a) are in a non-trival location and one has to manually specify this location in the Make.defs.local file used by Chombo. This can be done using an environment variable. When I install gcc using brew on Mac with Sonoma 14.5, I use:

.. code-block:: bash

   export FORT_LIB_DIR="/opt/homebrew/Cellar/gcc/14.1.0_2/lib/gcc/current"

HDF5
----

PICNIC uses the HDF5 format for parallel writing of particles and grid quantities. On LC, set the following environment variables for HDF5:

.. code-block:: bash

   module load hdf5-parallel/$hdf5_version
   H5DIFF_SUFFIX=/bin/h5diff
   H5DIFF_PATH=$(which h5diff)
   HDF5_DIR=${H5DIFF_PATH%$H5DIFF_SUFFIX}
   export HDF5_DIR_SERIAL="${HDF5_DIR}"
   export HDF5_DIR_PARALLEL="${HDF5_DIR}"

If installing on a Mac, HDF5 can be intalled using brew. We need the mpi compatible version (which requires first having open-mpi following instructions above):

.. code-block:: bash

   brew install hdf5-mpi

The default configuration can be viewed using

.. code-block:: bash

   h5pcc -showconfig

This command will show you the installation point. For me, it is /opt/homebrew/Cellar/hdf5-mpi/1.14.3. Use this path for the environment variables:

.. code-block:: bash

   HDF5_DIR=/opt/homebrew/Cellar/hdf5-mpi/1.14.3
   export HDF5_DIR_SERIAL="${HDF5_DIR}"
   export HDF5_DIR_PARALLEL="${HDF5_DIR}"

