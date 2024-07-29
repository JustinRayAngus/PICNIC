.. _install-environment:

.. only::html

Environment variables
=====

Building PICNIC requires setting environment variables for various compilers and paths. The environment variables below can be placed in your .bashrc file. 

Compilers
----

Set the following environment variables for compilers:

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

If installing on LC, these mpi and fortran compilers should be pre-installed. If installing on a personal mac, these compilers can be intalled using brew:

.. code-block:: bash

   brew install gcc
   brew install open-mpi

HDF5
----

PICNIC uses the HDF5 format for particles and grid quantities. On LC, set the following environment variables for HDF5:

.. code-block:: bash

   H5DIFF_SUFFIX=/bin/h5diff
   module load hdf5-parallel/$hdf5_version
   H5DIFF_PATH=$(which h5diff)
   HDF5_DIR=${H5DIFF_PATH%$H5DIFF_SUFFIX}
   export HDF5_DIR_SERIAL="${HDF5_DIR}"
   export HDF5_DIR_PARALLEL="${HDF5_DIR}"

If installing on a personal mac, hdf5 can be intalled using brew. We need the mpi compatible version (which requires first having open-mpi following instructions above):

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

Libraries
----

The following environment variables are need for Chombo and PETSc:

.. code-block:: bash

   export CHOMBO_DIR="path/to/chombo/"
   export PETSC_DIR="path/to/petsc/"
   export PETSC_ARCH="arch-opt"
