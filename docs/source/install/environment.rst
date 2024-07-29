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

Libraries
----

The following environment variables are need for Chombo and PETSc:

.. code-block:: bash

   export CHOMBO_DIR="path/to/chombo/"
   export PETSC_DIR="path/to/petsc/"
   export PETSC_ARCH="arch-opt"
