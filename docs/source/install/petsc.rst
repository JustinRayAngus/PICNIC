.. _install-petsc:

.. raw::html

   <style>
   .rst-content section>img {
       width: 30px;
       margin-bottom: 0;
       margin-top: 0;
       margin-right: 15px;
       margin-left: 15px;
       float: left;
   }
   </style>

.. only:: html

   .. image:: PETSc-TAO_RGB.svg

PETSc
=====

PICNIC employs several implicit time integrators that use JFNK. The PETSc library (https://petsc.org/release/) is used for both nonlinear (https://petsc.org/main/manual/snes/) and linear (https://petsc.org/main/manual/ksp/) solvers.

Obtaining PETSc
----

Checkout the main development branch of PETSc:

.. code-block:: bash

   git clone https://gitlab.com/petsc/petsc.git petsc

Set environment variables for PETSc
----

The following environment variables need to be set for PICNIC to use PETSc. PICNIC can work without PETSc, you just wont be able to use their efficient linear solvers. Simply to not define PETSC_DIR or undefine it prior to compiling PICNIC in order to compile without PETSc.

.. code-block:: bash

   export PETSC_DIR="path/to/petsc/"
   export PETSC_ARCH="arch-opt"

Configure PETSc
----

Configure PETSc with appropriate linear solvers needed by PICNIC and then build:

.. code-block:: bash

   cd $PETSC_DIR
   ./configure --with-batch --with-cc=mpicc --with-fc=mpif90 --with-cxx=mpicxx COPTFLAGS="-O2 -std=c99" FOPTFLAGS="-O2" CXXOPTFLAGS="-O2" --with-shared-libraries --with-debugging=0 --download-make --download-cmake --download-hypre --download-superlu --download-superlu_dist --download-parmetis --download-metis --with-cxx-dialect=C++11
   make all

