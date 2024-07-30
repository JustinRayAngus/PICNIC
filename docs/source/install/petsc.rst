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

PICNIC employs several implicit time integrators that use JFNK. The PETSc library is used for both nonlinear (SNES) and linear solvers.

Obtaining PETSc
----

To checkout PETSc:

.. code-block:: bash

   git clone -b https://gitlab.com/petsc/petsc.git petsc

Set environment variables for PETSc
----

.. code-block:: bash

   export PETSC_DIR="path/to/petsc/"
   export PETSC_ARCH="arch-opt"

Configure PETSc to with appropriate linear solvers:

Configure PETSc
----

.. code-block:: bash

   cd $PETSC_DIR
   ./configure --with-batch --with-cc=mpicc --with-fc=mpif90 --with-cxx=mpicxx COPTFLAGS="-O2 -std=c99" FOPTFLAGS="-O2" CXXOPTFLAGS="-O2" --with-shared-libraries --with-debugging=0 --download-make --download-cmake --download-hypre --download-superlu --download-superlu_dist --download-parmetis --download-metis --with-cxx-dialect=C++11

