--------------------
PICNIC
--------------------

PICNIC is a Particle-In-Cell code with Nicely Incorporated Collisions.

The code is written in c++/FORTRAN and uses the Chombo framework.

--------------------
Obtaining Chombo
--------------------

PICNIC uses the Chombo library for data containers and efficient MPI-handling.
Chombo is developed and maintained by LBL-ANAG. Information about Chombo and
instructions for dowloading can be found at https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations

--------------------
Compiling with PETSc
--------------------

See PETSc (https://gitlab.com/petsc/petsc) documentation
and repo for instructions on how to download and compile
it. A quick set of instructions is available here:
https://debog.github.io/codes/petsc.html

- PETSc's compilation requires setting the environment
variables PETSC_DIR and PETSC_ARCH. Make sure they are
set to meaningful values.

- Compile PICNIC as usual - it will use the PETSC_DIR variable
to find PETSc and compile with the PETSc interface.

Note: You can compile and run without PETSc. Just do not define
the PETSC_DIR environment variable.

--------------------
Build instructions
--------------------

See picnic/docs/picnic_user_manual_gitRepo.pptx

--------------------
Release
--------------------

LLNL-CODE-862459

