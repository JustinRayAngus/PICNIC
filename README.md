--------------------
PICNIC
--------------------

PICNIC is a Particle-In-Cell code with Nicely Incorporated Collisions.

The code is written in c++/FORTRAN and uses the Chombo framework.

--------------------
Compiling with PETSc
--------------------

See PETSc (https://gitlab.com/petsc/petsc) documentation
and repo for instructions on how to download and compile
it. A quick set of instructions is available here:
https://debog.github.io/codes/petsc.html

- PETSc's compilation requires setting the environment
variables PETSC_DIR. Make sure they are 
set to meaningful values.

- Compile PICNIC as usual - it will use the PETSC_DIR variable
to find PETSc and compile with the PETSc interface.
