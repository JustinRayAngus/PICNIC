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

- PICNIC requires setting the environment variables CHOMBO_DIR. It should be like this:

`export CHOMBO_DIR=path/to/chombo`

- No need to configure or install Chombo. That will be done at an initial compilation of PICNIC


--------------------
Compiling with PETSc
--------------------

See PETSc (https://gitlab.com/petsc/petsc) documentation
and repo for instructions on how to download and compile
it. A quick set of instructions is available here:
https://debog.github.io/codes/petsc.html
See picnic/docs/picnic_user_manual_gitRepo.pdf for the details.

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

- Setup environment variables according to picnic/docs/picnic_user_manual_gitRepo.pdf

- Go to **path_to_picnic/exec** directory;

- Check the **Make.defs.local** file, adjust the main compiler/linker options to match your system; 

- Copy the **Make.defs.local** file into **path_to_chombo/lib/mk**;

- Run the following command to build PICNIC:

1D: `make -j all MPI=TRUE DEBUG=FALSE OPT=TRUE DIM=1`

2D: `make -j all MPI=TRUE DEBUG=FALSE OPT=TRUE DIM=2`


--------------------
Release
--------------------

LLNL-CODE-862459

