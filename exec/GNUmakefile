# -*- Mode: Makefile -*- 

# trace the chain of included makefiles

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = picnic

# the location of the Chombo "lib" directory
#CHOMBO_HOME = ../../Chombo_3.2/lib
CHOMBO_HOME = ${CHOMBO_DIR}/lib

# names of Chombo libraries needed by this program, in order of search.
#
#LibNames = BaseTools
LibNames = BaseTools BoxTools AMRElliptic AMRTools
#LibNames = BaseTools AMRElliptic AMRTimeDependent AMRTools BoxTools HOAMRTools
#LibNames = BaseTools AMRElliptic AMRTimeDependent AMRTools BoxTools ParticleTools HOAMRTools

# the locations of the source code directories
base_dir = .
src_dirs = ../src
src_dirs += ../src/mesh
src_dirs += ../src/ibc
src_dirs += ../src/ibc/boundaries
src_dirs += ../src/ibc/grid_functions
src_dirs += ../src/ibc/time_functions
src_dirs += ../src/core
src_dirs += ../src/utils
src_dirs += ../src/particle_tools
src_dirs += ../src/species/pic
src_dirs += ../src/special_operators
src_dirs += ../src/scattering
src_dirs += ../src/fields
src_dirs += ../src/poisson
src_dirs += ../src/solver/linear_solvers
src_dirs += ../src/solver/nonlinear_solvers
src_dirs += ../src/solver/preconditioner
src_dirs += ../src/solver/utils
src_dirs += ../src/time

# input file for 'run' target
#INPUT = inputs

#########################################################################

# define picnic specific compile time flags here
#XTRACPPFLAGS += -DMASS_MATRIX_TEST
#XTRACPPFLAGS += -DMASS_MATRIX_COST_TEST
#XTRACPPFLAGS += -DRELATIVISTIC_PARTICLES

# if PETSc exists and PETSC_DIR is set in the environment,
# compile with PETSc

ifdef PETSC_DIR

ifneq ("$(wildcard $(PETSC_DIR))","")
$(info Found PETSc installation ($(PETSC_DIR)). Compiling with PETSc interface.)
XTRACPPFLAGS += -Dwith_petsc
include $(PETSC_DIR)/lib/petsc/conf/variables
XTRACPPFLAGS += $(PETSC_CC_INCLUDES)
XTRALIBFLAGS += $(PETSC_LIB) 
else
$(info Env var PETSC_DIR is set to $(PETSC_DIR); PETSc not found at this location. Compiling without PETSc interface.)
endif
else
$(info Environment var PETSC_DIR not set. Compiling without PETSc interface.)
endif

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets

