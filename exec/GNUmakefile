# -*- Mode: Makefile -*- 

# trace the chain of included makefiles

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = myPIC

export MAXDIM=2
export DIM=2
CH_SPACEDIM = 2

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../Chombo_3.2/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = BaseTools AMRElliptic AMRTimeDependent AMRTools BoxTools HOAMRTools
//LibNames = BaseTools AMRElliptic AMRTimeDependent AMRTools BoxTools ParticleTools HOAMRTools

# the locations of the source code directories
base_dir = .
src_dirs = ../src
src_dirs += ../src/mesh
src_dirs += ../src/core
src_dirs += ../src/utils
src_dirs += ../src/particle_tools
src_dirs += ../src/state_variables/species/pic

# input file for 'run' target
INPUT = inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
