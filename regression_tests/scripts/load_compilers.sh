#!/bin/bash

if [[ $HOSTNAME == "quartz"* ]]; then
  module load gcc/10.3.1-magic
  module load mvapich2/2.3.7
else
  if [[ $HOSTNAME == "surface"* ]]; then
    use gcc-4.7.4p
    use mvapich2-gnu-2.2
  fi
fi

export CC=$(which gcc)
export CXX=$(which g++)
export CPP=$(which cpp)
export F77=$(which gfortran)
export FC=$(which gfortran)
export MPICC=$(which mpicc)
export MPICXX=$(which mpicxx)
export MPIF77=$(which mpif77)
export MPIF90=$(which mpif90)

export HDF5_SERIAL_DIR=$HDF5_DIR_SERIAL
export HDF5_PARALLEL_DIR=$HDF5_DIR_PARALLEL

echo "CC is $CC"
echo "CXX is $CXX"
echo "CPP is $CPP"
echo "F77 is $F77"
echo "FC is $FC"
echo "MPICXX is $MPICXX"
echo "HDF5_SERIAL_DIR is $HDF5_SERIAL_DIR"
echo "HDF5_PARALLEL_DIR is $HDF5_PARALLEL_DIR"
echo ""
