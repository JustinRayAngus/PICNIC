#!/bin/bash

COMPILE_OPT=TRUE
COMPILE_DEBUG=FALSE

: ${PICNIC_TEST_DIM:?"Please set PICNIC_TEST_DIM in your environment or manually above"}
MAKE_COMMAND="make -j all OPT=$COMPILE_OPT DEBUG=$COMPILE_DEBUG DIM=$PICNIC_TEST_DIM" 

# Install the Chombo makefile
sed -e "s:MASTER_MPICXX:"$MPICXX":" \
    -e "s:MASTER_CXX:"$CXX":" \
    -e "s:MASTER_CC:"$CC":" \
    -e "s:MASTER_FC:"$FC":" \
    -e "s:MASTER_HDF5_SERIAL_DIR:"$HDF5_SERIAL_DIR":" \
    -e "s:MASTER_HDF5_PARALLEL_DIR:"$HDF5_PARALLEL_DIR":" \
    < $PICNIC_TEST_PICNIC_DIR/config/Chombo_Make.defs.local_template > $PICNIC_TEST_CHOMBO_DIR/lib/mk/Make.defs.local

# Build PICNIC
cd $PICNIC_TEST_PICNIC_DIR/exec
make realclean 2>&1 > make_realclean.log
echo "Building PICNIC ($MAKE_COMMAND)."
make_logname="make_${PICNIC_TEST_DIM}d.log"
$MAKE_COMMAND &> $make_logname
echo "  Done. Check $PICNIC_TEST_PICNIC_DIR/exec/${make_logname} for errors, if any."

