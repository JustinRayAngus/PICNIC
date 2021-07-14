#!/bin/bash

# user sets picnic and Chombo dirs used for testing
export PICNIC_TEST_PICNIC_DIR=/path/to/picnic
export PICNIC_TEST_CHOMBO_DIR=/path/to/Chombo

# Remember to set this to the location of the baselines
export PICNIC_TEST_BASELINES_LOC=$PWD/../baselines
export PICNIC_TEST_TOLERANCE=1.e-15

export PICNIC_TEST_SCRIPTS_PATH=$PICNIC_TEST_PICNIC_DIR/regression_tests/scripts

for i in {1..3}
do
  export PICNIC_TEST_DIM=$i
  export PICNIC_TEST_BASELINES_DIR=${PICNIC_TEST_BASELINES_LOC}/${i}d
  export PICNIC_TEST_INPUTS_DIR=$PICNIC_TEST_BASELINES_DIR
  if [ -d $PICNIC_TEST_BASELINES_DIR ]
  then
    echo ""
    echo "Running ${i}D tests."
    echo "Baselines are in $PICNIC_TEST_BASELINES_DIR"
    $PICNIC_TEST_SCRIPTS_PATH/build_picnic_simple.sh
    $PICNIC_TEST_SCRIPTS_PATH/test_solutions.sh "`ls $PICNIC_TEST_BASELINES_DIR`"
    echo "Done running ${i}D tests."
  fi
done
