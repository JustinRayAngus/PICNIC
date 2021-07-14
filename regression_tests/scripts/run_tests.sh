#!/bin/bash

ROOT_DIR=$PWD

export PICNIC_SOURCE_REPO="ssh://git@czgitlab.llnl.gov:7999/angus1/picnic.git"
export PICNIC_BRANCH="main"

export PICNIC_TEST_PICNIC_DIR=$ROOT_DIR/picnic
export PICNIC_TEST_CHOMBO_DIR=$ROOT_DIR/Chombo

export PICNIC_TEST_SCRIPTS_PATH=$PICNIC_TEST_PICNIC_DIR/regression_tests/scripts

# Remember to set this to the location of the baselines
export PICNIC_TESTS_BASELINES_LOC=$PWD/../baselines
export PICNIC_TEST_TOLERANCE=1.e-15

# Get PICNIC
./get_picnic.sh

# Load compilers
$PICNIC_TEST_SCRIPTS_PATH/load_compilers.sh
# Get Chombo
$PICNIC_TEST_SCRIPTS_PATH/get_chombo.sh

for i in {1..3}
do
  export PICNIC_TEST_DIM=$i
  export PICNIC_TEST_BASELINES_DIR=${PICNIC_TESTS_BASELINES_LOC}/${i}d
  export PICNIC_TEST_INPUTS_DIR=$PICNIC_TEST_BASELINES_DIR
  if [ -d $PICNIC_TEST_BASELINES_DIR ]
  then
    echo ""
    echo "Running ${i}D tests."
    echo "Baselines are in $PICNIC_TEST_BASELINES_DIR"
    $PICNIC_TEST_SCRIPTS_PATH/build_picnic.sh
    $PICNIC_TEST_SCRIPTS_PATH/test_solutions.sh "`ls $PICNIC_TEST_BASELINES_DIR`"
    echo "Done running ${i}D tests."
  fi
done
