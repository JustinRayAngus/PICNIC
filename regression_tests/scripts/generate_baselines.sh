#!/bin/bash

ROOT_DIR=$PWD

export PICNIC_SOURCE_REPO="ssh://git@czgitlab.llnl.gov:7999/angus1/picnic.git"
export PICNIC_BRANCH="main"

export PICNIC_TEST_PICNIC_DIR=$ROOT_DIR/picnic
export PICNIC_TEST_CHOMBO_DIR=$ROOT_DIR/Chombo

export PICNIC_TEST_SCRIPTS_PATH=$PICNIC_TEST_PICNIC_DIR/regression_tests/scripts

# Get PICNIC
./get_picnic.sh

# Load compilers
$PICNIC_TEST_SCRIPTS_PATH/load_compilers.sh
# Get Chombo
$PICNIC_TEST_SCRIPTS_PATH/get_chombo.sh

for i in {1..3}
do
  export PICNIC_TEST_DIM=$i
  export PICNIC_TEST_INPUTS_DIR=$PICNIC_TEST_PICNIC_DIR/regression_tests/${i}d
  if [ -d $PICNIC_TEST_INPUTS_DIR ]
  then
    echo ""
    echo "Generating baselines for ${i}D tests."
    $PICNIC_TEST_SCRIPTS_PATH/build_picnic.sh
    $PICNIC_TEST_SCRIPTS_PATH/generate_solutions.sh "`ls $PICNIC_TEST_INPUTS_DIR`"
    echo "Done generating baselines for ${i}D tests."
  fi
done
echo ""
