#!/bin/bash

# user specifies paths to picnic and Chombo directories
export PICNIC_TEST_PICNIC_DIR=/path/to/picnic
export PICNIC_TEST_CHOMBO_DIR=/path/to/Chombo

export PICNIC_TEST_SCRIPTS_PATH=$PICNIC_TEST_PICNIC_DIR/regression_tests/scripts

for i in {1..3}
do
  export PICNIC_TEST_DIM=$i
  export PICNIC_TEST_INPUTS_DIR=$PICNIC_TEST_PICNIC_DIR/regression_tests/${i}d
  if [ -d $PICNIC_TEST_INPUTS_DIR ]
  then
    echo ""
    echo "Generating baselines for ${i}D tests."
    $PICNIC_TEST_SCRIPTS_PATH/build_picnic_simple.sh
    $PICNIC_TEST_SCRIPTS_PATH/generate_solutions.sh "`ls $PICNIC_TEST_INPUTS_DIR`"
    echo "Done generating baselines for ${i}D tests."
  fi
done
echo ""
