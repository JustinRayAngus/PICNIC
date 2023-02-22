#!/bin/bash

: ${PICNIC_TEST_PICNIC_DIR:?"Please set PICNIC_TEST_PICNIC_DIR in your environment or manually above"}
: ${PICNIC_TEST_CHOMBO_DIR:?"Please set PICNIC_TEST_CHOMBO_DIR in your environment or manually above"}
: ${PICNIC_TEST_BASELINES_DIR:?"Please set PICNIC_TEST_BASELINES_DIR in your environment or manually above"}
: ${PICNIC_TEST_TOLERANCE:?"Please set PICNIC_TEST_TOLERANCE in your environment or manually above"}

# Location of h5diff
diff_prog=$(which h5diff)

timestamp=`date | sed -e 's/ /_/g' -e 's/:/./g'`
test_outputs=${PICNIC_BRANCH}_${timestamp}
rm -rf $test_outputs
mkdir $test_outputs
echo "Entering directory $test_outputs"
cd $test_outputs

case_list=$(ls $PICNIC_TEST_BASELINES_DIR/)
$PICNIC_TEST_SCRIPTS_PATH/generate_solutions_salloc.sh " $case_list "

echo "Checking results..."

# Loop over test problems again and check them
for prob_dir in $case_list; do
   if [ -d "${PICNIC_TEST_DIM}d/${prob_dir}" ]; then
      echo "  Checking ${PICNIC_TEST_DIM}d/${prob_dir}."
      echo ------------------------------- >> report
      $PICNIC_TEST_SCRIPTS_PATH/check_problem.sh ${prob_dir} $PICNIC_TEST_PICNIC_DIR $PICNIC_TEST_CHOMBO_DIR $PICNIC_TEST_BASELINES_DIR $diff_prog $PICNIC_TEST_TOLERANCE >> report
   else
      echo "${PICNIC_TEST_DIM}d/${prob_dir} doesn't exist."
   fi
done
echo ------------------------------- >> report

# Email the results
mail -s "PICNIC test results" $USER < report

