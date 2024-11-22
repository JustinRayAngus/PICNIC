: ${PICNIC_TEST_PICNIC_DIR:?"Please set PICNIC_TEST_PICNIC_DIR in your environment or manually above"}
: ${PICNIC_TEST_INPUTS_DIR:?"Please set PICNIC_TEST_INPUTS_DIR in your environment or manually above"}

# Get the PICNIC executable name
picnic_executable=$(ls $PICNIC_TEST_PICNIC_DIR/exec/picnic${PICNIC_TEST_DIM}d*.ex)
#picnic_executable=$(ls $PICNIC_TEST_PICNIC_DIR/exec_relativistic/picnic${PICNIC_TEST_DIM}d*.ex)
echo "PICNIC executable is $picnic_executable"

# create dir for baseline solutions
DIRNAME=${PICNIC_TEST_DIM}d
rm -rf $DIRNAME && mkdir $DIRNAME
cd $DIRNAME

# Loop over the test problems
for prob_dir in $1; do

   # Copy the baseline; determine the number of processors and timelimit
   mkdir $prob_dir
   cp -r $PICNIC_TEST_INPUTS_DIR/$prob_dir/* $prob_dir
   cd $prob_dir
   rm -rf mesh_data particle_data history.txt time.table* stdout pout*
   input_file=$(ls *.in)
   np=`grep TEST $input_file | sed -e 's/.*\(np=\)/\1/' -e 's/\(,\).*//' -e 's/np=//'`
   timelimit=`grep TEST $input_file | sed -e 's/.*\(timelimit=\)/\1/' -e 's/\(,\).*//' -e 's/timelimit=//'`

   # Submit the run
   export CH_OUTPUT_INTERVAL=${np}
   jobcmd="srun -n${np} -ppdebug -t ${timelimit} ${picnic_executable} ${input_file} 2>&1 > stdout"
   echo "Running test problem $prob_dir on $np processors"
   echo "  input file is $input_file"
   eval $jobcmd &> srun.log
   echo "  done. see $PWD for simulation results."
   cd ..
done

echo "All done. Bye."
cd ..
