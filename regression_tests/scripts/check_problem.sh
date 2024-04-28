#!/bin/bash

# This script compares the HDF5 files produced by a
# specifed test problem against its baseline using
# a prescribed absolute tolerance

# Name of the test problem to be checked
prob_dir=$1

PICNIC_root=$2
Chombo_root=$3

# Baseline directory.  The test name passed in $1 must
# match a directory in the baseline directory
base_dir=$4

# Location of h5diff
diff_prog=$5

# Absolute tolerance used for comparisons
#atol=$6

# Relative tolerance used for comparisons
rtol=$6

if [ -d "$base_dir/$prob_dir" ]; then
   echo "Checking $prob_dir"

   # Loop over the plot file subdirectories
   m=0
   success=1
   for base_plot_dir in $base_dir/$prob_dir/*; do
      if [ -d "$base_plot_dir" ]; then
         for base_file in $base_plot_dir/*.h5; do
            if [ -f "$base_file" ]; then
               plot_dir=`echo $base_plot_dir | sed 's/.*\(\/\)//'`
               file=`echo $base_file | sed 's/.*\(\/\)//'`
               #out=`$diff_prog --delta=$atol ${PICNIC_TEST_DIM}d/$prob_dir/$plot_dir/$file $base_file`
               out=`$diff_prog --relative=$rtol ${PICNIC_TEST_DIM}d/$prob_dir/$plot_dir/$file $base_file`
               if [ "$out" != "" ]; then
                  echo $file\:
                  echo "   " $out
                  success=0
               fi
               ((m++))
            fi
         done
         for base_plot_subdir in $base_plot_dir/*; do
            if [ -d "$base_plot_subdir" ]; then
               for base_file in $base_plot_subdir/*.h5; do
                  if [ -f "$base_file" ]; then
                     plot_dir=`echo $base_plot_dir | sed 's/.*\(\/\)//'`
                     plot_subdir=`echo $base_plot_subdir | sed 's/.*\(\/\)//'`
                     file=`echo $base_file | sed 's/.*\(\/\)//'`
                     #out=`$diff_prog --delta=$atol ${PICNIC_TEST_DIM}d/$prob_dir/$plot_dir/$plot_subdir/$file $base_file`
                     out=`$diff_prog --relative=$rtol ${PICNIC_TEST_DIM}d/$prob_dir/$plot_dir/$plot_subdir/$file $base_file`
                     if [ "$out" != "" ]; then
                        echo $file\:
                        echo "   " $out
                        success=0
                     fi
                     ((m++))
                  fi
               done
            fi
         done
      fi
   done
   if [ $m == 0 ]; then
      echo Found no files to compare
   else
      if (( success==1 )); then
         echo All $m file comparisons were successful
      fi
   fi
else
   echo "$base_dir/$prob_dir doesn't exist"
fi

cd ..


