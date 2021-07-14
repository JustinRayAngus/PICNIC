#!/bin/bash

: ${PICNIC_TEST_CHOMBO_DIR:?"Please set PICNIC_TEST_CHOMBO_DIR in your environment or manually above"}

if [ -d $PICNIC_TEST_CHOMBO_DIR ]; then
   echo "Updating Chombo. This may take a while."
   cd $PICNIC_TEST_CHOMBO_DIR
   svn up 2>&1 > chombo_svn_up.log
   echo "  Done. See $PICNIC_TEST_CHOMBO_DIR/chombo_svn_up.log."
   cd ..
else
   echo "Checking out Chombo. This may take a while."
   svn checkout https://anag-repo.lbl.gov/svn/Chombo/trunk $PICNIC_TEST_CHOMBO_DIR 2>&1 > chombo_svn_co.log
   echo "  Done. See chombo_svn_co.log."
fi
