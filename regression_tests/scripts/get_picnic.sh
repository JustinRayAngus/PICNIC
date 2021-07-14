#!/bin/bash

echo ""
echo "PICNIC source is $PICNIC_SOURCE_REPO ($PICNIC_BRANCH branch)"
: ${PICNIC_TEST_PICNIC_DIR:?"Please set PICNIC_TEST_PICNIC_DIR in your environment or manually above"}

# delete PICNIC if it exists, then check out a new copy
if [ -d $PICNIC_TEST_PICNIC_DIR ]; then
  rm -rf $PICNIC_TEST_PICNIC_DIR
fi
echo "Checking out PICNIC (${PICNIC_BRANCH} branch)..."
git clone --progress $PICNIC_SOURCE_REPO $PICNIC_TEST_PICNIC_DIR &> picnic_git.log
echo "   Done. See picnic_git.log."
cd $PICNIC_TEST_PICNIC_DIR
git checkout $PICNIC_BRANCH
cd ../
echo ""
