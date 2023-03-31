#!/usr/bin/env bash

FIT_DIR="./BruFitOutput"

rm -rf "${FIT_DIR}"
mkdir -p "${FIT_DIR}"
./fitMissingMassSquared.py "${FIT_DIR}" &> "${FIT_DIR}/fitMissingMassSquared.log"
RET=${?}
if [ ${RET} -ne 0 ]
then
  echo "Fitting script failed with exit code '${RET}'; exiting"
  exit 1
fi

./cleanFitDir.sh "${FIT_DIR}"
./cleanBruFitLogFile.py "${FIT_DIR}"
./plotFitResults.py "${FIT_DIR}" &> "${FIT_DIR}/plotFitResults.log"
./plotEfficiencies.py "${FIT_DIR}" &> "${FIT_DIR}/plotEfficiencies.log"
