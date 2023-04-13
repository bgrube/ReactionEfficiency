#!/usr/bin/env bash

FIT_DIR="./BruFitOutput"

rm -rfv "${FIT_DIR}"
mkdir -pv "${FIT_DIR}"
echo "Starting fits..."
./fitMissingMassSquared.py "${FIT_DIR}" &> "${FIT_DIR}/fitMissingMassSquared.log"
RET=${?}
if [ ${RET} -ne 0 ]
then
  echo "Fitting script failed with exit code '${RET}'; exiting"
  exit 1
fi

./cleanFitDir.sh "${FIT_DIR}"
./cleanBruFitLogFile.py "${FIT_DIR}"
echo "Plotting fit results..."
./plotFitResults.py "${FIT_DIR}" &> "${FIT_DIR}/plotFitResults.log"
echo "Plotting efficiencies..."
./plotEfficiencies.py "${FIT_DIR}" &> "${FIT_DIR}/plotEfficiencies.log"
