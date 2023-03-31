#!/usr/bin/env bash

FIT_DIR="./BruFitOutput"

rm -rf "${FIT_DIR}"
./fitMissingMassSquared.py &> "${FIT_DIR}/fitMissingMassSquared.log"

./cleanFitDir.sh "${FIT_DIR}"
./cleanBruFitLogFile.py "${FIT_DIR}"
./plotFitResults.py "${FIT_DIR}" &> "${FIT_DIR}/plotFitResults.log"
./plotEfficiencies.py "${FIT_DIR}" &> "${FIT_DIR}/plotEfficiencies.log"
