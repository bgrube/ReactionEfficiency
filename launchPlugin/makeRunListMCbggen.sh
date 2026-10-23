#!/usr/bin/bash

MSS_DIR="/mss/halld/gluex_simulations/REQUESTED_MC/2019-11_bggen2026_batch01_4906/hddm"
RUN_LIST="2019_11-ver01.MCbggen.runList"

rm -f "${RUN_LIST}"
for FILE in "${MSS_DIR}"/*.hddm.tar
do
  FILE_NAME=$(basename "${FILE}" ".hddm.tar")
  RUN_NUMBER=${FILE_NAME/dana_rest_bggen_/}
  echo "${RUN_NUMBER}" >> "${RUN_LIST}"
done

# see https://unix.stackexchange.com/a/369204
# to take only every 5th line from file list run
# awk 'NR % 5 == 0' input > output
