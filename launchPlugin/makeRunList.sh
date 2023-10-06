#!/usr/bin/bash

MSS_DIR="/mss/halld/RunPeriod-2017-01/recon/ver03/REST"
RUN_LIST="2017_01-ver03.runList"
# MSS_DIR="/mss/halld/RunPeriod-2018-01/recon/ver02/REST"
# RUN_LIST="2018_01-ver02.runList"
# MSS_DIR="/mss/halld/RunPeriod-2018-08/recon/ver02/REST"
# RUN_LIST="2018_08-ver02.runList"

rm -f "${RUN_LIST}"
for RUN in $(find "${MSS_DIR}" -type d -name '??????')
do
  basename "${RUN}" >> "${RUN_LIST}"
done

# see https://unix.stackexchange.com/a/369204
# to take only every 5th line from file list run
# awk 'NR % 5 == 0' input > output
