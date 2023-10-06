#!/usr/bin/bash

MSS_DIR="/mss/halld/RunPeriod-2017-01/recon/ver03/REST"
RUN_LIST="2017_01-ver03.part.runList"
# MSS_DIR="/mss/halld/RunPeriod-2018-01/recon/ver02/REST"
# RUN_LIST="2018_01-ver02.part.runList"
# MSS_DIR="/mss/halld/RunPeriod-2018-08/recon/ver02/REST"
# RUN_LIST="2018_08-ver02.part.runList"

CMD_LINE="makeJcacheFileList.py"
while read RUN_NMB; do
  CMD_LINE="${CMD_LINE} ${MSS_DIR}/${RUN_NMB}"
done <"${RUN_LIST}"
echo "${CMD_LINE}"
eval "${CMD_LINE}"
