#!/usr/bin/bash

PERIOD="2017_01-ver03"
# PERIOD="2018_01-ver02"
# PERIOD="2018_08-ver02"

RUN_LIST="../data/RD/${PERIOD}/${PERIOD}.part.runList"

if ! swif2 create "ReactionEfficiency_RD_${PERIOD}"
then
  echo "Failed to create workflow. Aborting."
  exit 1
fi

while read RUN_NMB
do
  "${HD_UTILITIES_HOME}/launch_scripts/launch/launch.py" --verbose=True ./jobs_ReactionEfficiency_RD.config "${RUN_NMB}" "${RUN_NMB}" &> "ReactionEfficiency_RD_${PERIOD}_${RUN_NMB}.log"
done <${RUN_LIST}

swif2 run -workflow "ReactionEfficiency_RD_${PERIOD}"
