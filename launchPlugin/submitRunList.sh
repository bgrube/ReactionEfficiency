#!/usr/bin/bash

PERIOD="2017_01-ver03"
# PERIOD="2018_01-ver02"
# PERIOD="2018_08-ver02"

RUN_LIST="${PERIOD}.part.runList"

swif2 create ReactionEfficiency_RD_${PERIOD}

while read RUN_NMB; do
  ./launch.orig.py --verbose=True jobs_ReactionEfficiency_RD.config "${RUN_NMB}" "${RUN_NMB}" &> "ReactionEfficiency_RD_${PERIOD}_${RUN_NMB}.log"
done <${RUN_LIST}

swif2 run -workflow ReactionEfficiency_RD_${PERIOD}
