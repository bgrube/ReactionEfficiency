#/usr/bin/bash

RUN_LIST="2018_01-ver02.part.runList"

swif2 create ReactionEfficiency_RD_2018_01-ver02

while read RUN_NMB; do
  ./launch.py --verbose=True jobs_ReactionEfficiency_RD.config "${RUN_NMB}" "${RUN_NMB}" &> "ReactionEfficiency_RD_2018_01-ver02_${RUN_NMB}.log"
done <${RUN_LIST}

swif2 run -workflow ReactionEfficiency_RD_2018_01-ver02
