#!/usr/bin/env bash
# removes data files

FIT_DIR=${1:-"./BruFitOutput"}

for PATTERN in '*tree*.root' '*weights*.root' 'boot*.root'
do
  find "${FIT_DIR}" -iname "${PATTERN}" -exec rm -vf {} \;
done
