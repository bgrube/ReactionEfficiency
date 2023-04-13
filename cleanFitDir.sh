#!/usr/bin/env bash
# removes data files

FIT_DIR=${1:-"./BruFitOutput"}

find "${FIT_DIR}" -iname '*tree*.root' -exec rm -vf {} \;
find "${FIT_DIR}" -iname '*weights*.root' -exec rm -vf {} \;
