#!/usr/bin/bash

# write stdout and stderr to logfile
LOG_FILE="./hd_root.log"
exec > "${LOG_FILE}"
exec 2>&1

source env_HALLD_MY_launch.sh
echo "HALLD_RECON_HOME    = '${HALLD_RECON_HOME}'"
echo "HALLD_RECON_VERSION = '${HALLD_RECON_VERSION}'"
echo "PATH = '${PATH}'"
type -a hd_root

echo "-------------------------------------------------------------------------------"
hd_root -Pprint --config=jana_ReactionEfficiency_MC.config /cache/halld/gluex_simulations/REQUESTED_MC/S2018_bggen_ver02_23_batch01_20220103074512pm/hddm/dana_rest_bggen_042559_000.hddm
