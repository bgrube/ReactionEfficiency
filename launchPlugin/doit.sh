#!/usr/bin/bash

# write stdout and stderr to logfile
LOG_FILE="./hd_root.log"
exec > "${LOG_FILE}"
exec 2>&1

echo "HALLD_RECON_HOME    = '${HALLD_RECON_HOME}'"
echo "HALLD_RECON_VERSION = '${HALLD_RECON_VERSION}'"  # needs at least version 5.15.1
echo "CCDB_CONNECTION     = '${CCDB_CONNECTION}'"
echo "JANA_CALIB_URL      = '${JANA_CALIB_URL}'"
echo "PATH = '${PATH}'"
type -a hd_root

echo "-------------------------------------------------------------------------------"
hd_root --loadconfigs ./jana_ReactionEfficiency_MC.config ./work/osgpool/halld/REQUESTEDMC_OUTPUT/2019-11_bggen2026_batch01_4906/hddm/dana_rest_bggen_071854_000.hddm
