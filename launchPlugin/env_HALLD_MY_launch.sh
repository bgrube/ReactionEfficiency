#!/bin/bash

# ALTERNATIVE HOME FOLDER
export HALLD_MY="/work/halld/bgrube/halld_my"
# SET FULL ENVIRONMENT
# source /group/halld/Software/build_scripts/gluex_env_jlab.sh /group/halld/www/halldweb/html/halld_versions/version_5.11.0.xml
# use private copy of halld-recon
source /group/halld/Software/build_scripts/gluex_env_jlab.sh /work/halld/bgrube/ProtonTrackEfficiency/ReactionEfficiency/launchPlugin/version.xml.halld_recon

# # quick fix for rootcling compilation
# unset CPLUS_INCLUDE_PATH

# #try geometry from ccdb
# export JANA_GEOMETRY_URL="ccdb:///GEOMETRY/main_HDDS.xml"

# use local file-based CDDB
export SQLITE_PATH="/work/halld/ccdb_sqlite/55/ccdb.sqlite"
export CCDB_CONNECTION="sqlite:///${SQLITE_PATH}"
export JANA_CALIB_URL="sqlite:///${SQLITE_PATH}"
