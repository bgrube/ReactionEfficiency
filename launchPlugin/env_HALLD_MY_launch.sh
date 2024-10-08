#!/bin/bash

# ALTERNATIVE HOME FOLDER
export HALLD_MY="/work/halld/bgrube/halld_my"
# SET FULL ENVIRONMENT
# private plugins should be compiled using the SConstruct file generated by `mkfactory_plugin <plugin name>`
# running `scons -j10 install` in the plugin directory will install the plugin into `${HALLD_MY}/Linux_<OS and gcc version>` which takes precedence over the standard installation directories
source /group/halld/Software/build_scripts/gluex_env_jlab.sh /group/halld/www/halldweb/html/halld_versions/version.xml
# use private fully compiled copy of halld-recon
# source /group/halld/Software/build_scripts/gluex_env_jlab.sh /work/halld/bgrube/ProtonTrackEfficiency/ReactionEfficiency/launchPlugin/version.xml.halld_recon

# # quick fix for rootcling compilation
# unset CPLUS_INCLUDE_PATH

# #try geometry from ccdb
# export JANA_GEOMETRY_URL="ccdb:///GEOMETRY/main_HDDS.xml"

# use local file-based CDDB
# see https://halldweb.jlab.org/wiki/index.php/SQLite-form_of_the_CCDB_database
# export SQLITE_PATH="/work/halld/ccdb_sqlite/55/ccdb.sqlite"
export SQLITE_PATH="/group/halld/www/halldweb/html/dist/ccdb.sqlite"
export CCDB_CONNECTION="sqlite:///${SQLITE_PATH}"
export JANA_CALIB_URL="sqlite:///${SQLITE_PATH}"
