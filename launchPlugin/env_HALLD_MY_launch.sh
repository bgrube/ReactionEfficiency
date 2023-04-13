#!/bin/bash

# ALTERNATIVE HOME FOLDER
export HALLD_MY="/work/halld/bgrube/halld_my"
# SET FULL ENVIRONMENT
source /group/halld/Software/build_scripts/gluex_env_jlab.sh /group/halld/www/halldweb/html/halld_versions/version_5.7.1.xml

# # quick fix for rootcling compilation
# unset CPLUS_INCLUDE_PATH

# #try geometry from ccdb
# export JANA_GEOMETRY_URL="ccdb:///GEOMETRY/main_HDDS.xml"
