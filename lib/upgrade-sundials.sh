#!/bin/bash
#===============================================================================
#
# This script is for the Math library maintainers to use upgrade Sundials.
# To update to a new Sundials version, download the sundials-<version>.tar.gz to
# the lib/ directory.
#   https://computation.llnl.gov/projects/sundials/sundials-software
# This script will also manage the git process so the different steps are in
# separate git commits.
#
# This script will:
# 1. Remove the old version of Sundials and add a new git commit.
# 2. Unpack the new Sundials version and add a new git commit.
# 3. Prune Sundials and add a new git commit.
# 4. Modify Sundials so it passes CRAN and add a new git commit.
#
# This script should be run from the lib/ folder and the argument should be the
# name of the sundials-*.tar.gz file.
#===============================================================================

usage() {
    cat <<HELP_USAGE

    Usage:
        $0 sundials-<version>.tar.gz

    Download the new Sundials version to sundials-<version>.tar.gz in this
    directory. Provide that filename as the only parameter to this script.

    This script will upgrade the Sundials library to the one provide,
    making modifications to the library as necessary. The modifications
    will be separate git commits for verification.

HELP_USAGE
    exit 0
}

[ -z "$1" ] && { usage; }

SUNDIALS_VERSION="4.0.1"

SUNDIALS_TAR=sundials-${SUNDIALS_VERSION}.tar.gz

tar xvzf $SUNDIALS_TAR
mv sundials-${SUNDIALS_VERSION} sundials_${SUNDIALS_VERSION}

cd sundials_${SUNDIALS_VERSION}

rm -rf INSTALL_GUIDE.pdf config/ doc/ examples/ test/ */cvode */ida */arkode */kinsol
find . -name CMakeLists.txt -exec rm {} \;

find src -name "*.c" -type f -exec sed -E -i _orig  's#[^sf]printf\(#STAN_SUNDIALS_PRINTF(#'g {} \;
find src -name "*.c" -type f -exec sed -E -i _orig  's#fprintf\(#STAN_SUNDIALS_FPRINTF(#'g {} \;

find src -name "*.c_orig" -exec rm {} \;

echo "In case of a large version change please consider updating sundials_config.h in sundials-config/include/sundials"
