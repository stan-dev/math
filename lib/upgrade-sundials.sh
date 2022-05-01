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
# This script will take these steps and add a git commit after each:
# 1. Remove the old version of Sundials.
# 2. Modify makefiles and README with new version numbers
# 3. Unpack the new Sundials version.
# 4. Prune Sundials.
# 5. Add cmake-generated config files to Sundials.
# 6. Modify Sundials by removing printf and fprintf (for CRAN)
# 7. Remove unused bits from Sundials which have compilation problems and are
#    not needed.
#
# This script should be run from the lib/ folder and the argument should be the
# name of the sundials-*.tar.gz file.
#===============================================================================
set -e
err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

usage() {
    cat <<HELP_USAGE

    Usage:
        $0 sundials-<version>.tar.gz

    Download the new Sundials version to sundials-<version>.tar.gz in this
    directory. Provide that filename as the only parameter to this script.

    This script will upgrade the Sundials library to the one provide,
    making modifications to the library as necessary. The modifications
    will be separate git commits for verification.

    If there are any unstaged modifications or staged modifications that
    haven't been git committed, this script will quit prior to adding any
    git commits.

HELP_USAGE
    exit 0
}

[ -z "$1" ] && { usage; }


sundials_filename=$1
## extract version number from argument
## ex: sundials-4.0.1.tar.gz -> version=4.0.1

sundials_version=`expr "$sundials_filename" : '^sundials-\(.*\)\.tar\.gz'`

if [[ -z "$sundials_version" ]] ; then
    cat <<SUNDIALS_VERSION_ERROR

    Can not parse the Sundials version from the filename.
    Expecting: sundials-<version>.tar.gz

    Is the filename correct? "$sundials_filename"

SUNDIALS_VERSION_ERROR
    exit 1
fi

## check git tree for modifications
if ! git diff-index --quiet HEAD --
then
    cat <<GIT_ERROR

    Please commit or stash any git changes prior to
    running this script. Quitting.

GIT_ERROR
    exit 1
fi

sundials_old_version=`expr "$(ls -d sundials_*/)" : '^sundials_\(.*\)/'`

echo "Old version:  $sundials_old_version"
echo "New version:  $sundials_version"

# 1. Remove the old version of Sundials.
git rm -r sundials_${sundials_old_version}/
rm -rf sundials_${sundials_old_version}/
git commit -m "upgrading to sundials v${sundials_version}; removing old sundials library"

# 2. Modify makefiles and README with new version number
sed -i -e "s|lib/sundials_${sundials_old_version}|lib/sundials_${sundials_version}|g" ../README.md ../make/*
sed -i -e "s|SUNDIALS (version ${sundials_old_version})|SUNDIALS (version ${sundials_version})|g" ../README.md
rm -f ../README.md*-e ../make/*-e
git add -u ../README.md ../make/*
git commit -m "upgrading to sundials v${sundials_version}; modifying with new version number"

# 3. Unpack the new Sundials version.
tar xvzf $sundials_filename
mv sundials-${sundials_version} sundials_${sundials_version}
git add sundials_${sundials_version}
git commit -m "upgrading to sundials v${sundials_version}; adding unmodified sundials library"

# 4. Prune Sundials.
cd sundials_${sundials_version}

# generate sundials_config.h and sundials_fconfig.h prior to removing these files
mkdir build
cd build
cmake ..
cd ..

git rm -rf cmake/ doc/ examples/ test/ */cvode */ida */arkode src/sundials/sundials_xbraid.c
find . -name CMakeLists.txt -exec git rm {} \;
git commit -m "upgrading to sundials v${sundials_version}; pruning files"

# 5. Add cmake-generated config files to Sundials.
cp build/include/sundials/*config.h include/sundials/
cp build/include/sundials/*export.h include/sundials/
# ... and  further files needed for compilation
cp -v ./src/sundials/sundials_debug.h include/

rm -rf build
git add include
git commit -m "upgrading to sundials v${sundials_version}; adding config files"

# 6. Modify Sundials by removing printf and fprintf (for CRAN)
cat >include/stan_sundials_printf_override.hpp <<EOF
#ifndef STAN_MATH_SUNDIALS_PRINTF_OVERRIDE_HPP
#define STAN_MATH_SUNDIALS_PRINTF_OVERRIDE_HPP

#ifdef WITH_SUNDIAL_PRINTF
#define STAN_SUNDIALS_PRINTF(...) printf(__VA_ARGS__)
#define STAN_SUNDIALS_FPRINTF(...) fprintf(__VA_ARGS__)
#else
#define STAN_SUNDIALS_PRINTF(...)
#define STAN_SUNDIALS_FPRINTF(...)
#endif

#endif
EOF

find src -name "*.c" -type f -exec sed -E -i _orig  's#[^sf]printf\(#STAN_SUNDIALS_PRINTF(#'g {} \;
find src -name "*.c" -type f -exec sed -E -i _orig  's#fprintf\(#STAN_SUNDIALS_FPRINTF(#'g {} \;
find src -name "*.c_orig" -exec rm {} \;
git add src include/stan_sundials_printf_override.hpp
git commit -m "upgrading to sundials v${sundials_version}; removing printf and fprintf for CRAN"

cat <<EOF





    Done upgrading Sundials from v${sundials_old_version} to v${sundials_version}.

    Please check the upgrade worked by running a test with CVODES linked.
    Example (from Math home directory):
      ./runTests.py test/unit/math/rev/functor/sho_ode_typed_test.cpp

EOF
