#!/bin/bash
#===============================================================================
#
# This script is for the Math library maintainers to use upgrade TBB.
# To update to a new TBB version, download the Source code tar tbb-<version>.tar.gz to
# the lib/ directory.
#   https://github.com/intel/tbb/releases
# This script will also manage the git process so the different steps are in
# separate git commits.
#
# This script will take these steps and add a git commit after each:
# 1. Remove the old version of TBB.
# 2. Modify makefiles and README with new version numbers
# 3. Unpack the new TBB version.
# 4. Prune the TBB folder
#
# This script should be run from the lib/ folder and the argument should be the
# in the form tbb-<version>.gz file.
#===============================================================================
set -e
err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

usage() {
    cat <<HELP_USAGE

    Usage:
        $0 tbb-<version>.tar.gz

    Download the new TBB version to tbb-<version>.tar.gz in this
    directory. Provide that filename as the only parameter to this script.

    This script will upgrade the TBB library to the one provide,
    making modifications to the library as necessary. The modifications
    will be separate git commits for verification.

    If there are any unstaged modifications or staged modifications that
    haven't been git committed, this script will quit prior to adding any
    git commits.

HELP_USAGE
    exit 0
}

[ -z "$1" ] && { usage; }


tbb_filename=$1
## extract version number from argument
## ex: 2019_U8.tar.gz -> version=2019_U8

tbb_version=`expr "$tbb_filename" : '^tbb-\(.*\)\.tar\.gz'`

## check git tree for modifications
if ! git diff-index --quiet HEAD --
then
    cat <<GIT_ERROR

    Please commit or stash any git changes prior to
    running this script. Quitting.

GIT_ERROR
    exit 1
fi

tbb_old_version=`expr "$(ls -d tbb_*/)" : '^tbb_\(.*\)/'`

echo "Old version:  $tbb_old_version"
echo "New version:  $tbb_version"

# 1. Remove the old version of TBB.
git rm -r tbb_${tbb_old_version}/
rm -rf tbb_${tbb_old_version}/
git commit -m "upgrading to TBB version ${tbb_version}; removing old TBB library"

# 2. Modify makefiles and README with new version number
# if the versions are the same we are only trimming files
if [ "$tbb_version" != "$tbb_old_version" ]; then

sed -i -e "s|lib/tbb_${tbb_old_version}|lib/tbb_${tbb_version}|g" ../README.md ../make/*
sed -i -e "s|TBB (version ${tbb_old_version})|TBB (version ${tbb_version})|g" ../README.md
rm -f ../README.md*-e ../make/*-e
git add -u ../README.md ../make/*
git commit -m "upgrading to TBB version ${tbb_version}; modifying with new version number"

fi
# 3. Unpack the new TBB version.
tar xvzf $tbb_filename
mv oneTBB-${tbb_version} tbb_${tbb_version}
git add tbb_${tbb_version}
git commit -m "upgrading to TBB version ${tbb_version}; adding unmodified TBB library"

# 4. Prune TBB
cd tbb_${tbb_version}

git rm -rf examples/*
git commit -m "upgrading to TBB version ${tbb_version}; pruning files"

cat <<EOF

    Please port the changes in the commit a38babf to the new version. This
    ensures that the -g flag is not set for release builds.

    Done upgrading TBB from version ${tbb_old_version} to version ${tbb_version}.

    Please check the upgrade worked by running a test that uses TBB
    Example (from Math home directory):
      ./runTests.py test/unit/math/prim/core/init_threadpool_tbb_test.cpp
EOF
