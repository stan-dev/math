#!/bin/bash
#===============================================================================
#
# This script is for the Math library maintainers to use upgrade Boost.
# To update to a new Boost version, download the boost_<version>.tar.gz to
# the lib/ directory.
#   https://www.boost.org/users/download/
# This script will also manage the git process so the different steps are in
# separate git commits.
#
# This script will take these steps and add a git commit after each:
# 1. Remove the old version of Boost.
# 2. Modify makefiles and README with new version numbers
# 3. Unpack the new Boost version.
# 4. Prune Boost
#
# This script should be run from the lib/ folder and the argument should be the
# name of the boost_*.tar.gz file.
#===============================================================================
set -e
err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

usage() {
    cat <<HELP_USAGE

    Usage:
        $0 boost_<version>.tar.gz

    Download the new Boost version to boost_<version>.tar.gz in this
    directory. Provide that filename as the only parameter to this script.

    This script will upgrade the Boost library to the one provide,
    making modifications to the library as necessary. The modifications
    will be separate git commits for verification.

    If there are any unstaged modifications or staged modifications that
    haven't been git committed, this script will quit prior to adding any
    git commits.

HELP_USAGE
    exit 0
}

[ -z "$1" ] && { usage; }


boost_filename=$1
## extract version number from argument
## ex: boost_1_69_0.tar.gz -> version=1.69.0

boost_version=`expr "$boost_filename" : '^boost_\(.*\)\.tar\.gz'`
boost_version_underscore=$boost_version
#replace underscores with dots
boost_version=${boost_version//_/.}
if [[ -z "$boost_version" ]] ; then
    cat <<BOOST_VERSION_ERROR

    Can not parse the Boost version from the filename.
    Expecting: boost_<version>.tar.gz

    Is the filename correct? "$boost_filename"

BOOST_VERSION_ERROR
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

boost_old_version=`expr "$(ls -d boost_*/)" : '^boost_\(.*\)/'`
boost_old_version=${boost_old_version//_/.}

echo "Old version:  $boost_old_version"
echo "New version:  $boost_version"

# 1. Remove the old version of Boost.
git rm -r boost_${boost_old_version}/
rm -rf boost_${boost_old_version}/
git commit -m "upgrading to boost v${boost_version}; removing old boost library"

# 2. Modify makefiles and README with new version number
# if the versions are the same we are only trimming files
if [ "$boost_version" != "$boost_old_version" ]; then

sed -i -e "s|lib/boost_${boost_old_version}|lib/boost_${boost_version}|g" ../README.md ../make/*
sed -i -e "s|Boost (version ${boost_old_version})|Boost (version ${boost_version})|g" ../README.md
rm -f ../README.md*-e ../make/*-e
git add -u ../README.md ../make/*
git commit -m "upgrading to boost v${boost_version}; modifying with new version number"

fi
# 3. Unpack the new Boost version.
tar xvzf $boost_filename
mv boost_${boost_version_underscore} boost_${boost_version}
git add boost_${boost_version}
git commit -m "upgrading to boost v${boost_version}; adding unmodified boost library"

# 4. Prune Boost
cd boost_${boost_version}

git rm -rf boost.css boost.png index.htm index.html INSTALL rst.css --ignore-unmatch
git rm -rf doc/ libs/*/test/ libs/*/example/ libs/*/doc/ libs/*/*/doc/ --ignore-unmatch
git rm -rf boost/leaf/ libs/leaf/ boost/nowide/ libs/nowide --ignore-unmatch
git rm -rf boost/pfr libs/pfr boost/json/ libs/json/ --ignore-unmatch
git rm -rf boost/static_string libs/static_string --ignore-unmatch
git rm -rf boost/stl_interfaces/ libs/stl_interfaces --ignore-unmatch
git rm -rf boost/phoenix/ libs/phoenix --ignore-unmatch
git rm -rf boost/msm/ libs/msm --ignore-unmatch
git rm -rf boost/redis/ libs/redis --ignore-unmatch
git rm -rf boost/mysql/ libs/mysql --ignore-unmatch
git rm -rf boost/log/ libs/log --ignore-unmatch
git rm -rf libs/*/examples libs/*/*/examples --ignore-unmatch
git rm -rf **/*.svg **/*.png **/*.jpg **/*.html **/*.htm **/*.gold **/*.json --ignore-unmatch
git rm -rf **/*.pdf **/*.manifest **/*.css **/*.md **/*.qbk **/*.rst **/*.txt --ignore-unmatch
git commit -m "upgrading to boost v${boost_version}; pruning files"

cat <<EOF

    Done upgrading Boost from v${boost_old_version} to v${boost_version}.

    Please check the upgrade worked by running a test that uses Boost
    Example (from Math home directory):
      ./runTests.py test/unit/math/prim/prob/multi_normal_test.cpp

EOF
