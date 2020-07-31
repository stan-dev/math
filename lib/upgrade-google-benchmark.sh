#!/bin/bash
#===============================================================================
#
# This script is for the Math library maintainers to use upgrade Google Benchmark
# which includes Google Test.
#
# To update to a new Google Benchmark version, provide the version number as an
# argument to the script.
#
# This script will also manage the git process so the different steps are in
# separate git commits.
#
# This script will take these steps and add a git commit after each:
# 1. Remove the old version of Google Benchmark.
# 2. Modify makefiles and README with new version numbers
# 3. Commit the new version of Google Benchmark
# 4. Commit the new version of Google Test
# 5. Modify google test's gtest_main.cc for MPI
#
# This script should be run from the lib/ folder and the argument should be the
# version of Google Benchmark: https://github.com/google/benchmark
#===============================================================================
set -e
err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

usage() {
    cat <<HELP_USAGE

    Usage:
        $0 <Google benchmark version>

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

## check git tree for modifications
if ! git diff-index --quiet HEAD --
then
    cat <<GIT_ERROR

    Please commit or stash any git changes prior to
    running this script. Quitting.

GIT_ERROR
    exit 1
fi

version=$1
old_version=`expr "$(ls -d benchmark_*/ 2>/dev/null)" : '^benchmark_\(.*\)/' || echo ''`

echo "Old version: $old_version"
echo "New version: $version"


# 1. Remove old version of benchmark

if [ ! -z "$old_version" ]
then
    git rm -r benchmark_${old_version}/
    rm -rf benchmark_${old_version}/
    git commit -m "upgrading to google benchmark v${version}; removing old library v${old_version}"
fi


# 2. Modify makefiles with new version number

sed -i -e "s|lib/benchmark_${old_version}|lib/benchmark_${version}|g" ../make/*
rm -f ../make/*-e
git add -u ../make/*
git commit -m "upgrading to google benchmark v${version}; modifying with new version number" || echo 'same google benchmark version'


# 3. Commit the new version of Google Benchmark
git clone --depth=1 --branch=v${version} https://github.com/google/benchmark.git benchmark_${version}
rm -rf benchmark_${version}/.git
git add benchmark_${version}
cp benchmark_${version}/LICENSE ../licenses/googlebenchmark-license.txt
git add ../licenses/googlebenchmark-license.txt
git commit -m "upgrading to google benchmark v${version}; adding google benchmark" || echo 'same google benchmark version'


# 4. Commit the new version of Google Test
git clone --depth=1 --branch=master https://github.com/google/googletest.git benchmark_${version}/googletest
gtest_version=$(cd benchmark_${version}/googletest && git log --pretty=tformat:"%h" -n1 .)
rm -rf benchmark_${version}/googletest/.git
git add -f benchmark_${version}/googletest
cp benchmark_${version}/googletest/LICENSE ../licenses/googletest-license.txt
git add ../licenses/googletest-license.txt
git commit -m "upgrading to google benchmark v${version}; adding google test version ${gtest_version}"

# 5. Modify google test's gtest_main.cc for MPI
#    Add two things:
#    1. include stan/math/prim/functor/mpi_cluster.hpp
#    2. start mpi cluster.

main=benchmark_${version}/googletest/googletest/src/gtest_main.cc

lineno=$(grep -n "#include" ${main} | tail -1)
lineno=${lineno%%:*}

sed -i -e "${lineno}a\\
\#include <stan/math/prim/functor/mpi_cluster.hpp>\\
" ${main}


lineno=$(grep -n "int main" ${main} | tail -1)
lineno=${lineno%%:*}

sed -i -e "${lineno}a\\
\\
#ifdef STAN_MPI\\
\ \ // for MPI testing we test with all workers in listen mode.\\
\ \ // No output is generated from the workers.\\
\ \ stan::math::mpi_cluster cluster;\\
\ \ cluster.listen();\\
\ \ if (cluster.rank_ != 0)\\
\ \ \ \ return 0;\\
#endif\\
\\
" ${main}

git add ${main}
git commit -m "upgrading to google benchmark v${version}; modifying ${main} for MPI"


cat <<EOF

    Done upgrading Google Benchmark from v${old_version} to v${version}.
    Google Test was also upgraded to ${gtest_version}.

    Please check the upgrade worked by running any test.
    Example (from Math home directory):
      ./runTests.py test/unit/math/prim/prob/multi_normal_test.cpp

EOF
