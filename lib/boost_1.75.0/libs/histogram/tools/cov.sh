#!/bin/sh
# must be executed in project root folder

# Copyright Hans Dembinski 2018-2019
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

if [ -z $GCOV ]; then
  # gcov-9, gcov-7, gcov-6 do not work
  for i in 8 5; do
    if test $(which gcov-$i); then
      GCOV=gcov-$i
      break;
    fi;
  done
fi

LCOV_VERSION="1.14"
LCOV_DIR="tools/lcov-${LCOV_VERSION}"

if [ ! -e $LCOV_DIR ]; then
  cd tools
  curl -L https://github.com/linux-test-project/lcov/releases/download/v${LCOV_VERSION}/lcov-${LCOV_VERSION}.tar.gz | tar zxf -
  cd ..
fi

# --rc lcov_branch_coverage=1 doesn't work on travis
# LCOV="${LCOV_DIR}/bin/lcov --gcov-tool=${GCOV} --rc lcov_branch_coverage=1"
LCOV="${LCOV_DIR}/bin/lcov --gcov-tool=${GCOV}"

# collect raw data
$LCOV --base-directory `pwd` \
  --directory `pwd`/../../bin.v2/libs/histogram/test \
  --capture --output-file coverage.info

# remove uninteresting entries
$LCOV --extract coverage.info "*/boost/histogram/*" --output-file coverage.info

if [ $1 ]; then
  # upload if on CI or when token is passed as argument
  which cpp-coveralls || echo "Error: you need to install cpp-coveralls"
  cpp-coveralls -l coverage.info -r ../.. -n -t $1
elif [ ! $CI ]; then
  # otherwise generate html report
  $LCOV_DIR/bin/genhtml coverage.info --demangle-cpp -o coverage-report
fi
