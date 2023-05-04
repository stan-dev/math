#!/bin/sh
#
# Copyright (c) 2023 alandefreitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# This code has been adapted from libs/json/fuzzing/fuzz.sh
# By Paul Dreik 2019-2020 for the boost json project
# License: Boost 1.0

set -e

# Current directory
fuzzdir=$(dirname $0)
me=$(basename $0)

cd $fuzzdir

# Find Clang
if [ -z $CLANG ]; then
  for clangver in -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -6.0 ""; do
    CLANG=clang++$clangver
    if which $CLANG >/dev/null; then
      break
    fi
  done
fi

if ! which $CLANG >/dev/null; then
  if ! -x $CLANG; then
    echo $me: sorry, could not find clang $CLANG
    exit 1
  fi
fi
echo "\$CLANG: $CLANG"

# Fuzzer parameters
MAXLEN="-max_len=4000" # max input size, to avoid blowing up the corpus size
JOBS= # Adjust this to utilize more of your cpu: JOBS="-jobs=32"
MAXTIME="-max_total_time=30" # seconds

# Compile and run each fuzzer
variants="parse"
for variant in $variants; do
  # Paths
  srcfile=fuzz_$variant.cpp
  fuzzer=./fuzzer_$variant
  seedcorpus=seedcorpus/$variant

  echo "Fuzz \"$variant\""
  echo "\$srcfile= $srcfile"
  echo "\$fuzzer= $fuzzer"
  echo "\$seedcorpus= $seedcorpus"

  # Compile
  if [ ! -e $fuzzer -o $srcfile -nt $fuzzer ]; then
    $CLANG \
      -std=c++11 \
      -O3 \
      -g \
      -fsanitize=fuzzer,address,undefined \
      -fno-sanitize-recover=undefined \
      -I../../include \
      -o $fuzzer \
      ../../src/src.cpp \
      $srcfile
  fi

  # Create initial corpus from old crashes in repo
  if [ -d old_crashes/$variant ]; then
    find old_crashes/$variant -type f -print0 | xargs -0 --no-run-if-empty $fuzzer
  fi
  if [ ! -d $seedcorpus ]; then
    mkdir -p $seedcorpus
    find ../test -name "*.json" -type f -print0 | xargs -0 --no-run-if-empty cp -f -t $seedcorpus/
  fi

  # Potentially merge file with old corpus
  if [ -e corpus.tar ]; then
    mkdir -p oldcorpus
    tar xf corpus.tar -C oldcorpus || echo "corpus.tar was broken! ignoring it"
    OLDCORPUS=oldcorpus/cmin/$variant
    mkdir -p $OLDCORPUS
  else
    OLDCORPUS=
  fi

  # Run fuzzer with corpus
  export UBSAN_OPTIONS="halt_on_error=1"
  mkdir -p out/$variant
  $fuzzer out/$variant $OLDCORPUS $seedcorpus/ $MAXTIME $MAXLEN $JOBS

  # Minimize the corpus
  mkdir -p cmin/$variant
  $fuzzer cmin/$variant $OLDCORPUS out/$variant $seedcorpus/ -merge=1 $MAXLEN
done
