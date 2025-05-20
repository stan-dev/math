#!/bin/bash
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

benchs=(
   "nopool-tcp"
   "nopool-tcpssl"
   "nopool-unix"
   "pool-tcp"
   "pool-tcpssl"
   "pool-unix"
)

outfile=private/benchmark-results.txt

echo "bench,ellapsed" > $outfile

for bench in ${benchs[@]}
do
   echo $bench
   for i in {1..10}
   do
      ellapsed=$(./__build/bench/boost_mysql_bench_connection_pool $bench localhost)
      echo "$bench,$ellapsed" | tee -a $outfile
   done
done