#!/bin/sh
# Copyright (C) 2022 John Maddock
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

# boostinspect:notab - Tabs are required for the Makefile.

# boost-root:
cd ../../..
./b2 -j3 tools/bcp//dist  || { echo './b2 -j3 tools/bcp failed' ; exit 1; }
rm -rf ../bcp-output
mkdir ../bcp-output  || { echo 'mkdir failed' ; exit 1; }
./dist/bin/bcp program_options build ../bcp-output  || { echo './dist/bin/bcp program_options build ../bcp-output failed' ; exit 1; }

# tests:
cd ../bcp-output
ls -l
./bootstrap.sh  || { echo './bootstrap.sh failed' ; exit 1; }
./b2 -j3 libs/program_options/build  || { echo './b2 -j3 libs/program_options/build failed' ; exit 1; }
./b2 -j3 --with-program_options || { echo './b2 -j3 --with-program_options failed' ; exit 1; }
./b2 -j3  || { echo './b2 -j3 failed' ; exit 1; }

# run inspect:
cd ../boost-root
echo "using xsltproc ;" | tee -a $HOME/user-config.jam
echo "using boostbook : /usr/share/xml/docbook/stylesheet/docbook-xsl : /usr/share/sgml/docbook/dtd/xml/4.2 ;" | tee -a $HOME/user-config.jam
cd tools/inspect/build
../../../b2 -j2 address-model=64 architecture=x86 toolset=gcc cxxstd=14 release dist-bin || { echo 'building inspect failed' ; exit 1; }
cd ../../bcp/doc
rm -rf html
../../../b2 -j2 address-model=64 architecture=x86 toolset=gcc cxxstd=14 release || { echo 'building docs failed' ; exit 1; }
../../../dist/bin/inspect -text || { echo 'inspect failed' ; exit 1; }

