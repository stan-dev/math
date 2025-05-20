#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

FROM ubuntu:23.04

RUN \
    apt-get update && \
    apt-get --no-install-recommends -y install \
        ca-certificates \
        clang-16 \
        llvm-16 \
        libclang-rt-16-dev \
        libc++-16-dev \
        libc++abi-16-dev \
        libssl-dev \
        git \
        ninja-build \
        cmake \
        python3 \
        python3-requests \
        python-is-python3 \
        mysql-client && \
    ln -s /usr/bin/clang++-16 /usr/bin/clang++ && \
    ln -s /usr/bin/clang-16 /usr/bin/clang

