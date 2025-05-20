#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

FROM ubuntu:24.04

RUN \
    apt-get update && \
    apt-get --no-install-recommends -y install \
        gcc-14 \
        g++-14 \
        libssl-dev \
        git \
        ca-certificates \
        python3 \
        python3-requests \
        python-is-python3 \
        cmake \
        ninja-build \
        valgrind \
        gpg \
        gpg-agent \
        mysql-client && \
    ln -s /usr/bin/g++-14 /usr/bin/g++ && \
    ln -s /usr/bin/gcc-14 /usr/bin/gcc
