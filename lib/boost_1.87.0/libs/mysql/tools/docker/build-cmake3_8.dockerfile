#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

FROM ubuntu:16.04

RUN \
    apt-get update && \
    apt-get --no-install-recommends -y install \
        software-properties-common \
        ca-certificates \
        git \
        python3 \
        python3-requests \
        g++ \
        gnupg-curl \
        libssl-dev \
        ninja-build \
        curl \
        mysql-client && \
    ln -s /usr/bin/python3 /usr/bin/python

RUN \
    mkdir -p ~/cmake && \
    CMAKE_VERSION=3.8.0 && \
    curl -OLs https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-SHA-256.txt && \
    curl -OLs https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    sha256sum -c --ignore-missing cmake-${CMAKE_VERSION}-SHA-256.txt && \
    curl -OLs https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-SHA-256.txt.asc && \
    gpg --keyserver hkps://keyserver.ubuntu.com --recv-keys C6C265324BBEBDC350B513D02D2CEF1034921684 && \
    gpg --verify cmake-${CMAKE_VERSION}-SHA-256.txt.asc cmake-${CMAKE_VERSION}-SHA-256.txt && \
    bash cmake-${CMAKE_VERSION}-Linux-x86_64.sh --skip-license && \
    cd ~ && \
    rm -rf ~/cmake && \
    cmake --version
