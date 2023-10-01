#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

FROM ubuntu:20.04

RUN \
    apt-get update && \
    apt-get --no-install-recommends -y install \
        ca-certificates \
        clang-11 \
        libssl-dev \
        git \
        python3 \
        python-is-python3 \
        mysql-client && \
    ln -s /usr/bin/clang++-11 /usr/bin/clang++ && \
    ln -s /usr/bin/clang-11 /usr/bin/clang

