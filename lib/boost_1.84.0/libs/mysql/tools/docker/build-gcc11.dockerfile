#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

FROM ubuntu:22.04

# DEBIAN_FRONTEND is needed because lcov pulls the package
# tzdata, which by default attempts to read interactively from stdin
RUN \
    export DEBIAN_FRONTEND=noninteractive && \
    apt-get update && \
    apt-get --no-install-recommends -y install \
        ca-certificates \
        gcc-11 \
        g++-11 \
        lcov \
        libssl-dev \
        wget \
        python3 \
        python-is-python3 \
        ninja-build \
        valgrind \
        git \
        gpg \
        cmake \
        gpg-agent \
        mysql-client && \
    unlink /usr/bin/gcc && \
    ln -s /usr/bin/g++-11 /usr/bin/g++ && \
    ln -s /usr/bin/gcc-11 /usr/bin/gcc && \
    wget https://keybase.io/codecovsecurity/pgp_keys.asc && \
    gpg --no-default-keyring --keyring trustedkeys.gpg --import pgp_keys.asc && \
    wget -q https://uploader.codecov.io/latest/linux/codecov && \
    wget -q https://uploader.codecov.io/latest/linux/codecov.SHA256SUM && \
    wget -q https://uploader.codecov.io/latest/linux/codecov.SHA256SUM.sig && \
    gpgv codecov.SHA256SUM.sig codecov.SHA256SUM && \
    shasum -a 256 -c codecov.SHA256SUM && \
    chmod +x codecov && \
    mv codecov /usr/bin/ && \
    rm codecov*
