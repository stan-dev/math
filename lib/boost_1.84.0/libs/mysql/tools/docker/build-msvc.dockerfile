#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

ARG BASE_IMAGE=cppalliance/dronevs2022:1
FROM ${BASE_IMAGE}

COPY tools/win-ci.cnf C:/my.cnf
COPY tools/ssl C:/ssl

RUN choco install --no-progress -y mysql --version 8.0.20

RUN powershell -Command \
    $ErrorActionPreference = 'Stop' ; \
    $OPENSSL_VERSION = '1.1.1.1900' ; \
    choco install --no-progress -y --force --forceX86 --version $OPENSSL_VERSION openssl ; \
    Copy-Item 'C:/Program Files (x86)/OpenSSL-Win32/' -Recurse -Destination 'C:/openssl-32/' ; \
    choco install --no-progress -y --force --version $OPENSSL_VERSION openssl ; \
    Copy-Item 'C:/Program Files/OpenSSL-Win64/' -Recurse -Destination 'C:/openssl-64/' ; \
    choco uninstall -y openssl
