#!/bin/bash
#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set -e
GLOBIGNORE=".:.."

# Install dependencies. rsync is needed by the GitHub action to upload the pages
apt-get update
apt-get install --no-install-recommends -y \
    doxygen \
    xsltproc \
    wget \
    ca-certificates \
    clang-11 \
    python3 \
    python-is-python3 \
    rsync
ln -s /usr/bin/clang++-11 /usr/bin/clang++
ln -s /usr/bin/clang-11 /usr/bin/clang

# Install saxonhe
mkdir -p /tmp/saxonhe
cd /tmp/saxonhe
wget -q -O saxonhe.zip https://sourceforge.net/projects/saxon/files/Saxon-HE/9.9/SaxonHE9-9-1-4J.zip/download
unzip -o -qq saxonhe.zip
mkdir -p /usr/share/java/
mv saxon9he.jar /usr/share/java/Saxon-HE.jar
cd
rm -rf /tmp/saxonhe

# Install docbook XSL stylesheets
mkdir -p $DOCBOOK_XSL_DIR
cd $DOCBOOK_XSL_DIR
wget -q https://github.com/docbook/xslt10-stylesheets/releases/download/snapshot%2F2020-06-03/docbook-xsl-snapshot.zip
unzip -o -qq docbook-xsl-snapshot.zip
mv docbook-xsl-snapshot/* .
rmdir docbook-xsl-snapshot
rm docbook-xsl-snapshot.zip

# Install docbook DTD
mkdir -p $DOCBOOK_DTD_DIR
cd $DOCBOOK_DTD_DIR
wget -q http://www.oasis-open.org/docbook/xml/4.2/docbook-xml-4.2.zip
unzip -o -qq docbook-xml-4.2.zip
rm docbook-xml-4.2.zip
