#!/bin/bash

set -xe

export REPONAME=$(basename -s .git `git config --get remote.origin.url`)
export BOOST_SRC_FOLDER=$(git rev-parse --show-toplevel)

# Determine if running in boost-root
PARENTNAME=$(basename -s .git `git --git-dir ${BOOST_SRC_FOLDER}/../../.git config --get remote.origin.url`) || true
if [ -n "${PARENTNAME}" -a "${PARENTNAME}" = "boost" ]; then
    echo "Starting out inside boost-root"
    BOOSTROOTLIBRARY="yes"
else
    echo "Not starting out inside boost-root"
    BOOSTROOTLIBRARY="no"
fi

if git rev-parse --abbrev-ref HEAD | grep master ; then BOOST_BRANCH=master ; else BOOST_BRANCH=develop ; fi

echo '==================================> INSTALL'

brew install doxygen
brew install wget
brew tap adoptopenjdk/openjdk
brew install --cask adoptopenjdk11
brew install gnu-sed
brew install docbook
brew install docbook-xsl

# an alternate set of packages from https://www.boost.org/doc/libs/develop/doc/html/quickbook/install.html
# sudo port install libxslt docbook-xsl docbook-xml-4.2

export PATH="/usr/local/opt/gnu-sed/libexec/gnubin:$PATH"

cd $BOOST_SRC_FOLDER
cd ..
mkdir -p tmp && cd tmp
if [ ! -f saxonhe.zip ]; then wget -O saxonhe.zip https://sourceforge.net/projects/saxon/files/Saxon-HE/9.9/SaxonHE9-9-1-4J.zip/download; fi
unzip -d saxonhe -o saxonhe.zip
cd saxonhe

sudo rm /Library/Java/Extensions/Saxon-HE.jar || true
sudo cp saxon9he.jar /Library/Java/Extensions/Saxon-HE.jar
sudo chmod 755 /Library/Java/Extensions/Saxon-HE.jar

cd $BOOST_SRC_FOLDER

if [ "${BOOSTROOTLIBRARY}" = "yes" ]; then
    cd ../..
    git checkout $BOOST_BRANCH
    git pull
    export BOOST_ROOT=$(pwd)
else
    cd ..
    if [ ! -d boost-root ]; then
        git clone -b $BOOST_BRANCH https://github.com/boostorg/boost.git boost-root --depth 1
        cd boost-root
        export BOOST_ROOT=$(pwd)
        rsync -av $BOOST_SRC_FOLDER/ libs/$REPONAME
    else
        cd boost-root
        git checkout $BOOST_BRANCH
        git pull
        export BOOST_ROOT=$(pwd)
        rsync -av $BOOST_SRC_FOLDER/ libs/$REPONAME
    fi
fi

git submodule update --init libs/context
git submodule update --init tools/boostbook
git submodule update --init tools/boostdep
git submodule update --init tools/docca
git submodule update --init tools/quickbook
git submodule update --init tools/build
sed -i 's~GLOB "/usr/share/java/saxon/"~GLOB "/Library/Java/Extensions/" "/usr/share/java/saxon/"~' tools/build/src/tools/saxonhe.jam

python3 tools/boostdep/depinst/depinst.py ../tools/quickbook
./bootstrap.sh
./b2 headers

echo '==================================> COMPILE'

echo "using doxygen ; using boostbook ; using saxonhe ;" > tools/build/src/user-config.jam

./b2 libs/$REPONAME/doc/
