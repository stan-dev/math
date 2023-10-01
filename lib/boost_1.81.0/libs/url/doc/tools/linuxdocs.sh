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

sudo apt-get update
sudo apt-get install -y docbook docbook-xml docbook-xsl xsltproc libsaxonhe-java default-jre-headless flex libfl-dev bison unzip rsync wget python3 cmake build-essential

cd $BOOST_SRC_FOLDER
cd ..
mkdir -p tmp && cd tmp

if which doxygen; then
    echo "doxygen found"
else
    echo "building doxygen"
    if [ ! -d doxygen ]; then git clone -b 'Release_1_8_15' --depth 1 https://github.com/doxygen/doxygen.git && echo "not-cached" ; else echo "cached" ; fi
    cd doxygen
    cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
    cd build
    sudo make install
    cd ../..
fi

if [ ! -f saxonhe.zip ]; then wget -O saxonhe.zip https://sourceforge.net/projects/saxon/files/Saxon-HE/9.9/SaxonHE9-9-1-4J.zip/download && echo "not-cached" ; else echo "cached" ; fi
unzip -d saxonhe -o saxonhe.zip
cd saxonhe
sudo rm /usr/share/java/Saxon-HE.jar || true
sudo cp saxon9he.jar /usr/share/java/Saxon-HE.jar

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
python3 tools/boostdep/depinst/depinst.py ../tools/quickbook
./bootstrap.sh
./b2 headers

echo '==================================> COMPILE'

echo "using doxygen ; using boostbook ; using saxonhe ;" > tools/build/src/user-config.jam

./b2 libs/$REPONAME/doc/
