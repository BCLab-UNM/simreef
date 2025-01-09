#!/usr/bin/env bash

# Check if an argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <debug|release|clean>"
    exit 1
fi

# Validate the argument
if [ "$1" != "debug" ] && [ "$1" != "release" ] && [ "$1" != "clean" ]; then
    echo "Invalid argument: '$1'"
    echo "Usage: $0 <debug|release|clean>"
    exit 1
fi

export UPCXX_THREADMODE=seq
export UPCXX_CODEMODE=opt
export UPCXX_NETWORK=ibv
module purge
module load gcc/12.1.0-crtl
module load cmake/3.11.4-qkyj
module load openmpi/4.1.3-j6zb
module load upcxx/2020.10.0-6eh2

if [ -n "$SIMREEF_BUILD_ENV" ]; then
    source $SIMREEF_BUILD_ENV
fi

upcxx_exec=`which upcxx`

if [ -z "$upcxx_exec" ]; then
    echo "upcxx not found. Please load the appropriate module."
    exit 1
fi

upcxx_exec_canonical=$(readlink -f $upcxx_exec)
if [ "$upcxx_exec_canonical" != "$upcxx_exec" ]; then
    echo "Found symlink for upcxx - using target at $upcxx_exec_canonical"
    export PATH=`dirname $upcxx_exec_canonical`:$PATH
fi

set -e

rootdir=`pwd`

INSTALL_PATH=${SIMREEF_INSTALL_PATH:=$rootdir/install}

rm -rf $INSTALL_PATH/bin/simreef

BUILD_PATH=build

if [ "$1" == "clean" ]; then
    echo Deleting $BUILD_PATH and $INSTALL_PATH
    rm -rf $BUILD_PATH
    rm -rf $INSTALL_PATH
    exit 0
else
    mkdir -p $rootdir/.build
    cd $rootdir/.build
    if [ "$1" == "debug" ] || [ "$1" == "release" ]; then # Seems like this will always be true
	
        #rm -rf *
        #rm -rf $INSTALL_PATH/cmake
	SECONDS=0
    	cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH -DCMAKE_CXX_COMPILER=/opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/openmpi-4.1.3-j6zbgs4rx7w7mb4imwl6fqk2wxvglehb/bin/mpicxx
	echo "Build took $((SECONDS))s"
    fi
    echo "Installing to $INSTALL_PATH"
    make -j install
fi



