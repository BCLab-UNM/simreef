#!/bin/bash

arg=$1

if [ -z $arg ]
then
    arg=release
fi

cd simforager
echo ./build_hopper.sh $arg
./build_hopper.sh $arg
cd ..
