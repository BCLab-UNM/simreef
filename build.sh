#!/bin/bash

# Default values
cluster="hopper"
type="release"

# Parse named arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --cluster)
            cluster="$2"
            shift 2
            ;;
        --type)
            type="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--cluster <cluster>] [--type <type>]"
            exit 1
            ;;
    esac
done

cd simforager || exit 1
echo "./build_${cluster}.sh $type"
./build_"${cluster}".sh "$type"
cd ..
