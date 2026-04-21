#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <debug|release|clean>"
    exit 1
fi

case "$1" in
    debug)
        BUILD_TYPE=Debug
        export UPCXX_CODEMODE=debug
        ;;
    release)
        BUILD_TYPE=Release
        export UPCXX_CODEMODE=opt
        ;;
    clean)
        BUILD_TYPE=""
        ;;
    *)
        echo "Usage: $0 <debug|release|clean>"
        exit 1
        ;;
esac

export UPCXX_THREADMODE=seq
export UPCXX_NETWORK=ibv

module purge
module load boost/1.88.0-65dn
module load openmpi/4.1.7-e7k3
module load upcxx/2023.9.0-27k4
module load cmake/3.31.6-auvg
module load opencv/4.9.0-ap5r

if [ -n "${SIMREEF_BUILD_ENV:-}" ]; then
    source "$SIMREEF_BUILD_ENV"
fi

if [ -z "${OPENCV_LIB64:-}" ] && [ -z "${OPENCV_LIB:-}" ]; then
    echo "OpenCV module did not set OPENCV_LIB64 or OPENCV_LIB"
    exit 1
fi

export LD_LIBRARY_PATH="${OPENCV_LIB64:-}:${OPENCV_LIB:-}:${CARC_LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"

gcc_exec="$(which gcc || true)"
mpicxx_exec="$(which mpic++ || which mpicxx || true)"
upcxx_exec="$(which upcxx || true)"

if [ -z "$gcc_exec" ]; then
    echo "gcc not found after loading modules"
    exit 1
fi

if [ -z "$mpicxx_exec" ]; then
    echo "mpic++ not found after loading modules"
    exit 1
fi

if [ -z "$upcxx_exec" ]; then
    echo "upcxx not found after loading modules"
    exit 1
fi

# Keep wrapper compilers as wrappers. Do not resolve mpic++ to opal_wrapper.
gcc_exec="$(readlink -f "$gcc_exec")"

export CC="$gcc_exec"
export CXX="$mpicxx_exec"

rootdir="$(pwd)"
cluster="$(hostname | sed 's/[0-9]*$//')"
BUILD_PATH="$rootdir/.build_${cluster}"
INSTALL_PATH="${SIMREEF_INSTALL_PATH:-$rootdir/install_${cluster}}"

if [ "$1" = "clean" ]; then
    echo "Deleting $BUILD_PATH and $INSTALL_PATH"
    rm -rf "$BUILD_PATH" "$INSTALL_PATH" "$rootdir/src/.build"
    exit 0
fi

mkdir -p "$BUILD_PATH"
cd "$BUILD_PATH"

echo "Cluster:            $cluster"
echo "Build type:         $BUILD_TYPE"
echo "C compiler:         $gcc_exec"
echo "C++ compiler:       $mpicxx_exec"
echo "UPCXX wrapper:      $upcxx_exec"
echo "Build directory:    $BUILD_PATH"
echo "Install directory:  $INSTALL_PATH"
echo "UPCXX_CODEMODE:     $UPCXX_CODEMODE"
echo "UPCXX_THREADMODE:   $UPCXX_THREADMODE"
echo "UPCXX_NETWORK:      $UPCXX_NETWORK"

SECONDS=0
cmake "$rootdir" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_PATH" \
    -DCMAKE_C_COMPILER="$gcc_exec" \
    -DCMAKE_CXX_COMPILER="$mpicxx_exec"

echo "Configure took ${SECONDS}s"
echo "Installing to $INSTALL_PATH"

make -j"$(nproc)" install
