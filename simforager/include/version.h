
// version.h — Build Version Metadata
// ----------------------------------
// Contains extern C declarations for metadata describing the current build.
// Also defines compilation context macros for CUDA support (or lack thereof).

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// These strings are defined at build time and may be substituted by build system
extern const char *SIMFORAGER_VERSION;
extern const char *SIMFORAGER_VERSION_DATE;
extern const char *SIMFORAGER_BUILD_DATE;
extern const char *SIMFORAGER_BRANCH;

#define MAX_BUILD_KMER_STR "Maximum kmer len=MAX_BUILD_KMER"

// CUDA capability flag: allows conditional compilation
#ifdef USE_CUDA
#define WITH_CUDA "CUDA + CPU SW"
#else
#define WITH_CUDA "CPU SW Only"
#endif

#ifdef __cplusplus
}
#endif
