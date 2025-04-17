
// utils.hpp — Utility Function Declarations
// -----------------------------------------
// Provides shared functions and constants used throughout the simulation,
// including:
// - BMP image parsing
// - Debugging image formats
// - System utilities (e.g., thread pinning)
// - Parallel file output helpers

#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>
#include <stdint.h>
#include <sys/types.h>

#include <upcxx/upcxx.hpp>
#include <upcxx/atomic_domain.hpp>

using std::string;
using std::vector;

// Forward declarations
class IntermittentTimer;

// Reads a 24-bit BMP file and returns a 2D grid of values based on RGB pixel color
vector<vector<uint8_t>> readBMPColorMap(const string& file_name);

// Re-validates the decoded map against the BMP file it was read from
void debugColorMapData(const string& file_name, const vector<vector<uint8_t>>& color_map);

// Pins a given process/thread ID to a specific core index
int pin_thread(pid_t pid, int cid);

// Writes string content to a shared file across ranks using atomic offset updates
void dump_single_file(const string &fname, const string &out_str);

#endif
