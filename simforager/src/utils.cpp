
// utils.cpp — Simulation Utility Functions
// ----------------------------------------
// Contains helper routines for:
// - Loading BMP maps of the simulation grid (e.g., ecosystem structure)
// - Debugging and validating BMP color mappings
// - Thread pinning and parallel-safe output to a single file using UPC++ atomics


#include <fstream>
#include <string_view>
#include "upcxx_utils.hpp"
using std::ifstream;
using std::string_view;
using namespace upcxx_utils;

#include "utils.hpp"

using namespace upcxx_utils;

using std::min;
using std::string;
using std::string_view;
using std::to_string;

// Reads a 24-bit BMP file and converts its RGB values into a 2D map of
// integer codes, corresponding to simulation substrate types or terrain.
// Expected values:
//   (255, 0, 0)   → 1 (Red)
//   (0, 255, 0)   → 2 (Green)
//   (0, 0, 255)   → 3 (Blue)
//   (0, 0, 0)     → 0 (Black)
// Any deviation from these exact pixel values will raise an error.
vector<vector<uint8_t>> readBMPColorMap(const string& file_name) {
    ifstream file(file_name, ios::binary);
    if (!file) {
        DIE("Failed to open file ", file_name, "\n");
        return {};
    }

    // Standard BMP header is 54 bytes
    uint8_t header[54];
    file.read(reinterpret_cast<char*>(header), 54);

    // Extract metadata from the header
    int width = *(int*)&header[18];
    int height = *(int*)&header[22];
    int data_offset = *(int*)&header[10];
    int row_padded = (width * 3 + 3) & (~3); // pad rows to 4-byte alignment

    vector<vector<uint8_t>> color_map(height, vector<uint8_t>(width, 0));
    file.seekg(data_offset);  // move to pixel data

    // BMP stores rows from bottom to top
    for (int i = height - 1; i >= 0; --i) {
        vector<uint8_t> row(row_padded);
        file.read(reinterpret_cast<char*>(row.data()), row_padded);

        for (int j = 0; j < width; ++j) {
            uint8_t blue  = row[j * 3 + 0];
            uint8_t green = row[j * 3 + 1];
            uint8_t red   = row[j * 3 + 2];

            if (red == 255 && green == 0 && blue == 0) {
                color_map[i][j] = 1;
            } else if (green == 255 && red == 0 && blue == 0) {
                color_map[i][j] = 2;
            } else if (blue == 255 && red == 0 && green == 0) {
                color_map[i][j] = 3;
            } else if (red == 0 && green == 0 && blue == 0) {
                color_map[i][j] = 0;
            } else {
                DIE("Unexpected pixel value at (", i, ", ", j, "): ",
                    "R=", red, " G=", green, " B=", blue, "\n");
            }
        }
    }

    return color_map;
}

// Debug utility to verify that the RGB → color_map encoding is correct
void debugColorMapData(const string& file_name, const vector<vector<uint8_t>>& color_map) {
    ifstream file(file_name, ios::binary);
    if (!file) {
        DIE("Failed to reopen file ", file_name, "\n");
        return;
    }

    uint8_t header[54];
    file.read(reinterpret_cast<char*>(header), 54);

    int width = *(int*)&header[18];
    int height = *(int*)&header[22];
    int data_offset = *(int*)&header[10];
    int row_padded = (width * 3 + 3) & (~3);

    if (color_map.size() != height || color_map[0].size() != width) {
        DIE("ERROR: Image dimensions don't match BMP header.\n");
        return;
    }

    file.seekg(data_offset);
    int count_red = 0, count_green = 0, count_blue = 0, count_black = 0;
    int mismatches = 0;
    uint64_t total_original_pixel_value = 0;
    uint64_t total_mapped_value_sum = 0;

    for (int i = height - 1; i >= 0; --i) {
        vector<uint8_t> row(row_padded);
        file.read(reinterpret_cast<char*>(row.data()), row_padded);

        for (int j = 0; j < width; ++j) {
            uint8_t blue  = row[j * 3 + 0];
            uint8_t green = row[j * 3 + 1];
            uint8_t red   = row[j * 3 + 2];
            total_original_pixel_value += red + green + blue;

            uint8_t expected_value = 0;
            if (red == 255 && green == 0 && blue == 0) {
                expected_value = 1; count_red++;
            } else if (green == 255 && red == 0 && blue == 0) {
                expected_value = 2; count_green++;
            } else if (blue == 255 && red == 0 && green == 0) {
                expected_value = 3; count_blue++;
            } else if (red == 0 && green == 0 && blue == 0) {
                expected_value = 0; count_black++;
            } else {
                DIE("Unexpected pixel value at (", i, ", ", j, ")\n");
            }

            total_mapped_value_sum += color_map[i][j];

            if (color_map[i][j] != expected_value) {
                mismatches++;
                if (mismatches <= 10) {
                    DIE("Mismatch at (", i, ", ", j, "): ",
                        "Expected ", (int)expected_value,
                        " Got ", (int)color_map[i][j], "\n");
                }
            }
        }
    }

    SLOG("BMP Color Map Statistics:\n",
         "  Red Pixels (1):   ", count_red, "\n",
         "  Green Pixels (2): ", count_green, "\n",
         "  Blue Pixels (3):  ", count_blue, "\n",
         "  Black Pixels (0): ", count_black, "\n",
         "  Mismatches: ", mismatches, "\n");

    if (mismatches == 0) {
        SLOG("✅ Color map matches original BMP pixel data.\n");
    } else {
        DIE("❌ Found ", mismatches, " mismatches.\n");
    }
}

// Binds a Linux thread to a specific CPU core
int pin_thread(pid_t pid, int cid) {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  CPU_SET(cid, &cpu_set);
  if (sched_setaffinity(pid, sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), pid);
    return -1;
  }
  return 0;
}

// Uses UPC++ atomics to allow parallel ranks to write disjoint file sections
void dump_single_file(const string &fname, const string &out_str) {
  auto sz = out_str.length();
  upcxx::atomic_domain<size_t> ad({upcxx::atomic_op::fetch_add, upcxx::atomic_op::load});
  upcxx::global_ptr<size_t> fpos = nullptr;
  if (!upcxx::rank_me()) fpos = upcxx::new_<size_t>(sz);
  fpos = upcxx::broadcast(fpos, 0).wait();
  size_t my_fpos = 0;
  if (upcxx::rank_me()) my_fpos = ad.fetch_add(fpos, sz, std::memory_order_relaxed).wait();
  upcxx::barrier();
  int fileno = -1;
  size_t fsize = 0;
  if (!upcxx::rank_me()) {
    fsize = ad.load(fpos, std::memory_order_relaxed).wait();
    fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) WARN("Error creating file ", fname, ": ", strerror(errno), "\n");
    if (ftruncate(fileno, fsize) == -1)
      WARN("Could not truncate ", fname, " to ", fsize, " bytes\n");
    close(fileno);
  }
  upcxx::barrier();
  fileno = open(fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) DIE("Cannot open file ", fname, ": ", strerror(errno), "\n");
  auto bytes_written = pwrite(fileno, out_str.data(), sz, my_fpos);
  if (bytes_written != sz)
    DIE("Could not write all bytes to ", fname, "; only wrote ", bytes_written);
  close(fileno);
  ad.destroy();
  if (!upcxx::rank_me()) upcxx::delete_(fpos);
}
