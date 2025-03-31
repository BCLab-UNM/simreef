#include "utils.hpp"

using namespace upcxx_utils;

using std::min;
using std::string;
using std::string_view;
using std::to_string;

// Function to read a 24-bit BMP file into a 2D vector where:
//  0 = Black pixel (0, 0, 0)
//  1 = Red-only pixel (255, 0, 0)
//  2 = Green-only pixel (0, 255, 0)
//  3 = Blue-only pixel (0, 0, 255)
// Assumes all pixel values are exactly 0 or 255.
vector<vector<uint8_t>> readBMPColorMap(const string& file_name) {
    ifstream file(file_name, ios::binary);
    if (!file) {
        DIE("Failed to open file ", file_name, "\n");
        return {};
    }

    // Read standard BMP header (54 bytes)
    uint8_t header[54];
    file.read(reinterpret_cast<char*>(header), 54);

    // Extract image metadata
    int width = *(int*)&header[18];
    int height = *(int*)&header[22];
    int data_offset = *(int*)&header[10];
    int row_padded = (width * 3 + 3) & (~3);  // Each row is padded to a 4-byte boundary

    // Prepare the output 2D vector to hold the encoded pixel values
    vector<vector<uint8_t>> color_map(height, vector<uint8_t>(width, 0));

    // Seek to the beginning of pixel data
    file.seekg(data_offset);

    // BMP stores rows bottom-up, so iterate in reverse row order
    for (int i = height - 1; i >= 0; --i) {
        vector<uint8_t> row(row_padded);
        file.read(reinterpret_cast<char*>(row.data()), row_padded);

        for (int j = 0; j < width; ++j) {
            uint8_t blue  = row[j * 3 + 0];
            uint8_t green = row[j * 3 + 1];
            uint8_t red   = row[j * 3 + 2];

            // Map RGB values to encoded colour class
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

// Function to validate that a BMP was correctly read into a 2D vector color map.
// Confirms that the encoded pixel values (0–3) match the original RGB pixels.
void debugColorMapData(const string& file_name, const vector<vector<uint8_t>>& color_map) {
    ifstream file(file_name, ios::binary);
    if (!file) {
        DIE("Failed to reopen file ", file_name, "\n");
        return;
    }

    // Read BMP header
    uint8_t header[54];
    file.read(reinterpret_cast<char*>(header), 54);

    // Extract image info
    int width = *(int*)&header[18];
    int height = *(int*)&header[22];
    int data_offset = *(int*)&header[10];
    int row_padded = (width * 3 + 3) & (~3);

    // Verify dimensions match
    if (color_map.size() != height || color_map[0].size() != width) {
        DIE("ERROR: Image dimensions don't match BMP header.\n");
        return;
    }

    // Seek to pixel data
    file.seekg(data_offset);

    // Stats
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

            // Accumulate total pixel brightness
            total_original_pixel_value += red + green + blue;

            // Determine the expected encoding
            uint8_t expected_value = 0;
            if (red == 255 && green == 0 && blue == 0) {
                expected_value = 1;
                count_red++;
            } else if (green == 255 && red == 0 && blue == 0) {
                expected_value = 2;
                count_green++;
            } else if (blue == 255 && red == 0 && green == 0) {
                expected_value = 3;
                count_blue++;
            } else if (red == 0 && green == 0 && blue == 0) {
                expected_value = 0;
                count_black++;
            } else {
                DIE("Unexpected pixel value at (", i, ", ", j, "): ",
                    "R=", red, " G=", green, " B=", blue, "\n");
            }

            // Compare against the decoded color map
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

    // Print summary statistics
    SLOG("BMP Color Map Statistics:\n",
         "  Dimensions: ", width, " x ", height, "\n",
         "  Red Pixels (1):   ", count_red, "\n",
         "  Green Pixels (2): ", count_green, "\n",
         "  Blue Pixels (3):  ", count_blue, "\n",
         "  Black Pixels (0): ", count_black, "\n",
         "  Total Mismatches: ", mismatches, "\n",
         "  Total raw pixel value sum:     ", total_original_pixel_value, "\n",
         "  Total mapped value sum:        ", total_mapped_value_sum, "\n");

    if (mismatches == 0) {
        SLOG("✅ Color map matches original BMP pixel data.\n");
    } else {
        DIE("❌ Found ", mismatches, " mismatches in pixel data.\n");
    }

    // Verify that value sum matches what we expect from encoding
    if (total_mapped_value_sum == count_red * 1 + count_green * 2 + count_blue * 3) {
        SLOG("✅ Encoded value sum matches expected pixel encoding.\n");
    } else {
        DIE("❌ Mismatch in total encoded value sum.\n");
    }
}


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

void dump_single_file(const string &fname, const string &out_str) {
  auto sz = out_str.length();
  upcxx::atomic_domain<size_t> ad({upcxx::atomic_op::fetch_add, upcxx::atomic_op::load});
  upcxx::global_ptr<size_t> fpos = nullptr;
  // always give rank 0 the first chunk so it can write header info if needed
  if (!upcxx::rank_me()) fpos = upcxx::new_<size_t>(sz);
  fpos = upcxx::broadcast(fpos, 0).wait();
  size_t my_fpos = 0;
  if (upcxx::rank_me()) my_fpos = ad.fetch_add(fpos, sz, std::memory_order_relaxed).wait();
  // wait until all ranks have updated the global counter
  upcxx::barrier();
  int fileno = -1;
  size_t fsize = 0;
  if (!upcxx::rank_me()) {
    fsize = ad.load(fpos, std::memory_order_relaxed).wait();
    // rank 0 creates the file and truncates it to the correct length
    fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) WARN("Error trying to create file ", fname, ": ", strerror(errno), "\n");
    if (ftruncate(fileno, fsize) == -1)
      WARN("Could not truncate ", fname, " to ", fsize, " bytes\n");
  }
  upcxx::barrier();
  ad.destroy();
  // wait until rank 0 has finished setting up the file
  if (rank_me()) fileno = open(fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) WARN("Error trying to open file ", fname, ": ", strerror(errno), "\n");
  auto bytes_written = pwrite(fileno, out_str.c_str(), sz, my_fpos);
  close(fileno);
  if (bytes_written != sz)
    DIE("Could not write all ", sz, " bytes; only wrote ", bytes_written, "\n");
  upcxx::barrier();
  auto tot_bytes_written = upcxx::reduce_one(bytes_written, upcxx::op_fast_add, 0).wait();
  SLOG_VERBOSE("Successfully wrote ", get_size_str(tot_bytes_written), " to ", fname, "\n");
}

std::shared_ptr<Random> _rnd_gen;
