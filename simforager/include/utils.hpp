#pragma once
#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <random>

#include <opencv2/opencv.hpp>

#include "upcxx_utils/log.hpp"

using namespace upcxx_utils;

using std::min;
using std::string;
using std::string_view;
using std::to_string;

#ifdef USE_BYTELL
#include "bytell_hash_map.hpp"
#define HASH_TABLE ska::bytell_hash_map
#else
#include <unordered_map>
#define HASH_TABLE std::unordered_map
#endif

int pin_thread(pid_t pid, int cid);

void dump_single_file(const string &fname, const string &out_str);

cv::Mat render_frame(
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &fish_points,
    int scale  // <- This makes it flexible and backward-compatible
);

// Write a full frame to the MP4 video file
void write_full_frame_to_video(
    const std::string &video_path,
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &fish_points);

// Finalize and close the video file
void finalize_video_writer();

class Random {
 private:
  std::mt19937_64 generator;

  double get_prob(double max_val = 1.0) {
    return std::uniform_real_distribution<>(0, max_val)(generator);
  }

 public:
  Random(unsigned seed)
      : generator(seed) {}

  int get(int64_t begin, int64_t end) {
    return std::uniform_int_distribution<int64_t>(begin, end - 1)(generator);
  }

  bool trial_success(double thres) {
    assert(thres >= 0);
    if (thres > 1) return true;
    if (thres == 0) return false;
    return (get_prob() <= thres);
  }

  int get_normal(vector<int> dist_params) {
    return (int)std::normal_distribution<float>(dist_params[0], dist_params[1])(generator);
  }

  int get_poisson(int avg) { return (int)std::poisson_distribution<int>(avg)(generator); }
};

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <string>

using namespace std;

// Reads a 24-bit BMP file and returns a 2D vector of encoded pixel values:
// 0 = Black, 1 = Red, 2 = Green, 3 = Blue.
// Assumes all RGB values in the BMP are either 0 or 255.
std::vector<std::vector<uint8_t>> readBMPColorMap(const std::string& file_name);

// Verifies that the encoded color map matches the original BMP file.
// Logs detailed statistics and errors if mismatches are found.
void debugColorMapData(const std::string& file_name,
                       const std::vector<std::vector<uint8_t>>& color_map);

// Writes a 2D vector of encoded pixel values (0 = Black, 1 = Red, 2 = Green, 3 = Blue)
// to a 24-bit BMP file. Output BMP will use only RGB values 0 or 255.
// The image is written bottom-up to match BMP format.
void writeBMPColorMap(const std::string& file_name,
                      const std::vector<std::vector<uint8_t>>& substrate_array);

extern std::shared_ptr<Random> _rnd_gen;

// Test function to see if we can write an MP4
void write_test_video(const std::string& output_path);
