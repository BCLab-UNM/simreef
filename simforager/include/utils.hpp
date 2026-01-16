#pragma once
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <string>

using namespace std;
#include <algorithm>
#include <random>
#include "vonmises.hpp"

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
    int width, int height,
    const std::vector<std::tuple<int,int,cv::Scalar>> &coral_w_algae_points,
    const std::vector<std::tuple<int,int,cv::Scalar>> &coral_no_algae_points,
    const std::vector<std::tuple<int,int,cv::Scalar>> &sand_w_algae_points,
    const std::vector<std::tuple<int,int,cv::Scalar>> &sand_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar, float, float, float, int>> &fish_points, // <-- 6 elements
    int scale = 1);

/**
 * @brief Writes a full frame to the MP4 video using substrate and fish data.
 *
 * @param video_path Path to the output video file (MP4).
 * @param width Width of the frame in pixels.
 * @param height Height of the frame in pixels.
 * @param coral_points Vector of coral points: (x, y, color).
 * @param algae_points Vector of algae points: (x, y, color).
 * @param sand_points Vector of sand points: (x, y, color).
 * @param fish_points Vector of fish points: (x, y, color, thickness).
 *        - thickness = -1 → filled circle
 *        - thickness > 0 → outline with given thickness
 * @param scale Optional downscale factor for rendering (default = 1).
 */
void write_full_frame_to_video(
    const std::string &video_path,
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar, float, float, float, int>> &fish_points, // <-- 5 elements
    int scale);

double sigmoid(double x, double midpoint, double steepness);

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

void log_grazer_step(const std::string& fish_id,int timestep, int64_t x, int64_t y, int64_t z, int substrate, float kappa);

void finalize_grazer_logs();

void log_algae_and_grazer_stats(int64_t total_timesteps,
                                int64_t grazer_count,
                                int64_t coral_w_algae_steps,
                                int64_t coral_no_algae_steps,
                                int64_t sand_w_algae_steps,
                                int64_t sand_no_algae_steps,
                                double initial_algae,
                                double final_algae);


// ------------------------------------------------------------
// 3x3 Chamfer Distance Transform (integer weights)
// ------------------------------------------------------------
//
// Computes distance-to-target for every cell in a W x H grid.
//
// Weights (chamfer units):
//   orthogonal step = 3
//   diagonal step   = 4
//
// Approximate Euclidean distance in cell units:
//   distance_cells ≈ D / 3.0f
//
// Notes:
// - This is an approximation to Euclidean distance - will be off by a couple of %
// - If the target class does not appear in the map, the output will remain INF everywhere.
// ------------------------------------------------------------

namespace utils {

static inline int grid_index(int x, int y, int W) {
  return y * W + x;
}

// Safe saturating add for uint16_t (prevents overflow on INF).
static inline uint16_t sat_add_u16(uint16_t a, uint16_t b, uint16_t INF) {
  if (a >= INF) return INF;
  uint32_t s = static_cast<uint32_t>(a) + static_cast<uint32_t>(b);
  return (s >= INF) ? INF : static_cast<uint16_t>(s);
}

template <typename LabelT>
inline void chamfer_3x3_distance_transform(
    const std::vector<LabelT>& grid,  // length W*H, labels per cell
    int W,
    int H,
    const LabelT& target,
    std::vector<uint16_t>& D_out,
    uint16_t INF = static_cast<uint16_t>(std::numeric_limits<uint16_t>::max() / 4))
{
  constexpr uint16_t w_orth = 3;
  constexpr uint16_t w_diag = 4;

  D_out.assign(static_cast<size_t>(W) * static_cast<size_t>(H), INF);

  // Initialise: 0 for target cells, INF for all others.
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      int i = grid_index(x, y, W);
      if (grid[static_cast<size_t>(i)] == target) {
        D_out[static_cast<size_t>(i)] = 0;
      }
    }
  }

  // Forward pass (top-left -> bottom-right): uses W, N, NW, NE
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      int i = grid_index(x, y, W);
      uint16_t d = D_out[static_cast<size_t>(i)];

      if (x > 0) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x-1, y, W))], w_orth, INF));
      }
      if (y > 0) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x, y-1, W))], w_orth, INF));
      }
      if (x > 0 && y > 0) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x-1, y-1, W))], w_diag, INF));
      }
      if (x + 1 < W && y > 0) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x+1, y-1, W))], w_diag, INF));
      }

      D_out[static_cast<size_t>(i)] = d;
    }
  }

  // Backward pass (bottom-right -> top-left): uses E, S, SE, SW
  for (int y = H - 1; y >= 0; --y) {
    for (int x = W - 1; x >= 0; --x) {
      int i = grid_index(x, y, W);
      uint16_t d = D_out[static_cast<size_t>(i)];

      if (x + 1 < W) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x+1, y, W))], w_orth, INF));
      }
      if (y + 1 < H) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x, y+1, W))], w_orth, INF));
      }
      if (x + 1 < W && y + 1 < H) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x+1, y+1, W))], w_diag, INF));
      }
      if (x > 0 && y + 1 < H) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x-1, y+1, W))], w_diag, INF));
      }

      D_out[static_cast<size_t>(i)] = d;
    }
  }
}

// Convert chamfer units to approximate Euclidean distance in cells.
static inline float chamfer_to_cells(uint16_t D) {
  return static_cast<float>(D) / 3.0f;
}

} // namespace simreef_utils
