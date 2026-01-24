#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <unistd.h>
#include <vector>

#include <opencv2/opencv.hpp>

#include "upcxx_utils/log.hpp"
#include "vonmises.hpp"

#ifdef USE_BYTELL
  #include "bytell_hash_map.hpp"
  #define HASH_TABLE ska::bytell_hash_map
#else
  #include <unordered_map>
  #define HASH_TABLE std::unordered_map
#endif

// -----------------------------------------------------------------------------
// Forward declarations
// -----------------------------------------------------------------------------
class Options;
class Reef;
enum class SubstrateType;

// -----------------------------------------------------------------------------
// Logging namespace
// -----------------------------------------------------------------------------
using namespace upcxx_utils;

// -----------------------------------------------------------------------------
// Misc utilities
// -----------------------------------------------------------------------------
int pin_thread(pid_t pid, int cid);

void dump_single_file(const std::string &fname, const std::string &out_str);

// -----------------------------------------------------------------------------
// Video rendering helpers
// -----------------------------------------------------------------------------
cv::Mat render_frame(
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar, float, float, float, int>> &fish_points,
    int scale = 1);

/**
 * @brief Writes a full frame to the MP4 video using substrate and fish data.
 *
 * @param video_path Path to the output video file (MP4).
 * @param width Width of the frame in pixels.
 * @param height Height of the frame in pixels.
 * @param coral_w_algae_points Vector of points: (x, y, colour).
 * @param coral_no_algae_points Vector of points: (x, y, colour).
 * @param sand_w_algae_points Vector of points: (x, y, colour).
 * @param sand_no_algae_points Vector of points: (x, y, colour).
 * @param fish_points Vector of fish points:
 *        (x, y, colour, kappa, step_length, density, thickness)
 * @param scale Optional downscale factor for rendering.
 */
void write_full_frame_to_video(
    const std::string &video_path,
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar, float, float, float, int>> &fish_points,
    int scale);

double sigmoid(double x, double midpoint, double steepness);

// Finalise and close the video file
void finalize_video_writer();

// -----------------------------------------------------------------------------
// Random helper
// -----------------------------------------------------------------------------
class Random {
 private:
  std::mt19937_64 generator;

  double get_prob(double max_val = 1.0) {
    return std::uniform_real_distribution<double>(0.0, max_val)(generator);
  }

 public:
  explicit Random(unsigned seed)
      : generator(seed) {}

  int get(int64_t begin, int64_t end) {
    return static_cast<int>(std::uniform_int_distribution<int64_t>(begin, end - 1)(generator));
  }

  bool trial_success(double thres) {
    assert(thres >= 0.0);
    if (thres > 1.0) return true;
    if (thres == 0.0) return false;
    return (get_prob() <= thres);
  }

  int get_normal(std::vector<int> dist_params) {
    return static_cast<int>(std::normal_distribution<float>(dist_params[0], dist_params[1])(generator));
  }

  int get_poisson(int avg) {
    return static_cast<int>(std::poisson_distribution<int>(avg)(generator));
  }
};

extern std::shared_ptr<Random> _rnd_gen;

// Test function to see if we can write an MP4
void write_test_video(const std::string &output_path);

// -----------------------------------------------------------------------------
// BMP helpers
// -----------------------------------------------------------------------------
std::vector<std::vector<uint8_t>> readBMPColorMap(const std::string &file_name);

void debugColorMapData(
    const std::string &file_name,
    const std::vector<std::vector<uint8_t>> &color_map);

void writeBMPColorMap(
    const std::string &file_name,
    const std::vector<std::vector<uint8_t>> &substrate_array);

// -----------------------------------------------------------------------------
// Grazer trajectory logging
// -----------------------------------------------------------------------------
void log_grazer_step(
    const std::string &fish_id,
    int timestep,
    int64_t x,
    int64_t y,
    int64_t z,
    int substrate,
    float kappa);

void finalize_grazer_logs();

void log_algae_and_grazer_stats(
    int64_t total_timesteps,
    int64_t grazer_count,
    int64_t coral_w_algae_steps,
    int64_t coral_no_algae_steps,
    int64_t sand_w_algae_steps,
    int64_t sand_no_algae_steps,
    double initial_algae,
    double final_algae);


namespace utils {
std::pair<int64_t, int64_t> nearest_grid_point(double x, double y, const std::shared_ptr<GridCoords>& grid_size);
}

// -----------------------------------------------------------------------------
// Chamfer distance transform utilities
// -----------------------------------------------------------------------------
namespace utils {

static inline int grid_index(int x, int y, int W) {
  return y * W + x;
}
  
static inline uint16_t sat_add_u16(uint16_t a, uint16_t b, uint16_t INF) {
  if (a >= INF) return INF;
  const uint32_t s = static_cast<uint32_t>(a) + static_cast<uint32_t>(b);
  return (s >= INF) ? INF : static_cast<uint16_t>(s);
}

/**
 * 3x3 Chamfer Distance Transform (integer weights)
 *
 * Computes distance to a target label for every cell in a W x H grid.
 *
 * Weights (chamfer units):
 *   orthogonal step = 3
 *   diagonal step   = 4
 *
 * Approximate Euclidean distance in cell units:
 *   distance_cells ≈ D / 3.0f
 */
template <typename LabelT>
inline void chamfer_3x3_distance_transform(
    const std::vector<LabelT> &grid,   // length W*H, labels per cell
    int W,
    int H,
    const LabelT &target,
    std::vector<uint16_t> &D_out,
    uint16_t INF = static_cast<uint16_t>(std::numeric_limits<uint16_t>::max() / 4))
{
  constexpr uint16_t w_orth = 3;
  constexpr uint16_t w_diag = 4;

  D_out.assign(static_cast<size_t>(W) * static_cast<size_t>(H), INF);

  // Initialise: 0 for target cells, INF for all others.
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const int i = grid_index(x, y, W);
      if (grid[static_cast<size_t>(i)] == target) {
        D_out[static_cast<size_t>(i)] = 0;
      }
    }
  }

  // Forward pass (top left to bottom right): uses W, N, NW, NE
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const int i = grid_index(x, y, W);
      uint16_t d = D_out[static_cast<size_t>(i)];

      if (x > 0) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x - 1, y, W))], w_orth, INF));
      }
      if (y > 0) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x, y - 1, W))], w_orth, INF));
      }
      if (x > 0 && y > 0) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x - 1, y - 1, W))], w_diag, INF));
      }
      if (x + 1 < W && y > 0) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x + 1, y - 1, W))], w_diag, INF));
      }

      D_out[static_cast<size_t>(i)] = d;
    }
  }

  // Backward pass (bottom right to top left): uses E, S, SE, SW
  for (int y = H - 1; y >= 0; --y) {
    for (int x = W - 1; x >= 0; --x) {
      const int i = grid_index(x, y, W);
      uint16_t d = D_out[static_cast<size_t>(i)];

      if (x + 1 < W) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x + 1, y, W))], w_orth, INF));
      }
      if (y + 1 < H) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x, y + 1, W))], w_orth, INF));
      }
      if (x + 1 < W && y + 1 < H) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x + 1, y + 1, W))], w_diag, INF));
      }
      if (x > 0 && y + 1 < H) {
        d = std::min(d, sat_add_u16(D_out[static_cast<size_t>(grid_index(x - 1, y + 1, W))], w_diag, INF));
      }

      D_out[static_cast<size_t>(i)] = d;
    }
  }
}

// Convert chamfer units to approximate Euclidean distance in cells.
static inline float chamfer_to_cells(uint16_t D) {
  return static_cast<float>(D) / 3.0f;
}

/**
 * Write a debug PNG showing random per cell distance annotations for substrate targets.
 * Definition should be in utils.cpp as: void utils::showSubstrateDistances(...)
 */
void showSubstrateDistances(
    const Reef &reef,
    const Options &opt,
    double fraction_to_show,
    const std::string &out_png,
    int cell_px = 10,
    uint32_t rng_seed = 314159u);
  
void showSubstrateDistanceContours(
    const Reef& reef,
    const Options& opt,
    SubstrateType target,
    const std::string& out_png,
    float contour_step_cells,
    unsigned int rng_seed);

} // namespace utils
