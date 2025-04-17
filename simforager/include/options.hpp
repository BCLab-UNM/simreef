
// options.hpp — Runtime Configuration
// ------------------------------------
// Declares the Options struct that holds simulation parameters,
// including grid dimensions, diffusion constants, timing controls,
// and input/output settings.

#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <string>
#include <vector>
#include <array>
#include <memory>

using std::string;
using std::vector;
using std::array;
using std::shared_ptr;

enum class ViewObject {
  SUBSTRATE,
  FISH,
  CHEMOKINE,
  ALGAE
};

struct Options {
  // Simulation size and behavior
  array<int, 3> dimensions = {0, 0, 0};      // x, y, z size
  int max_block_dim = 1;                    // Used for block partitioning logic

  // Timing constants for biological processes
  int incubation_period = 5;
  int expressing_period = 5;
  int apoptosis_period = 5;

  // Binding behavior
  double max_binding_prob = 1.0;
  double infectivity = 0.1;
  double infectivity_multiplier = 2.0;

  // Algae production and chemokine generation
  double floating_algae_production = 0.1;
  double floating_algae_production_multiplier = 2.0;
  double chemokine_production = 0.1;

  // Diffusion settings
  double chemokine_diffusion_coef = 0.2;
  double floating_algae_diffusion_coef = 0.3;
  double chemokine_decay_rate = 0.1;
  double floating_algae_clearance_rate = 0.1;

  double min_chemokine = 0.01;

  // Fish simulation
  int fish_generation_rate = 5;
  int fish_binding_period = 5;
  int fish_vascular_period = 10;
  int initial_infection = 5;
  bool fishes_follow_gradient = true;
  double kappa = 1.0;  // Von Mises concentration

  // File input/output
  string output_dir;
  string ecosystem_model_dir;
  int sample_resolution = 10;

  vector<array<int, 4>> infection_coords;
};

extern shared_ptr<Options> _options;

#endif
