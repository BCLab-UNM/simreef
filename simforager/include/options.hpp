#pragma once

#include <sys/stat.h>

#include <iostream>
#include <regex>
#include <algorithm>
#include <upcxx/upcxx.hpp>

#include "CLI11.hpp"
#include "version.h"
#include "utils.hpp"

using std::array;
using std::cout;
using std::endl;
using std::sort;
using std::vector;

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/timers.hpp"

using namespace upcxx_utils;

#define YES_NO(X) ((X) ? "YES" : "NO")

class Options {
  vector<string> splitter(string in_pattern, string &content) {
    vector<string> split_content;
    std::regex pattern(in_pattern);
    copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1),
         std::sregex_token_iterator(), back_inserter(split_content));
    return split_content;
  }

  template <typename T>
  string vec_to_str(const vector<T> &vec, const string &delimiter = ",") {
    std::ostringstream oss;
    for (auto elem : vec) {
      oss << elem;
      if (elem != vec.back()) oss << delimiter;
    }
    return oss.str();
  }

  void setup_output_dir() {
    if (!upcxx::rank_me()) {
      // create the output directory and stripe it
      if (mkdir(output_dir.c_str(), S_IRWXU) == -1) {
        // could not create the directory
        if (errno == EEXIST) {
          cerr << KLRED << "WARNING: " << KNORM << "Output directory " << output_dir
               << " already exists. May overwrite existing files\n";
        } else {
          ostringstream oss;
          oss << KLRED << "ERROR: " << KNORM << " Could not create output directory " << output_dir
              << ": " << strerror(errno) << endl;
          throw std::runtime_error(oss.str());
        }
      } else {
        // created the directory - now stripe it if possible
        auto status = std::system("which lfs 2>&1 > /dev/null");
        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
          string cmd = "lfs setstripe -c 24 " + output_dir;
          auto status = std::system(cmd.c_str());
          if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
            cout << "Set Lustre striping on the output directory\n";
          else
            cout << "Failed to set Lustre striping on output directory: " << WEXITSTATUS(status)
                 << endl;
        }
      }
    }
    upcxx::barrier();
  }
    // now setup a samples subdirectory
    /*if (!upcxx::rank_me()) {
      // create the output directory and stripe it
      string samples_dir = output_dir + "/samples";
      if (mkdir(samples_dir.c_str(), S_IRWXU) == -1) {
        // could not create the directory
        if (errno == EEXIST) {
          cerr << KLRED << "WARNING: " << KNORM
               << "Samples directory already exists. May overwrite existing files\n";
        } else {
          ostringstream oss;
          oss << KLRED << "ERROR: " << KNORM
              << " Could not create samples directory : " << strerror(errno) << endl;
          throw std::runtime_error(oss.str());
        }
      }
    }
    */
    //upcxx::barrier();
  

  void setup_log_file() {
    if (!upcxx::rank_me()) {
      string log_fname = output_dir + "/simforager.log";
      // check to see if simforager.log exists. If so, rename it
      if (file_exists(log_fname)) {
        string new_log_fname = output_dir + "/simforager-" + get_current_time(true) + ".log";
        cerr << KLRED << "WARNING: " << KNORM << log_fname << " exists. Renaming to "
             << new_log_fname << endl;
        if (rename(log_fname.c_str(), new_log_fname.c_str()) == -1)
          DIE("Could not rename ", log_fname, ": ", strerror(errno));
      }
    }
    upcxx::barrier();
  }

  void set_random_infections(int num) {
    for (int i = 0; i < num; i++) {
      if (i % rank_n() != rank_me()) continue;
      infection_coords.push_back(
          {_rnd_gen->get(dimensions[0] * 0.1, dimensions[0] * 0.9),
           _rnd_gen->get(dimensions[1] * 0.1, dimensions[1] * 0.9),
           _rnd_gen->get(dimensions[2] > 1 ? dimensions[2] * 0.1 : 0,
                         dimensions[2] > 1 ? dimensions[2] * 0.9 : dimensions[2]),
           0});
    }
  }

  void set_uniform_infections(int num) {
    vector<array<int, 3>> infections =
        get_uniform_infections(num, dimensions[0], dimensions[1], dimensions[2]);
    for (int i = 0; i < infections.size(); i++) {
      if (i % rank_n() != rank_me()) continue;
      infection_coords.push_back({infections[i][0], infections[i][1], infections[i][2], 0});
    }
  }

  vector<array<int, 3>> get_uniform_infections(int num, int64_t dim_x, int64_t dim_y,
                                               int64_t dim_z) {
    vector<array<int, 3>> infections;
    int x_splits = 1, y_splits = 1, z_splits = 1;
    while (x_splits * y_splits * z_splits < num) {
      double x_ratio = (double)dim_x / x_splits;
      double y_ratio = (double)dim_y / y_splits;
      double z_ratio = (double)dim_z / z_splits;
      double ratios[] = {x_ratio, y_ratio, z_ratio};
      sort(ratios, ratios + 3);
      if (ratios[2] == x_ratio) {
        x_splits++;
      } else if (ratios[2] == y_ratio) {
        y_splits++;
      } else {
        z_splits++;
      }
    }
    int x_spacing = (double)dim_x / (x_splits + 1);
    int y_spacing = (double)dim_y / (y_splits + 1);
    int z_spacing = (double)dim_z / (z_splits + 1);
    if (dim_z == 1) {
      for (int i = x_spacing; i < dim_x - 1; i += x_spacing) {
        for (int j = y_spacing; j < dim_y - 1; j += y_spacing) {
          if (infections.size() == num) return infections;
          infections.push_back({i, j, 0});
        }
      }
    } else {
      for (int i = x_spacing; i < dim_x - 1; i += x_spacing) {
        for (int j = y_spacing; j < dim_y - 1; j += y_spacing) {
          for (int k = z_spacing; k < dim_z - 1; k += z_spacing) {
            if (infections.size() == num) return infections;
            infections.push_back({i, j, k});
          }
        }
      }
    }
    return infections;
  }

  bool parse_infection_coords(vector<string> &coords_strs) {
    auto get_locations_count = [](const string &s, const string &name) -> int {
      int num = 0;
      if (s.compare(0, name.length(), name) == 0) {
        string num_str = s.substr(name.length());
        try {
          num = std::stoi(num_str);
        } catch (std::invalid_argument arg) {
          num = 0;
        }
        if (num < 1) return 0;
      }
      return num;
    };

    if (coords_strs.size() == 1) {
      int num = get_locations_count(coords_strs[0], "random:");
      if (num > 0) {
        set_random_infections(num);
        return true;
      }
      num = get_locations_count(coords_strs[0], "uniform:");
      if (num > 0) {
        set_uniform_infections(num);
        return true;
      }
    }
    for (int i = 0; i < coords_strs.size(); i++) {
      if (i % rank_n() != rank_me()) continue;
      auto coords_and_time = splitter(",", coords_strs[i]);
      if (coords_and_time.size() == 4) {
        try {
          infection_coords.push_back({std::stoi(coords_and_time[0]), std::stoi(coords_and_time[1]),
                                      std::stoi(coords_and_time[2]),
                                      std::stoi(coords_and_time[3])});
        } catch (std::invalid_argument arg) {
          coords_and_time.clear();
        }
      }
      if (coords_and_time.size() != 4) {
        ostringstream oss;
        oss << KLRED << "ERROR: " << KNORM << "incorrect specification of infection coords in "
            << "string \"" << coords_strs[i] << "\"\n"
            << " - should be four comma-separated (x,y,z,t) values or random:N or uniform:N\n";
        cerr << oss.str();
        return false;
      }
    }
    return true;
  }

  vector<int> get_model_dims(const string &fname) {
    ifstream f(fname, std::ios::in | std::ios::binary);
    if (!f) {
      cerr << "Couldn't open file " << fname << endl;
      return {};
    }
    vector<int> dims(3);
    if (!f.read(reinterpret_cast<char *>(&(dims[0])), 3 * sizeof(int))) {
      cerr << "Couldn't read dims in " << fname << endl;
      return {};
    }
    return dims;
  }

 public:
  vector<int> dimensions{300, 300, 1};
  vector<int> ecosystem_dims{48000, 40000, 20000};
  // each time step should be about 1 minute, so one day = 1440 time steps
  int num_timesteps = 20160;

  // x,y,z location and timestep
  vector<array<int, 4>> infection_coords;
  int initial_infection = 1000;

  string ecosystem_model_dir = "";

  // these periods are normally distributed with mean and stddev
  int incubation_period = 480;
  int apoptosis_period = 180;
  int expressing_period = 2286;

  string substrate_bitmap_path = "rgb.bmp";
  int fish_generation_rate = 100000;
  int fish_initial_delay = 10080;
  int fish_vascular_period = 5760;
  int fish_reef_period = 1440;
  int fish_binding_period = 10;
  double max_binding_prob = 1.0;

  double infectivity = 0.02;
  double infectivity_multiplier = 1.0;
  double floating_algae_production = 35;
  double floating_algae_production_multiplier = 1.0;
  double floating_algae_clearance_rate = 0.002;
  double floating_algae_diffusion_coef = 1.0;

  double chemokine_production = 1.0;
  double chemokine_decay_rate = 0.01;
  double chemokine_diffusion_coef = 1.0;
  double min_chemokine = 1e-6;

  double antibody_factor = 1;
  int antibody_period = 5760;

  unsigned rnd_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  string output_dir = "simforager-results-n" + to_string(upcxx::rank_n()) + "-N" +
                      to_string(upcxx::rank_n() / upcxx::local_team().rank_n());
  int sample_period = 0;
  int sample_resolution = 1;
  int max_block_dim = 10;

  bool fishes_follow_gradient = false;

  int num_fish = 1;
  int predator_ratio = 0;

  int predator_detection_radius = 0;
  int grazer_detection_radius = 0;

  // Grazers: κ values with predator 
  int kappa_grazer_w_predator_coral_w_algae = 0;
  int kappa_grazer_w_predator_coral_no_algae = 0;
  int kappa_grazer_w_predator_sand_w_algae = 0;
  int kappa_grazer_w_predator_sand_no_algae = 0;
  
  // Grazers: κ values without predator
  int kappa_grazer_wo_predator_coral_w_algae = 0;
  int kappa_grazer_wo_predator_coral_no_algae = 0;
  int kappa_grazer_wo_predator_sand_w_algae = 0;
  int kappa_grazer_wo_predator_sand_no_algae = 0;
  
  // Grazers: step lengths with predator
  int step_len_grazer_w_predator_coral_w_algae = 0;
  int step_len_grazer_w_predator_coral_no_algae = 0;
  int step_len_grazer_w_predator_sand_w_algae = 0;
  int step_len_grazer_w_predator_sand_no_algae = 0;
  
  // Grazers: step lengths without predator
  int step_len_grazer_wo_predator_coral_w_algae = 0;
  int step_len_grazer_wo_predator_coral_no_algae = 0;
  int step_len_grazer_wo_predator_sand_w_algae = 0;
  int step_len_grazer_wo_predator_sand_no_algae = 0;
  
  
  // Predators: κ values with grazer
  int kappa_predator_w_grazer_coral_w_algae = 0;
  int kappa_predator_w_grazer_coral_no_algae = 0;
  int kappa_predator_w_grazer_sand_w_algae = 0;
  int kappa_predator_w_grazer_sand_no_algae = 0;
  
  // Predators: κ values without grazer
  int kappa_predator_wo_grazer_coral_w_algae = 0;
  int kappa_predator_wo_grazer_coral_no_algae = 0;
  int kappa_predator_wo_grazer_sand_w_algae = 0;
  int kappa_predator_wo_grazer_sand_no_algae = 0;
  
  // Predators: step lengths with grazer
  int step_len_predator_w_grazer_coral_w_algae = 0;
  int step_len_predator_w_grazer_coral_no_algae = 0;
  int step_len_predator_w_grazer_sand_w_algae = 0;
  int step_len_predator_w_grazer_sand_no_algae = 0;
  
  // Predators: step lengths without grazer
  int step_len_predator_wo_grazer_coral_w_algae = 0;
  int step_len_predator_wo_grazer_coral_no_algae = 0;
  int step_len_predator_wo_grazer_sand_w_algae = 0;
  int step_len_predator_wo_grazer_sand_no_algae = 0;
  

  int log_grazer_tracks = 0;
  int log_population_stat = 0;

  // Amount of attached algae to seed on every ALGAE cell
  double algae_init_count = 100.0;

  // How much a grazer eats from the current cell per timestep
  double algae_grazing_rate = 1.0;

  // Optional: when algae is fully eaten, convert the cell to SAND
  //bool algae_turns_to_coral_when_depleted = false;



  

  bool show_progress = false;
  bool verbose = false;

  bool load(int argc, char **argv) {
    // SIMCOV version v0.1-a0decc6-master (Release) built on 2020-04-08T22:15:40 with g++
    string full_version_str = "SimForager version " + string(SIMFORAGER_VERSION) + "-" +
                              string(SIMFORAGER_BRANCH) + " built on " + string(SIMFORAGER_BUILD_DATE);
    vector<string> infection_coords_strs;

    CLI::App app(full_version_str);
    app.add_option("-d,--dim", dimensions, "Dimensions: x y z")
        ->delimiter(',')
        ->expected(3)
        ->capture_default_str();
    app.add_option("--ecosystem-dim", ecosystem_dims, "Ecosystem dimensions: x y z")
        ->delimiter(',')
        ->expected(3)
        ->capture_default_str();
    app.add_option("-t,--timesteps", num_timesteps, "Number of timesteps")
        ->check(CLI::Range(1, 1000000))
        ->capture_default_str();
    app.add_option(
           "--infection-coords", infection_coords_strs,
           "Location of multiple initial infections, of form \"x1,y1,z1,t1 x2,y2,z2,t2...\"\n"
           "where x,y,z are grid coords and t is a timestep, or\n"
           "\"uniform:N\" for N uniformly distributed points or\n"
           "\"random:N\" for N randomly distributed points")
        ->delimiter(' ')
        ->capture_default_str();
    app.add_option("--initial-infection", initial_infection,
                   "Number of floating_algaes at initial infection locations")
        ->capture_default_str();
    app.add_option("--substrate-bitmap", substrate_bitmap_path,
                   "Path to substrate bitmap. Red=coral, blue=sand, green=algae.")
        ->capture_default_str();
    app.add_option("--incubation-period", incubation_period,
                   "Average number of time steps to expressing floating_algaes after cell is infected")
        ->capture_default_str();
    app.add_option("--apoptosis-period", apoptosis_period,
                   "Average number of time steps to death after apoptosis is induced")
        ->capture_default_str();
    app.add_option("--expressing-period", expressing_period,
                   "Average number of time steps to death after a cell starts expresssing")
        ->capture_default_str();
    app.add_option("--infectivity", infectivity,
                   "Factor multiplied by number of floating_algaes to determine probability of infection")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--infectivity-multiplier", infectivity_multiplier,
                   "Multiplier to reduce infectivity rate")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--floating_algae-production", floating_algae_production,
                   "Number of floating_algaes produced by expressing cell each time step")
        ->capture_default_str();
    app.add_option("--floating_algae-production-multiplier", floating_algae_production_multiplier,
                   "Multiplier to reduce floating_algae production rate")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--floating_algae-clearance", floating_algae_clearance_rate,
                   "Fraction by which floating_algae count drops each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--floating_algae-diffusion", floating_algae_diffusion_coef,
                   "Fraction of floating_algaes that diffuse into all neighbors each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--chemokine-production", chemokine_production,
                   "Amount of chemokine produced by expressing cells each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--chemokine-decay", chemokine_decay_rate,
                   "Amount by which chemokine concentration drops each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--chemokine-diffusion", chemokine_diffusion_coef,
                   "Fraction of chemokine concentration that diffuses into all neighbors "
                   "each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--min-chemokine", min_chemokine,
                   "Minimum chemokine concentration that triggers a fish")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--antibody-factor", antibody_factor,
                   "Impact of antibodies; multiplier for floating_algae clearance")
        ->capture_default_str();
    app.add_option("--antibody-period", antibody_period,
                   "Number of time steps before antibodies start to be produced")
        ->capture_default_str();
    app.add_option("--fish-generation-rate", fish_generation_rate,
                   "Number of fishes generated at each timestep for the whole ecosystem")
        ->capture_default_str();
    app.add_option("--fish-initial-delay", fish_initial_delay,
                   "Number of time steps before fishes start to be produced")
        ->capture_default_str();
    app.add_option("--fish-vascular-period", fish_vascular_period,
                   "Average number of time steps to death for a fish in the vasculature")
        ->capture_default_str();
    app.add_option("--fish-reef-period", fish_reef_period,
                   "Average number of time steps to death after a fish extravasates")
        ->capture_default_str();
    app.add_option("--fish-binding-period", fish_binding_period,
                   "Number of time steps a fish is bound to an epithelial cell when inducing "
                   "apoptosis")
        ->capture_default_str();
    app.add_option("--max-binding-prob", max_binding_prob,
                   "Max probability of a fish binding to an infected cell in one time step")
        ->capture_default_str();
    app.add_flag("--fishes-follow-gradient", fishes_follow_gradient,
                 "fishes in reef follow the chemokine gradient")
        ->capture_default_str();
    app.add_option("-r,--seed", rnd_seed, "Random seed")->capture_default_str();
    app.add_option("--sample-period", sample_period,
                   "Number of timesteps between samples (set to 0 to disable sampling)")
        ->capture_default_str();
    app.add_option("--sample-resolution", sample_resolution, "Resolution for sampling")
        ->check(CLI::Range(1, 10000))
        ->capture_default_str();
    app.add_option("--max-block-dim", max_block_dim,
                   "Max. block dimension - larger means more locality but worse load balance. Set "
                   "to 0 for largest possible")
        ->capture_default_str();
    app.add_option("-o,--output", output_dir, "Output directory")->capture_default_str();
    app.add_option("--ecosystem-model", ecosystem_model_dir, "Directory containing files for ecosystem model")
        ->capture_default_str();
    app.add_option("-n, --num-fish", num_fish, "Number of fish to generate, default = 1")
    ->capture_default_str();
    app.add_option("--predator-ratio", predator_ratio, "Ratio of fish that are predators, default = 0")
    ->capture_default_str();

    app.add_option("--predator-detect-grazer-radius", predator_detection_radius, "radius at which a predator can detect the presence of a grazer, default = 0")
    ->capture_default_str();
        app.add_option("--grazer-detect-predator-radius", grazer_detection_radius, "radius at which a grazer can detect the presence of a predator, default = 0")
    ->capture_default_str();
    
    app.add_option("--kappa_grazer_wo_predator_coral_w_algae", kappa_grazer_wo_predator_coral_w_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_grazer_wo_predator_coral_no_algae", kappa_grazer_wo_predator_coral_no_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_grazer_wo_predator_sand_w_algae", kappa_grazer_wo_predator_sand_w_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_grazer_wo_predator_sand_no_algae", kappa_grazer_wo_predator_sand_no_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();

    // Grazer without predator, step lengths
        app.add_option("--step_length_grazer_wo_predator_coral_w_algae", step_len_grazer_wo_predator_coral_w_algae, "Step length for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--step_length_grazer_wo_predator_coral_no_algae", step_len_grazer_wo_predator_coral_no_algae, "Step length for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--step_length_grazer_wo_predator_sand_w_algae", step_len_grazer_wo_predator_sand_w_algae, "Step length for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--step_length_grazer_wo_predator_sand_no_algae", step_len_grazer_wo_predator_sand_no_algae, "Step length value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    
    
    app.add_option("--kappa_grazer_w_predator_coral_w_algae", kappa_grazer_w_predator_coral_w_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_grazer_w_predator_coral_no_algae", kappa_grazer_w_predator_coral_no_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_grazer_w_predator_sand_w_algae", kappa_grazer_w_predator_sand_w_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_grazer_w_predator_sand_no_algae", kappa_grazer_w_predator_sand_no_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();

    app.add_option("--kappa_predator_wo_grazer_coral_w_algae", kappa_predator_wo_grazer_coral_w_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_predator_wo_grazer_coral_no_algae", kappa_predator_wo_grazer_coral_no_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_predator_wo_grazer_sand_w_algae", kappa_predator_wo_grazer_sand_w_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_predator_wo_grazer_sand_no_algae", kappa_predator_wo_grazer_sand_no_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();


    app.add_option("--kappa_predator_w_grazer_coral_w_algae", kappa_predator_w_grazer_coral_w_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_predator_w_grazer_coral_no_algae", kappa_predator_w_grazer_coral_no_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_predator_w_grazer_sand_w_algae", kappa_predator_w_grazer_sand_w_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();
    app.add_option("--kappa_predator_w_grazer_sand_no_algae", kappa_predator_w_grazer_sand_no_algae, "Kappa value for Von Mises correlated random walk over sand, default = 0")
    ->capture_default_str();

    app.add_option("--algae-initial-count", algae_init_count, "Amount of attached algae to seed on every ALGAE grid point, default = 100")
    ->capture_default_str();

    app.add_option("--algae-decomp-rate-from-grazing", algae_grazing_rate, "Grazer consumption from the current cell per timestep")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();

    //app.add_option("--algae-turns-to-coral-when-depleted", algae_turns_to_coral_when_depleted, "when algae is fully eaten, convert the cell to SAND/something else, default = false")
    //->capture_default_str();

    
    app.add_option("--log_grazer_tracks", log_grazer_tracks, "Log the grazers' path for tracking: 0 = no , 1 = yes")->capture_default_str();
    app.add_option("--log_population_stat", log_population_stat, "Log the grazer population statistics at the end : 0 = no , 1 = yes")->capture_default_str();



    app.add_flag("--progress", show_progress, "Show progress");
    app.add_flag("-v, --verbose", verbose, "Verbose output");

    auto *cfg_opt = app.set_config("--config", "", "Load options from a configuration file");

    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
      if (upcxx::rank_me() == 0) app.exit(e);
      return false;
    }

    upcxx::barrier();

    _rnd_gen = make_shared<Random>(rnd_seed + rank_me());
    
    if (!ecosystem_model_dir.empty()) {
      auto model_dims = get_model_dims(ecosystem_model_dir + "/alveolus.dat");
      if (model_dims.size() < 3) return false;
      for (int i = 0; i < 3; i++) {
        if (model_dims[i] != dimensions[i]) {
          dimensions = model_dims;
          if (!rank_me())
            cerr << KLRED << "WARNING: " << KNORM
                 << "Setting dimensions to model data: " << dimensions[0] << ", " << dimensions[1]
                 << ", " << dimensions[2] << endl;
          break;
        }
      }
    }

    if (floating_algae_clearance_rate * antibody_factor > 1.0) {
      if (!rank_me())
        cerr << "Invalid parameter settings: floating_algae-clearance * antibody_factor > 1.\n"
             << "Reduce either or both of those settings\n";
      return false;
    }

    if (!infection_coords_strs.empty() && !parse_infection_coords(infection_coords_strs))
      return false;

    if (!max_block_dim) {
      max_block_dim = min(dimensions[0], dimensions[1]);
      if (dimensions[2] > 1) max_block_dim = min(dimensions[2], max_block_dim);
    }

    if (dimensions[0] % sample_resolution || dimensions[1] % sample_resolution ||
        (dimensions[2] > 1 && dimensions[2] % sample_resolution)) {
      if (!rank_me())
        cerr << "Error: sample period " << sample_resolution
             << " must be a factor of all the dimensions\n";
      return false;
    }

    for (int i = 0; i < 3; i++) {
      if (dimensions[i] > ecosystem_dims[i]) {
        if (!rank_me()) cerr << "Dimensions must be <= whole ecosystem dimensions\n";
        return false;
      }
    }
    setup_output_dir();
    setup_log_file();

    init_logger(output_dir + "/simforager.log", verbose);

#ifdef DEBUG
    open_dbg("debug");
#endif

    SLOG(KLBLUE, "SimForager version ", full_version_str, KNORM, "\n");

    if (upcxx::rank_me() == 0) {
      // print out all compiler definitions
      SLOG_VERBOSE(KLBLUE, "_________________________", KNORM, "\n");
      SLOG_VERBOSE(KLBLUE, "Compiler definitions:", KNORM, "\n");
      std::istringstream all_defs_ss(ALL_DEFNS);
      vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)),
                              std::istream_iterator<string>());
      for (auto &def : all_defs) SLOG_VERBOSE("  ", def, "\n");
      SLOG_VERBOSE("_________________________\n");
      SLOG(KLBLUE, "Options:", KNORM, "\n");
      auto all_opts_strings = app.config_to_str_vector(true, false);
      for (auto &opt_str : all_opts_strings) SLOG(KLBLUE, opt_str, KNORM, "\n");
      SLOG(KLBLUE, "_________________________", KNORM, "\n");
    }
    auto num_nodes = upcxx::rank_n() / upcxx::local_team().rank_n();
    SLOG("Starting run with ", upcxx::rank_n(), " processes on ", num_nodes, " node",
         (num_nodes > 1 ? "s" : ""), " at ", get_current_time(), "\n");
#ifdef DEBUG
    SWARN("Running low-performance debug mode");
#endif
    if (!upcxx::rank_me()) {
      // write out configuration file for restarts
      ofstream ofs(output_dir + "/simforager.config");
      ofs << app.config_to_str(true, true);
    }
    upcxx::barrier();
    return true;
  }
};
