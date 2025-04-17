
// main.cpp — Entry Point for SimForager
// -------------------------------------
// This file sets up and runs the core simulation loop.
// It initializes the reef, generates agents (fishes), seeds infections,
// performs updates to substrates and chemokines, and logs statistics.
// Uses UPC++ for distributed parallelism.
//
// Authors:
//   Steven Hofmeyr, LBNL (2020)
//   Matthew Fricke, UNM (2025)

#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <upcxx/upcxx.hpp>
#include <filesystem> // For output path management

using namespace std;

#include "options.hpp"
#include "reef.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"
#include "vonmises.hpp"

using namespace upcxx;
using namespace upcxx_utils;

// Short alias for getting current timestamp
#define NOW chrono::high_resolution_clock::now
#define STATS_COL_WIDTH 11

// Class to manage global simulation statistics
class SimStats {
 private:
  ofstream log_file;

 public:
  // Various counters tracked at each time step
  int64_t incubating = 0;
  int64_t expressing = 0;
  int64_t apoptotic = 0;
  int64_t dead = 0;
  int64_t fishes_vasculature = 0;
  int64_t fishes_reef = 0;
  float chemokines = 0;
  int64_t num_chemo_pts = 0;
  float floating_algaes = 0;

  // Called once to open output file
  void init() {
    if (!rank_me()) {
      log_file.open(_options->output_dir + "/simforager.stats");
      log_file << "# time" << header(0) << endl;
    }
  }

  // Constructs column headers
  string header(int width) {
    vector<string> columns = {"incb", "expr", "apop", "dead", "tvas", "ttis", "chem", "algae", "chempts", "%infct"};
    ostringstream oss;
    oss << left;
    for (auto column : columns) {
      if (width)
        oss << setw(width) << column;
      else
        oss << '\t' << column;
    }
    return oss.str();
  }

  // Collects data and formats current stats into tabular string
  string to_str(int width) {
    vector<int64_t> totals;
    totals.push_back(reduce_one(incubating, op_fast_add, 0).wait());
    totals.push_back(reduce_one(expressing, op_fast_add, 0).wait());
    totals.push_back(reduce_one(apoptotic, op_fast_add, 0).wait());
    totals.push_back(reduce_one(dead, op_fast_add, 0).wait());
    totals.push_back(reduce_one(fishes_vasculature, op_fast_add, 0).wait());
    totals.push_back(reduce_one(fishes_reef, op_fast_add, 0).wait());
    vector<float> totals_d;
    totals_d.push_back(reduce_one(chemokines, op_fast_add, 0).wait() / get_num_grid_points());
    totals_d.push_back(reduce_one(floating_algaes, op_fast_add, 0).wait());
    auto all_chem_pts = reduce_one(num_chemo_pts, op_fast_add, 0).wait();
    totals_d.push_back(all_chem_pts + totals[0] + totals[1] + totals[2] + totals[3]);
    auto perc_infected = 100.0 * (float)(totals[0] + totals[1] + totals[2] + totals[3]) / get_num_grid_points();

    ostringstream oss;
    oss << left;
    for (auto tot : totals) {
      if (width)
        oss << setw(width) << tot;
      else
        oss << '\t' << tot;
    }
    oss << fixed << setprecision(2) << scientific;
    for (auto tot : totals_d) {
      if (width)
        oss << setw(width) << tot;
      else
        oss << '\t' << tot;
    }
    oss << fixed << setprecision(2) << showpoint;
    if (width)
      oss << setw(width) << perc_infected;
    else
      oss << '\t' << perc_infected;
    return oss.str();
  }

  // Log output at current time step
  void log(int time_step) {
    string s = to_str(0);
    if (!rank_me()) log_file << time_step << s << endl;
  }
};

// Global instances used throughout the simulation
ofstream _logstream;
bool _verbose = false;
SimStats _sim_stats;
shared_ptr<Options> _options;

// Timers used to profile different parts of the simulation
IntermittentTimer generate_fish_timer(__FILENAME__ + string(":") + "generate fishes");
IntermittentTimer update_circulating_fishes_timer(__FILENAME__ + string(":") + "update circulating fishes");
IntermittentTimer update_fish_timer(__FILENAME__ + string(":") + "update fishes");
IntermittentTimer update_substrate_timer(__FILENAME__ + string(":") + "update substrates");
IntermittentTimer update_concentration_timer(__FILENAME__ + string(":") + "update concentrations");
IntermittentTimer compute_updates_timer(__FILENAME__ + string(":") + "compute updates");
IntermittentTimer accumulate_concentrations_timer(__FILENAME__ + string(":") + "dispatch updates");
IntermittentTimer add_new_actives_timer(__FILENAME__ + string(":") + "add new actives");
IntermittentTimer set_active_points_timer(__FILENAME__ + string(":") + "erase inactive");
IntermittentTimer sample_timer(__FILENAME__ + string(":") + "sample");
IntermittentTimer sample_write_timer(__FILENAME__ + string(":") + "sample write");
IntermittentTimer log_timer(__FILENAME__ + string(":") + "log");

// Additional simulation logic follows in the file...

// Note: This file is long, and the annotated version is delivered in parts.

