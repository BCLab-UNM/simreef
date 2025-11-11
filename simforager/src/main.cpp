// SimForager
//
// Steven Hofmeyr, LBNL May 2020, Matthew Fricke (mfricke@unm.edu) 8th March 2025

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

#include <filesystem> // For current directory

using namespace std;

#include "options.hpp"
#include "reef.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"

// For wrap_index function
#include <cmath>
#include <cstdint>

// For nearest neighbourhood
#include <utility>
#include <memory>

using namespace upcxx;
using namespace upcxx_utils;

#define NOW chrono::high_resolution_clock::now
#define STATS_COL_WIDTH 11

//static inline int64_t wrap_index(int64_t v, int64_t maxv) {
  // a true mathematical modulo for negatives
//  int64_t r = v % maxv;
//  return (r < 0) ? (r + maxv) : r;
//}

// Map continuous (x, y) to the nearest valid grid cell indices
static inline std::pair<int64_t, int64_t>
nearest_grid_point(double x, double y, const std::shared_ptr<GridCoords>& grid_size)
{
    // Round to nearest integer cell centre
    int64_t gx = static_cast<int64_t>(std::floor(x + 0.5));
    int64_t gy = static_cast<int64_t>(std::floor(y + 0.5));

    // Wrap indices into valid [0, size) range using modular arithmetic
    gx = (gx % grid_size->x + grid_size->x) % grid_size->x;
    gy = (gy % grid_size->y + grid_size->y) % grid_size->y;

    return {gx, gy};
}

static inline double wrap_index(double v, int64_t maxv) {
    double maxd = static_cast<double>(maxv);
    double r = std::fmod(v, maxd);
    if (r < 0.0) r += maxd;
    if (r >= maxd) r -= maxd;  // guard against rounding to exactly maxd
    return r;
}

// Random number generator
boost::random::mt19937 gen;

class SimStats {
 private:
  ofstream log_file;

 public:
  int64_t incubating = 0;
  int64_t expressing = 0;
  int64_t apoptotic = 0;
  int64_t dead = 0;
  int64_t fishes_vasculature = 0;
  int64_t fishes_reef = 0;
  float chemokines = 0;
  int64_t num_chemo_pts = 0;
  float floating_algaes = 0;
  float algae_on_substrate = 0;
  int64_t grazer_steps_total = 0;
  int64_t grazer_steps_on_coral_w_algae = 0;
  int64_t grazer_steps_on_coral_no_algae = 0;
  int64_t grazer_steps_on_sand_w_algae = 0;
  int64_t grazer_steps_on_sand_no_algae = 0;

  void init() {
    if (!rank_me()) {
      log_file.open(_options->output_dir + "/simforager.stats");
      log_file << "# time" << header(0) << endl;
    }
  }

  string header(int width) {
    vector<string> columns = {"incb", "expr", "apop", "dead",    "tvas",
                              "ttis", "chem", "algae", "chempts", "%infct"};
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
    //totals_d.push_back(reduce_one(floating_algaes, op_fast_add, 0).wait());  // / get_num_grid_points());
    totals_d.push_back(reduce_one(algae_on_substrate,  op_fast_add, 0).wait()); 
    auto all_chem_pts = reduce_one(num_chemo_pts, op_fast_add, 0).wait();
    totals_d.push_back(all_chem_pts + totals[0] + totals[1] + totals[2] + totals[3]);
    auto perc_infected =
        100.0 * (float)(totals[0] + totals[1] + totals[2] + totals[3]) / get_num_grid_points();

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

  void log(int time_step) {
    string s = to_str(0);
    if (!rank_me()) log_file << time_step << s << endl;
  }
};

ofstream _logstream;
bool _verbose = false;
SimStats _sim_stats;
shared_ptr<Options> _options;

IntermittentTimer generate_fish_timer(__FILENAME__ + string(":") + "generate fishes");
IntermittentTimer update_circulating_fishes_timer(__FILENAME__ + string(":") +
                                                  "update circulating fishes");
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





void seed_infection(Reef &reef, int time_step) {
  // pfft! This might be algae one day
  /*
  // _options->infection_coords contains the coords assigned just to rank_me()
  for (auto it = _options->infection_coords.begin(); it != _options->infection_coords.end(); it++) {
    auto infection_coords = *it;
    if (infection_coords[3] == time_step) {
      GridCoords coords({infection_coords[0], infection_coords[1], infection_coords[2]});
      auto coords_1d = coords.to_1d();
      int64_t num_tries = 0;
      while (true) {
        GridCoords new_coords(coords_1d);
        if (reef.set_initial_infection(coords_1d)) {
          // WARN("Time step ", time_step, ": SUCCESSFUL initial infection at ", new_coords.str() + " after ", num_tries, " tries");
          break;
        }
        num_tries++;
        coords_1d++;
        if (coords_1d >= get_num_grid_points()) {
          WARN("Could not find substrate to match uniform initial infection coord at ", coords.str());
          break;
        }
      }
      _options->infection_coords.erase(it--);
    }
  }
  */
  barrier();
  reef.add_new_actives(add_new_actives_timer);
  barrier();
}

// Generates the amount of fish specified at time step 0
void generate_fish(Reef &reef, int num_fish) {
    generate_fish_timer.start();

    int local_num = num_fish / rank_n();
    int rem = num_fish % rank_n();
    if (rank_me() < rem) local_num++;    

    if (rank_me() == 0) {
        SLOG("Generating fish: Total=", num_fish, ", Local=", local_num, "\n");
    }


    // Try 10 times to find an open location
    for ( int j = 0; j < num_fish; j++ ) {
      GridCoords coords(_rnd_gen);
      if (reef.try_add_new_reef_fish(coords.to_1d())) {
	_sim_stats.fishes_reef++;
	SLOG("Generated fish at ", coords.str(), "\n");
	
      }
    }
    
    // reef.change_num_circulating_fishes(num_fish);

    //#ifdef DEBUG
    // Validation: Ensure total fish matches the expected number
    //auto all_num = reduce_one(local_num, op_fast_add, 0).wait();
    //if (!rank_me() && all_num != num_fish) {
    //    DIE("Mismatch in generated fish: Total generated=", all_num, ", Expected=", num_fish);
    //}
    //#endif

    generate_fish_timer.stop();
}

void generate_fishes(Reef &reef, int time_step) {
  return;
  
  generate_fish_timer.start();
  int local_num = _options->fish_generation_rate / rank_n();
  int rem = _options->fish_generation_rate - local_num * rank_n();
  if (rank_me() < rem) local_num++;
  //if (time_step == 1) WARN("rem ", rem, " local num ", local_num, "\n");
  reef.change_num_circulating_fishes(local_num);
#ifdef DEBUG
  auto all_num = reduce_one(local_num, op_fast_add, 0).wait();
  if (!rank_me() && all_num != _options->fish_generation_rate)
    DIE("num fishes generated ", all_num, " != generation rate ", _options->fish_generation_rate);
#endif
  generate_fish_timer.stop();
}

int64_t get_rnd_coord(int64_t x, int64_t max_x) {
  int64_t new_x = x + _rnd_gen->get(0, 3) - 1;
  if (new_x < 0) new_x = 0;
  if (new_x >= max_x) new_x = max_x - 1;
  return new_x;
}

void update_circulating_fishes(int time_step, Reef &reef, double extravasate_fraction) {
  return;

  //update_circulating_fishes_timer.start();
  //auto num_circulating = reef.get_num_circulating_fishes();
  // fishes prob of dying in vasculature is 1/vascular_period
  //double portion_dying = (double)num_circulating / _options->fish_vascular_period;
  //int num_dying = floor(portion_dying);
  //if (_rnd_gen->trial_success(portion_dying - num_dying)) num_dying++;
  //reef.change_num_circulating_fishes(-num_dying);
  //_sim_stats.fishes_vasculature -= num_dying;
  //num_circulating = reef.get_num_circulating_fishes();
  //double portion_xtravasing = (time_step <= 15) ? extravasate_fraction * num_circulating : 0;
  //while (time_step < 20) double portion_xtravasing = ? extravasate_fraction * num_circulating : 0;
  //double portion_xtravasing = ? extravasate_fraction * num_circulating;
  //int num_xtravasing = floor(portion_xtravasing);
  //if (_rnd_gen->trial_success(portion_xtravasing - num_xtravasing)) num_xtravasing++;
  //for (int i = 0; i < num_xtravasing; i++) {
    progress();
    if (time_step != 0) return;

    // Try 10 times to find an open location
    for ( int i = 0; i < 10; i++ ) {
      GridCoords coords(_rnd_gen);
      if (reef.try_add_new_reef_fish(coords.to_1d())) {
	_sim_stats.fishes_reef++;
	SLOG(time_step, " fish extravasates at ", coords.str(), "\n");
	break;
      }
    }
    
    //_sim_stats.fishes_vasculature = num_circulating;
    update_circulating_fishes_timer.stop();
}

// --- Add near the top of main.cpp (or in a shared header) ---
static inline const char* substrate_color_name(SubstrateType t) {
  switch (t) {
    case SubstrateType::CORAL_WITH_ALGAE:  return "red";    // code 4
    case SubstrateType::CORAL_NO_ALGAE:    return "blue";   // code 2
    case SubstrateType::SAND_WITH_ALGAE:   return "green";  // code 3
    case SubstrateType::SAND_NO_ALGAE:     return "yellow"; // code 5
    case SubstrateType::NONE:
    default:                               return "black";  // code 1 / empty
  }
}

void update_reef_fish(int time_step, Reef &reef, GridPoint *grid_point, vector<int64_t> &nbs,
                         HASH_TABLE<int64_t, float> &chemokines_cache) {
  update_fish_timer.start();
  Fish *fish = grid_point->fish;

  // Log grazer position only once per timestep (before any early returns or movement)
  //if (fish->type == FishType::GRAZER)
     // log_grazer_step(fish->id, time_step, fish->x, fish->y, fish->z);

  //log grazer position, color, kappa

  if (fish->type == FishType::GRAZER) {
    const char* color = "black";
    if (grid_point->substrate) {
        color = substrate_color_name(grid_point->substrate->type);
    }

    // Determine the kappa value currently used (example shown; adapt as needed)
    double kappa_value = 0.0;
    if (fish->alert) {  // predator nearby
        switch (grid_point->substrate->type) {
            case SubstrateType::CORAL_WITH_ALGAE:  kappa_value = _options->kappa_grazer_w_predator_coral_w_algae; break;
            case SubstrateType::CORAL_NO_ALGAE:    kappa_value = _options->kappa_grazer_w_predator_coral_no_algae; break;
            case SubstrateType::SAND_WITH_ALGAE:   kappa_value = _options->kappa_grazer_w_predator_sand_w_algae; break;
            case SubstrateType::SAND_NO_ALGAE:     kappa_value = _options->kappa_grazer_w_predator_sand_no_algae; break;
            default:                               kappa_value = 0.0; break;
        }
    } else {  // no predator
        switch (grid_point->substrate->type) {
            case SubstrateType::CORAL_WITH_ALGAE:  kappa_value = _options->kappa_grazer_wo_predator_coral_w_algae; break;
            case SubstrateType::CORAL_NO_ALGAE:    kappa_value = _options->kappa_grazer_wo_predator_coral_no_algae; break;
            case SubstrateType::SAND_WITH_ALGAE:   kappa_value = _options->kappa_grazer_wo_predator_sand_w_algae; break;
            case SubstrateType::SAND_NO_ALGAE:     kappa_value = _options->kappa_grazer_wo_predator_sand_no_algae; break;
            default:                               kappa_value = 0.0; break;
        }
    }

    //log_grazer_step(fish->id, time_step, fish->x, fish->y, fish->z, , kappa_value);
}

  //count grazer time on substrate

  if (fish->type == FishType::GRAZER && grid_point->substrate){
     _sim_stats.grazer_steps_total++;
     switch (grid_point->substrate->type) {
          case SubstrateType::SAND_WITH_ALGAE: _sim_stats.grazer_steps_on_sand_w_algae ++; break;
          case SubstrateType::CORAL_WITH_ALGAE: _sim_stats.grazer_steps_on_coral_w_algae ++; break;
          case SubstrateType::CORAL_NO_ALGAE: _sim_stats.grazer_steps_on_coral_no_algae ++; break;
          case SubstrateType::SAND_NO_ALGAE: _sim_stats.grazer_steps_on_sand_no_algae ++; break;
          case SubstrateType::NONE: 
          defualt: break;
     }

     // Log grazer position once per timestep (before any early returns or movement)
     //log_grazer_step(fish->id, time_step, fish->x, fish->y, fish->z, static_cast<int>(grid_point->substrate->type), fish->kappa);

  }

  if (fish->moved) {
    // don't update fishes that were added this time step
    fish->moved = false;
    update_fish_timer.stop();
    return;
  }
  // Grazers consume algae on substrate on their current cell 
  if (fish->type == FishType::GRAZER && grid_point->substrate &&
      (grid_point->substrate->type == SubstrateType::CORAL_WITH_ALGAE || grid_point->substrate->type == SubstrateType::SAND_WITH_ALGAE)
      && grid_point->algae_on_substrate > 0) {

        float consumed = static_cast<float>(_options->algae_grazing_rate);
        float initial_algae = grid_point->algae_on_substrate;
        //grid_point->algae_on_substrate= std::max(0, grid_point->algae_on_substrate - consumed);

        if (initial_algae > consumed) 
          grid_point->algae_on_substrate = initial_algae - consumed;
        else 
          grid_point->algae_on_substrate = 0;

      

        //SLOG(_sim_stats.algae_on_substrate,"\n");
        //replacing with sand    
        //if (grid_point->algae_on_substrate <= 0 && _options->algae_turns_to_coral_when_depleted) 
          //grid_point->substrate->type = SubstrateType::CORAL;
    
        
        // debug log
        //SLOG("Grazer ", fish->id, " ate ", (initial_algae - grid_point->algae_on_substrate),
              //" algae at ", grid_point->coords.str(), " left=", grid_point->algae_on_substrate, "\n");
  }



  // Get count of neighbouring substrate types so fish don't bounce of regions of substrate by changing behaviour right on the edge
  float coral_w_algae_fraction = reef.count_neighbour_substrate(grid_point, SubstrateType::CORAL_WITH_ALGAE, 1, RadiusMetric::Chebyshev)/9.0;
  float coral_no_algae_fraction = reef.count_neighbour_substrate(grid_point, SubstrateType::CORAL_NO_ALGAE, 1, RadiusMetric::Chebyshev)/9.0;
  float sand_w_algae_fraction = reef.count_neighbour_substrate(grid_point, SubstrateType::SAND_WITH_ALGAE, 1, RadiusMetric::Chebyshev)/9.0;
  float sand_no_algae_fraction = reef.count_neighbour_substrate(grid_point, SubstrateType::SAND_NO_ALGAE, 1, RadiusMetric::Chebyshev)/9.0;
  
  // Sample von Mises with given kappa, mu = 0 
  double turning_angle = 0;
  if (fish->type == FishType::GRAZER){

    // Is there a predator nearby?
    int detection_radius = _options->grazer_detection_radius;
    fish->alert = reef.detect_neighbour_fish( grid_point, FishType::PREDATOR, detection_radius, RadiusMetric::Chebyshev);
    
    if ( fish->alert ){
      switch( grid_point->substrate->type ) {
      case SubstrateType::CORAL_WITH_ALGAE: turning_angle = coral_w_algae_fraction*sample_vonmises(0.0, _options->kappa_grazer_w_predator_coral_w_algae, gen);
	fish->kappa =  _options->kappa_grazer_w_predator_coral_w_algae;
        //SLOG("Type=Grazer Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	break;
      case SubstrateType::CORAL_NO_ALGAE:  turning_angle = coral_no_algae_fraction*sample_vonmises(0.0, _options->kappa_grazer_w_predator_coral_no_algae, gen);
	fish->kappa = _options->kappa_grazer_w_predator_coral_no_algae;
	//SLOG("Type=Grazer Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	break;
      case SubstrateType::SAND_WITH_ALGAE:  turning_angle = sand_w_algae_fraction*sample_vonmises(0.0, _options->kappa_grazer_w_predator_sand_w_algae, gen);
	fish->kappa = _options->kappa_grazer_w_predator_sand_w_algae;
	//SLOG("Type=Grazer Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	break;
      case SubstrateType::SAND_NO_ALGAE:  turning_angle = sand_no_algae_fraction*sample_vonmises(0.0, _options->kappa_grazer_w_predator_sand_no_algae, gen);
	fish->kappa = _options->kappa_grazer_w_predator_sand_no_algae;
  //SLOG("Type=Grazer Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
  break;
      default: WARN("In update_reef_fish unknown substrate type ", to_string(grid_point->substrate->type));
      }
    } else { // No predator nearby
      switch( grid_point->substrate->type ) {
      case SubstrateType::CORAL_WITH_ALGAE: turning_angle = coral_w_algae_fraction*sample_vonmises(0.0, _options->kappa_grazer_wo_predator_coral_w_algae, gen);
	fish->kappa =  _options->kappa_grazer_wo_predator_coral_w_algae;
	//SLOG("Type=Grazer Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	break;
      case SubstrateType::CORAL_NO_ALGAE:  turning_angle =  coral_no_algae_fraction*sample_vonmises(0.0, _options->kappa_grazer_wo_predator_coral_no_algae, gen);
	//SLOG("Type=Grazer Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	fish->kappa = _options->kappa_grazer_wo_predator_coral_no_algae;
	break;
      case SubstrateType::SAND_WITH_ALGAE: turning_angle = sand_w_algae_fraction*sample_vonmises(0.0, _options->kappa_grazer_wo_predator_sand_w_algae, gen);
	//SLOG("Type=Grazer Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	fish->kappa = _options->kappa_grazer_wo_predator_sand_w_algae;
	break;
      case SubstrateType::SAND_NO_ALGAE: turning_angle = sand_no_algae_fraction*sample_vonmises(0.0, _options->kappa_grazer_wo_predator_sand_no_algae, gen);
	fish->kappa = _options->kappa_grazer_wo_predator_sand_no_algae;
  //SLOG("Type=Grazer Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
  break;
      default: WARN("In update_reef_fish unknown substrate type ", to_string(grid_point->substrate->type));
      }
    } // End predator nearby conditional
  } else { // Fish is a predator
    // Is there a grazer nearby?
    int detection_radius = _options->predator_detection_radius;;
    fish->alert = reef.detect_neighbour_fish( grid_point, FishType::GRAZER, detection_radius, RadiusMetric::Chebyshev);
    if (fish->alert)
      {
	switch( grid_point->substrate->type ) {
	case SubstrateType::CORAL_WITH_ALGAE: turning_angle = coral_w_algae_fraction*sample_vonmises(0.0, _options->kappa_predator_w_grazer_coral_w_algae, gen);
	  //SLOG("Type=Predator Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	  fish->kappa = _options->kappa_predator_w_grazer_coral_w_algae;
	  break;
	case SubstrateType::CORAL_NO_ALGAE:  turning_angle =  coral_no_algae_fraction*sample_vonmises(0.0, _options->kappa_predator_w_grazer_coral_no_algae, gen);
	  //SLOG("Type=Predator Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	  fish->kappa = _options->kappa_predator_w_grazer_coral_no_algae;
	  break;
	case SubstrateType::SAND_WITH_ALGAE: turning_angle = sand_w_algae_fraction*sample_vonmises(0.0, _options->kappa_predator_w_grazer_sand_w_algae, gen);
	  //SLOG("Type=Predator Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	  fish->kappa = _options->kappa_predator_w_grazer_sand_w_algae;
	  break;
	case SubstrateType::SAND_NO_ALGAE: turning_angle = sand_no_algae_fraction*sample_vonmises(0.0, _options->kappa_predator_w_grazer_sand_no_algae, gen);
	  //SLOG("Type=Predator Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	  fish->kappa = _options->kappa_predator_w_grazer_sand_no_algae;
	  break;
	default: WARN("In update_reef_fish unknown substrate type ", to_string(grid_point->substrate->type));
	} // End of substrate conditional
      } else { // No prey nearby
      switch( grid_point->substrate->type ) {
      case SubstrateType::CORAL_WITH_ALGAE: turning_angle = coral_w_algae_fraction*sample_vonmises(0.0, _options->kappa_predator_wo_grazer_coral_w_algae, gen);
	//SLOG("Type=Predator Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	fish->kappa = _options->kappa_predator_wo_grazer_coral_w_algae;
	break;
      case SubstrateType::CORAL_NO_ALGAE:  turning_angle =  coral_no_algae_fraction*sample_vonmises(0.0, _options->kappa_predator_wo_grazer_coral_no_algae, gen);
	//SLOG("Type=Predator Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	fish->kappa = _options->kappa_predator_wo_grazer_coral_no_algae;
	break;
      case SubstrateType::SAND_WITH_ALGAE: turning_angle = sand_w_algae_fraction*sample_vonmises(0.0, _options->kappa_predator_wo_grazer_sand_w_algae, gen);
	//SLOG("Type=Predator Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	fish->kappa = _options->kappa_predator_wo_grazer_sand_w_algae;
	break;
      case SubstrateType::SAND_NO_ALGAE: turning_angle = sand_no_algae_fraction*sample_vonmises(0.0, _options->kappa_predator_wo_grazer_sand_no_algae, gen);
	//SLOG("Type=Predator Substrate=", to_string(grid_point->substrate->type), " Turning Angle = ", turning_angle,"\n");
	fish->kappa = _options->kappa_predator_wo_grazer_sand_no_algae;
	break;
      default: WARN("In update_reef_fish unknown substrate type ", to_string(grid_point->substrate->type));
      } // End of substrate conditional
    } // End of prey nearby conditional
  } // End of fish type conditional

  // Update fish orientation
fish->angle += turning_angle;

// Compute continuous displacement using step length and updated angle
auto [dx, dy] = polar_to_cartesian(fish->step_length, fish->angle, _grid_size);

// Compute continuous wrapped positions (toroidal world)
double new_xf = wrap_index(fish->x + dx, _grid_size->x);
double new_yf = wrap_index(fish->y + dy, _grid_size->y);
double new_zf = fish->z;  // assuming 2D motion for now

// Convert to nearest integer grid indices for cell mapping
auto [new_x, new_y] = nearest_grid_point(new_xf, new_yf, _grid_size);
int64_t new_z = static_cast<int64_t>(floor(new_zf + 0.5)) % _grid_size->z;

// Get 1D index for the target cell
int64_t selected_grid_i = GridCoords(new_x, new_y, new_z).to_1d();

// Attempt to move fish into the new location
for (int i = 0; i < 5; i++) {
    if (reef.try_add_reef_fish(selected_grid_i, *fish)) {

        // Update fish’s continuous position
        fish->x = new_xf;
        fish->y = new_yf;
        fish->z = new_zf;

        // Update occupancy bookkeeping
        delete grid_point->fish;
        grid_point->fish = nullptr;
        grid_point->substrate->status = SubstrateStatus::NO_FISH;

        // (Optional logging)
        // SLOG(time_step, "fish ", fish->id,
        //     " moves from ", grid_point->coords.str(),
        //     " to ", GridCoords(selected_grid_i).str(),
        //     " angle ", turning_angle, "\n");

        break;
    }
}

  
  update_fish_timer.stop();
}

void update_substrate(int time_step, Reef &reef, GridPoint *grid_point) {
  update_substrate_timer.start();
  if (grid_point->substrate->status == SubstrateStatus::DEAD)
    {
      update_substrate_timer.stop();
      return;
    }
  //std::srand(std::time(nullptr)); // Seed the random number generator

    if (std::rand() % 100 < 5) { // 5% chance to execute
      grid_point->substrate->infect();
      _sim_stats.incubating++;
      update_substrate_timer.stop();
      return;
    }
 
  /*
  if (!grid_point->substrate->infectable || grid_point->substrate->status == SubstrateStatus::DEAD) {
    update_substrate_timer.stop();
    return;
  }
  if (grid_point->substrate->status != SubstrateStatus::HEALTHY)
  DBG(time_step, " substrate ", grid_point->substrate->str(), "\n");
  bool produce_floating_algaes = false;
  switch (grid_point->substrate->status) {  
    case SubstrateStatus::HEALTHY: {
      
      double local_infectivity = _options->infectivity;
      if (grid_point->chemokine > 0) {
        local_infectivity *= _options->infectivity_multiplier;
      }
      if (grid_point->floating_algaes > 0) {
        if (_rnd_gen->trial_success(local_infectivity * grid_point->floating_algaes)) {
          grid_point->substrate->infect();
          _sim_stats.incubating++;
        }
      }
      break;
      }
    case SubstrateStatus::INCUBATING:
      if (grid_point->substrate->transition_to_expressing()) {
        _sim_stats.incubating--;
        _sim_stats.expressing++;
      }
      break;
      //case SubstrateStatus::EXPRESSING:
      case SubstrateStatus::HEALTHY:
	WARN("Substrate Died");
      if (grid_point->substrate->infection_death()) {
        _sim_stats.dead++;
        _sim_stats.expressing--;
      } else {
        produce_floating_algaes = true;
      }
      break;
    case SubstrateStatus::APOPTOTIC:
      if (grid_point->substrate->apoptosis_death()) {
        _sim_stats.dead++;
        _sim_stats.apoptotic--;
      } else if (grid_point->substrate->was_expressing()) {
        produce_floating_algaes = true;
      }
      break;
    default: break;
  }
  if (produce_floating_algaes) {
    double local_floating_algae_production = _options->floating_algae_production;
    if (grid_point->chemokine > 0) {
        local_floating_algae_production *= _options->floating_algae_production_multiplier;
    }
    grid_point->floating_algaes += local_floating_algae_production;
    grid_point->chemokine = min(grid_point->chemokine + _options->chemokine_production, 1.0);
  }
  update_substrate_timer.stop();
  */
}

void update_chemokines(GridPoint *grid_point, vector<int64_t> &nbs,
                       HASH_TABLE<int64_t, float> &chemokines_to_update) {
  update_concentration_timer.start();
  // Concentrations diffuse, i.e. the concentration at any single grid point tends to the average
  // of all the neighbors. So here we tell each neighbor what the current concentration is and
  // later those neighbors will compute their own averages. We do it in this "push" manner because
  // then we don't need to check the neighbors from every single grid point, but just push from
  // ones with concentrations > 0 (i.e. active grid points)
  if (grid_point->chemokine > 0) {
    grid_point->chemokine *= (1.0 - _options->chemokine_decay_rate);
    if (grid_point->chemokine < _options->min_chemokine) grid_point->chemokine = 0;
  }
  if (grid_point->chemokine > 0) {
    for (auto &nb_grid_i : nbs) {
      chemokines_to_update[nb_grid_i] += grid_point->chemokine;
    }
  }
  update_concentration_timer.stop();
}

void update_floating_algaes(GridPoint *grid_point, vector<int64_t> &nbs,
                    HASH_TABLE<int64_t, float> &floating_algaes_to_update) {
  update_concentration_timer.start();
  grid_point->floating_algaes = grid_point->floating_algaes * (1.0 - _options->floating_algae_clearance_rate);
  assert(grid_point->floating_algaes >= 0);
  if (grid_point->floating_algaes > 0) {
    for (auto &nb_grid_i : nbs) {
      floating_algaes_to_update[nb_grid_i] += grid_point->floating_algaes;
    }
  }
  update_concentration_timer.stop();
}

void diffuse(float &conc, float &nb_conc, double diffusion, int num_nbs) {
  // set to be average of neighbors plus self
  // amount that diffuses
  float conc_diffused = diffusion * conc;
  // average out diffused amount across all neighbors
  float conc_per_point = (conc_diffused + diffusion * nb_conc) / (num_nbs + 1);
  conc = conc - conc_diffused + conc_per_point;
  if (conc > 1.0) conc = 1.0;
  if (conc < 0) DIE("conc < 0: ", conc, " diffused ", conc_diffused, " pp ", conc_per_point);
  nb_conc = 0;
}

void spread_floating_algaes(float &floating_algaes, float &nb_floating_algaes, double diffusion, int num_nbs) {
  float floating_algaes_diffused = floating_algaes * diffusion;
  float floating_algaes_left = floating_algaes - floating_algaes_diffused;
  float avg_nb_floating_algaes = (floating_algaes_diffused + nb_floating_algaes * diffusion) / (num_nbs + 1);
  floating_algaes = floating_algaes_left + avg_nb_floating_algaes;
  nb_floating_algaes = 0;
}

void set_active_grid_points(Reef &reef) {
  set_active_points_timer.start();
  vector<GridPoint *> to_erase = {};
  // iterate through all active local grid points and set changes
   //_sim_stats.algae_on_substrate= 0;
  for (auto grid_point = reef.get_first_active_grid_point(); grid_point;
       grid_point = reef.get_next_active_grid_point()) {
    auto nbs = reef.get_neighbors(grid_point->coords);
    diffuse(grid_point->chemokine, grid_point->nb_chemokine, _options->chemokine_diffusion_coef,
            nbs->size());
    spread_floating_algaes(grid_point->floating_algaes, grid_point->nb_floating_algaes, _options->floating_algae_diffusion_coef,
                   nbs->size());
    if (grid_point->chemokine < _options->min_chemokine) grid_point->chemokine = 0;
    // only count up chemokine in healthy substrates or empty spaces
    // this will be added to the total number of infected and dead substrates to get cumulative
    // chemokine spread
    if (grid_point->chemokine > 0 &&
        (!grid_point->substrate || grid_point->substrate->status == SubstrateStatus::HEALTHY))
      _sim_stats.num_chemo_pts++;
    if (grid_point->floating_algaes > MAX_FLOATING_ALGAE) grid_point->floating_algaes = MAX_FLOATING_ALGAE;
    if (grid_point->floating_algaes < MIN_FLOATING_ALGAE) grid_point->floating_algaes = 0;
    if (grid_point->fish) grid_point->fish->moved = false;
    _sim_stats.chemokines += grid_point->chemokine;
    _sim_stats.floating_algaes += grid_point->floating_algaes;
    //_sim_stats.algae_on_substrate += grid_point->algae_on_substrate;

    if (!grid_point->is_active()) to_erase.push_back(grid_point);
  }
  for (auto grid_point : to_erase) reef.erase_active(grid_point);
  set_active_points_timer.stop();

  _sim_stats.algae_on_substrate = 0;
  for (auto gp = reef.get_first_local_grid_point(); gp; gp = reef.get_next_local_grid_point()) {
    _sim_stats.algae_on_substrate += gp->algae_on_substrate;
  }
}

void sample(int time_step,
            std::vector<SampleData> &samples,
            int64_t start_id,
            ViewObject view_object,
            Reef& reef)
{
  // Frame dimensions derived from grid & sampling
  int x_dim = _options->dimensions[0] / _options->sample_resolution;
  int y_dim = _options->dimensions[1] / _options->sample_resolution;
  int z_dim = _options->dimensions[2] / _options->sample_resolution;
  if (z_dim == 0) z_dim = 1;

  std::string video_path = std::filesystem::current_path().string() + "/reef.mp4";

  // Accumulators for one timestep
  static std::vector<std::tuple<int, int, cv::Scalar, float, int>> fish_points; // x,y,color,κ,thickness
  static std::vector<std::tuple<int, int, cv::Scalar>> coral_w_algae_points;
  static std::vector<std::tuple<int, int, cv::Scalar>> coral_no_algae_points;
  static std::vector<std::tuple<int, int, cv::Scalar>> sand_w_algae_points;
  static std::vector<std::tuple<int, int, cv::Scalar>> sand_no_algae_points;

  unsigned int n_fish = 0;

  for (int64_t i = 0; i < (int64_t)samples.size(); i++) {
    auto &sample = samples[i];
    int64_t index = start_id + i;
    auto [x, y, z] = GridCoords::to_3d(index);

    // FISH (white = grazer, yellow = predator)
    if (sample.fishes > 0) {
      n_fish++;
      cv::Scalar colour = (sample.fish_type == FishType::GRAZER)
                            ? cv::Scalar(255, 255, 255)    // Grazer white
                            : cv::Scalar(0, 0, 255); // Predator blue
      int thickness = (sample.fish_alert ? -1 : 2);
      float fish_kappa = sample.fish_kappa;

      // Keep (y, x) ordering to match your renderer
      fish_points.emplace_back(y, x, colour, fish_kappa, thickness);
    }

    if (sample.substrate_type == SubstrateType::CORAL_WITH_ALGAE) {
      cv::Scalar colour(0, 180, 165);
      
      coral_w_algae_points.emplace_back(y, x, colour);
    
    }
    if (sample.substrate_type == SubstrateType::CORAL_NO_ALGAE) {
      cv::Scalar colour(0, 0, 255);
      coral_no_algae_points.emplace_back(y, x, colour);
    }
    if (sample.substrate_type == SubstrateType::SAND_WITH_ALGAE) {
      cv::Scalar colour(0, 255, 0); 
      sand_w_algae_points.emplace_back(y, x, colour);
    }
    if (sample.substrate_type == SubstrateType::SAND_NO_ALGAE) {
      cv::Scalar colour(0, 255, 255);
      sand_no_algae_points.emplace_back(y, x, colour);
    }
      
    // After the final local sample, assemble and write the frame
    if (i == (int64_t)samples.size() - 1) {
      int width  = x_dim;
      int height = y_dim;

      write_full_frame_to_video(
        video_path, width, height,
        coral_w_algae_points,
        coral_no_algae_points,
        sand_w_algae_points,
        sand_no_algae_points,
        fish_points,
        /*scale=*/1
      ); // utils.cpp

      // Clear for next timestep
      fish_points.clear();
      coral_w_algae_points.clear();
      coral_no_algae_points.clear();
      sand_w_algae_points.clear();
      sand_no_algae_points.clear();

      SLOG("Rank ", upcxx::rank_me(),
           " at time step ", time_step,
           " wrote MP4 frame to: ", video_path, "\n");

      n_fish = 0;
    }
  }

  // --- BMP substrate export at the end of the whole run (keep this) ---
  if (view_object == ViewObject::SUBSTRATE
      && time_step == _options->num_timesteps - 1
      && upcxx::rank_me() == 0) {

    // Infer grid size from ecosystem_cells
    int max_x = 0, max_y = 0;
    for (int64_t id = 0; id < (int64_t)reef.get_ecosystem_cells().size(); ++id) {
      auto [gx, gy, gz] = GridCoords::to_3d(id);
      max_x = std::max(max_x, gx);
      max_y = std::max(max_y, gy);
    }
    int width  = max_x + 1;
    int height = max_y + 1;

    std::vector<std::vector<uint8_t>> bmp_array(height, std::vector<uint8_t>(width, 0));

    for (int64_t id = 0; id < (int64_t)reef.get_ecosystem_cells().size(); ++id) {
      const SubstrateType &type = reef.get_ecosystem_cells()[id];
      auto [gx, gy, gz] = GridCoords::to_3d(id);

      if (gx < 0 || gx >= width || gy < 0 || gy >= height) continue;

      uint8_t code = 0;
      switch (type) {
        case SubstrateType::SAND_NO_ALGAE:    code = 5; break;
        case SubstrateType::CORAL_WITH_ALGAE: code = 4; break;
        case SubstrateType::SAND_WITH_ALGAE:  code = 3; break;
        case SubstrateType::CORAL_NO_ALGAE:   code = 2; break;
        case SubstrateType::NONE:             code = 1; break;
      }
      bmp_array[gy][gx] = code;
    }

    std::string bmp_filename = "substrate.bmp";
    writeBMPColorMap(bmp_filename, bmp_array);
    SLOG("Rank ", upcxx::rank_me(), " wrote BMP input substrate to: ",
         std::filesystem::current_path().string(), "/", bmp_filename, "\n");
  }

  sample_write_timer.stop();
  upcxx::barrier();
}


// This function has each process iterate over the grid points points assigned to it and populates a sample datum. The nested for loops are so summary information can be created since a sample resolution greater than 1 need to combine the information for a region of cells into one sample.
int64_t get_samples(Reef &reef, vector<SampleData> &samples) {
  int64_t num_points =
      get_num_grid_points() / (_options->sample_resolution * _options->sample_resolution);
  if (_grid_size->z > 1) num_points /= _options->sample_resolution;
  int64_t num_points_per_rank = ceil((double)num_points / rank_n());
  int64_t start_id = rank_me() * num_points_per_rank;
  int64_t end_id = min((rank_me() + 1) * num_points_per_rank, num_points);
  samples.clear();
  if (end_id > start_id) {
    samples.reserve(end_id - start_id);
    int64_t i = 0;
    bool done = false;
    int block_size = _options->sample_resolution * _options->sample_resolution;
    if (_grid_size->z > 1) block_size *= _options->sample_resolution;
    vector<SampleData> block_samples;
    for (int x = 0; x < _grid_size->x && !done; x += _options->sample_resolution) {
      for (int y = 0; y < _grid_size->y && !done; y += _options->sample_resolution) {
        for (int z = 0; z < _grid_size->z; z += _options->sample_resolution) {
          if (i >= end_id) {
            done = true;
            break;
          }
          if (i >= start_id) {
            progress();
#ifdef AVERAGE_SUBSAMPLE
            float floating_algaes = 0;
            float chemokine = 0;
            int num_fishes = 0;
            bool substrate_found = false;
	    bool fish_found = false;
            array<int, 5> substrate_counts{0};
            block_samples.clear();
            bool done_sub = false;
            SubstrateType substrate_type = SubstrateType::NONE;
	    FishType fish_type = FishType::NONE;
	    bool fish_alert = false;
	    float fish_kappa = -1;
            for (int subx = x; subx < x + _options->sample_resolution; subx++) {
              if (subx >= _grid_size->x) break;
              for (int suby = y; suby < y + _options->sample_resolution; suby++) {
                if (suby >= _grid_size->y) break;
                for (int subz = z; subz < z + _options->sample_resolution; subz++) {
                  if (subz >= _grid_size->z) break;
                  auto sub_sd =
                      reef.get_grid_point_sample_data(GridCoords::to_1d(subx, suby, subz));
                  num_fishes += sub_sd.fishes;

		  if (sub_sd.has_fish)
		    {
		      fish_found = true;
		      switch (sub_sd.fish_type){
		      case FishType::NONE: //SLOG("Main.cpp:get_samples(): Fish type should not be NONE if the grid has a fish.\n");
			break;
		      case FishType::GRAZER:
			//SLOG("Main.cpp:get_samples(): Grazer.\n");
			fish_type = FishType::GRAZER;
			break;
		      case FishType::PREDATOR:
			//SLOG("Main.cpp:get_samples(): Predator.\n");
			fish_type = FishType::PREDATOR;
			break;
		      }

		      fish_kappa = sub_sd.fish_kappa;
		      fish_alert = sub_sd.fish_alert;
		      
		    }
		  
                  if (sub_sd.has_substrate) {
                    substrate_found = true;
                    /*
                    switch (sub_sd.substrate_status) {
                      case SubstrateStatus::HEALTHY: substrate_counts[0]++; break;
                      case SubstrateStatus::INCUBATING: substrate_counts[1]++; break;
                      case SubstrateStatus::EXPRESSING: substrate_counts[2]++; break;
                      case SubstrateStatus::APOPTOTIC: substrate_counts[3]++; break;
                      case SubstrateStatus::DEAD: substrate_counts[4]++; break;
                    }
                    */
                    switch (sub_sd.substrate_type) {
                      case SubstrateType::CORAL_WITH_ALGAE: substrate_type = SubstrateType::CORAL_WITH_ALGAE; break;
                      case SubstrateType::SAND_WITH_ALGAE: substrate_type = SubstrateType::SAND_WITH_ALGAE; break;
                      case SubstrateType::CORAL_NO_ALGAE: substrate_type = SubstrateType::CORAL_NO_ALGAE; break;
                      case SubstrateType::SAND_NO_ALGAE: substrate_type = SubstrateType::SAND_NO_ALGAE; break;
                      case SubstrateType::NONE: substrate_type = SubstrateType::NONE; break;
                    }
                  }
                  chemokine += sub_sd.chemokine;
                  floating_algaes += sub_sd.floating_algaes;

	     
                }
              }
            }
            
            SampleData sd = {.fishes = (double)num_fishes / block_size,
                             .has_substrate = substrate_found,
			     .has_fish = fish_found,
                             .substrate_type = substrate_type,
			     .fish_type = fish_type,
			     .fish_alert = fish_alert,
                             .floating_algaes = floating_algaes / block_size,
                             .chemokine = chemokine / block_size,
	                     .fish_kappa = fish_kappa};
#else
            auto sd = reef.get_grid_point_sample_data(GridCoords::to_1d(x, y, z));
#endif
            samples.push_back(sd);
          }
          i++;
        }
        if (done) break;
      }
      if (done) break;
    }
  }
  barrier();
  auto samples_written = reduce_one(samples.size(), op_fast_add, 0).wait();
  if (num_points != samples_written)
    SWARN("Number of point ", num_points, " != ", samples_written, " samples written");
  SLOG_VERBOSE("Number of samples written ", samples_written, "\n");
  return start_id;
}

void run_sim(Reef &reef) {
  BarrierTimer timer(__FILEFUNC__);

  SLOG("sample_period = ", _options->sample_period, "\n");
  
  auto start_t = NOW();
  auto curr_t = start_t;
  // TODO Allow for 1 timestep
  //auto five_perc = (_options->num_timesteps >= 50) ? _options->num_timesteps / 50 : 1;
  _sim_stats.init();
  int64_t ecosystem_volume = (int64_t)_options->ecosystem_dims[0] *
                              (int64_t)_options->ecosystem_dims[1] *
                              (int64_t)_options->ecosystem_dims[2];
  auto sim_volume = get_num_grid_points();
  double extravasate_fraction = (double)sim_volume / ecosystem_volume;
  SLOG("Fraction of circulating fishes extravasating is ", extravasate_fraction, "\n");
  SLOG("# datetime                    step    ", _sim_stats.header(STATS_COL_WIDTH),
       "<%active  lbln>\n");
  // store the total concentration increment updates for target grid points
  // chemokine, floating_algaes
  HASH_TABLE<int64_t, float> chemokines_to_update;
  HASH_TABLE<int64_t, float> chemokines_cache;
  HASH_TABLE<int64_t, float> floating_algaes_to_update;
  bool warned_boundary = false;
  vector<SampleData> samples;

  // Generate fish in sim
  generate_fish(reef, _options->num_fish);

  // Compute initial algae-on-substrate (before any grazing occurs)
  double local_algae_init = 0.0;
  int64_t local_green_on_coral = 0;
  int64_t local_green_on_sand = 0;


  //calculate for both algea
  for (auto gp = reef.get_first_local_grid_point(); gp; gp = reef.get_next_local_grid_point()) {
    local_algae_init += gp->algae_on_substrate;
    if (gp->substrate->type == SubstrateType::SAND_WITH_ALGAE) local_green_on_sand++;
    if (gp->substrate->type == SubstrateType::CORAL_WITH_ALGAE) local_green_on_coral++;
  }
 
  if (!rank_me()) {
    SLOG("Initial algae_on_substrate = ", local_algae_init, "\n");
  }

  SLOG("🚀 run_sim() has begun execution\n");
  
  for (int time_step = 0; time_step < _options->num_timesteps; time_step++) {
    //SLOG("Time step ", time_step, "\n");
    //SLOG("Grazer_timestep_total = ", _sim_stats.grazer_steps_total,"\n");
    //SLOG("Grazer_timestep_on_coral_w_algae = ", _sim_stats.grazer_steps_on_coral_w_algae,"\n");
    //SLOG("Grazer_timestep_on_coral_no_algae = ", _sim_stats.grazer_steps_on_coral_no_algae,"\n");
    //SLOG("Grazer_timestep_on_sand_w_algae = ", _sim_stats.grazer_steps_on_sand_w_algae,"\n");
    //SLOG("Grazer_timestep_on_sand_no_algae = ", _sim_stats.grazer_steps_on_sand_no_algae,"\n");

    seed_infection(reef, time_step);
    barrier();
    if (time_step == _options->antibody_period)
      _options->floating_algae_clearance_rate *= _options->antibody_factor;
    chemokines_to_update.clear();
    floating_algaes_to_update.clear();
    chemokines_cache.clear();
    if (time_step > _options->fish_initial_delay) {
      // generate_fishes(reef, time_step);
      barrier();
    }
    compute_updates_timer.start();
    update_circulating_fishes(time_step, reef, extravasate_fraction);
    // iterate through all active local grid points and update
    for (auto grid_point = reef.get_first_active_grid_point(); grid_point;
         grid_point = reef.get_next_active_grid_point()) {
      if (grid_point->chemokine > 0)
        DBG("chemokine\t", time_step, "\t", grid_point->coords.x, "\t", grid_point->coords.y, "\t",
            grid_point->coords.z, "\t", grid_point->chemokine, "\n");
      if (grid_point->floating_algaes > 0)
        DBG("floating_algaes\t", time_step, "\t", grid_point->coords.x, "\t", grid_point->coords.y, "\t",
            grid_point->coords.z, "\t", grid_point->floating_algaes, "\n");
      if (!warned_boundary && grid_point->substrate &&
          grid_point->substrate->status != SubstrateStatus::HEALTHY) {
        if (!grid_point->coords.x || grid_point->coords.x == _grid_size->x - 1 ||
            !grid_point->coords.y || grid_point->coords.y == _grid_size->y - 1 ||
            (_grid_size->z > 1 &&
             (!grid_point->coords.z || grid_point->coords.z == _grid_size->z - 1))) {
          //WARN("Hit boundary at ", grid_point->coords.str(), " ", grid_point->substrate->str(),
          //     " floating_algaes ", grid_point->floating_algaes, " chemokine ", grid_point->chemokine);
          warned_boundary = true;
        }
      }
      // DBG("updating grid point ", grid_point->str(), "\n");
      upcxx::progress();
      auto nbs = reef.get_neighbors(grid_point->coords);
      // the fishes are moved (added to the new list, but only cleared out at the end of all
      // updates)
      if (grid_point->fish)
        update_reef_fish(time_step, reef, grid_point, *nbs, chemokines_cache);
      //if (grid_point->substrate) update_substrate(time_step, reef, grid_point);
      //update_chemokines(grid_point, *nbs, chemokines_to_update);
      //update_floating_algaes(grid_point, *nbs, floating_algaes_to_update);
      if (grid_point->is_active()) reef.set_active(grid_point);
    }

    barrier();

    // Log fish locations
    for (auto grid_point = reef.get_first_active_grid_point(); grid_point; grid_point = reef.get_next_active_grid_point())
      {
	// Log grazer position once per timestep (before any early returns or movement)
	if ( grid_point->fish )
	  {
	    Fish* fish = grid_point->fish;
	    log_grazer_step(fish->id, time_step, fish->x, fish->y, fish->z, static_cast<int>(grid_point->substrate->type), fish->kappa);
	  }
      }
    
    barrier();
    compute_updates_timer.stop();
    //reef.accumulate_chemokines(chemokines_to_update, accumulate_concentrations_timer);
    //reef.accumulate_floating_algaes(floating_algaes_to_update, accumulate_concentrations_timer);
    barrier();
    //if (time_step % five_perc == 0 || time_step == _options->num_timesteps - 1) {
    if (time_step == _options->num_timesteps - 1) {
      auto num_actives = reduce_one(reef.get_num_actives(), op_fast_add, 0).wait();
      auto perc_actives = 100.0 * num_actives / get_num_grid_points();
      auto max_actives = reduce_one(reef.get_num_actives(), op_fast_max, 0).wait();
      auto load_balance = max_actives ? (double)num_actives / rank_n() / max_actives : 1;
      chrono::duration<double> t_elapsed = NOW() - curr_t;
      curr_t = NOW();
      SLOG("[", get_current_time(), " ", setprecision(2), fixed, setw(7), right, t_elapsed.count(),
           "s]: ", setw(8), left, time_step, _sim_stats.to_str(STATS_COL_WIDTH), setprecision(3),
           fixed, "< ", perc_actives, " ", load_balance, " >\n");
    }
    barrier();
    reef.add_new_actives(add_new_actives_timer);
    barrier();

    _sim_stats.floating_algaes = 0;
    _sim_stats.chemokines = 0;
    _sim_stats.num_chemo_pts = 0;
    _sim_stats.algae_on_substrate= 0;
    set_active_grid_points(reef);
    barrier();

    //if (_options->sample_period > 0 &&
    //    (time_step % _options->sample_period == 0 || time_step == _options->num_timesteps - 1)) {

     // SLOG("🔍 Checking sampling condition at timestep ", time_step, 
     //" (sample_period=", _options->sample_period, 
     //", num_timesteps=", _options->num_timesteps, ")\n");
      
      // Sample if the sample period is evenly dividible by the current time step count 
      if (time_step % _options->sample_period == 0)
	{
	  auto start = std::chrono::high_resolution_clock::now();
	  sample_timer.start();
	  samples.clear();
	  int64_t start_id = get_samples(reef, samples);

	  // Commenting out for now to avoid double mp4 frames
	  //sample(time_step, samples, start_id, ViewObject::SUBSTRATE, reef);


	  sample(time_step, samples, start_id, ViewObject::FISH, reef);
	  sample_timer.stop();
	  // End timer
	  auto end = std::chrono::high_resolution_clock::now();

	  // Compute elapsed time in seconds (double)
	  std::chrono::duration<double> elapsed = end - start;
	  
	  // Print with 2 decimals
	  std::cout << "Sampling took " << std::fixed << std::setprecision(2)
		    << elapsed.count() << " seconds\n";
	}

    log_timer.start();
    _sim_stats.log(time_step);
    barrier();
    log_timer.stop();

#ifdef DEBUG
    DBG("check actives ", time_step, "\n");
    reef.check_actives(time_step);
    barrier();
#endif
  }


  //code to count the final algea count
  double local_algae_final = 0.0;
  for (auto gp = reef.get_first_local_grid_point(); gp; gp = reef.get_next_local_grid_point()) {
    local_algae_final += gp->algae_on_substrate;
  }
  
  double pct = 0.0;  // default when init == 0
  if (local_algae_init > 0.0) {
      pct = (local_algae_init - local_algae_final) / local_algae_init * 100.0;
  }

  if (!rank_me()) SLOG("Total cells with Algae: ", local_green_on_sand + local_green_on_coral, "\n","On sand:", local_green_on_sand,"\n","On Coral", local_green_on_coral,"\n");
  if (!rank_me()) {
    SLOG(" Initial Algea Count = ", local_algae_init, "\n Final Algea Count = ", local_algae_final, " \n difference = ",  local_algae_init-local_algae_final,
      "\n reduction=", std::fixed, std::setprecision(2), pct, "%\n");
  }


  
  generate_fish_timer.done_all();
  //update_circulating_fishes_timer.done_all();
  update_fish_timer.done_all();
  //update_substrate_timer.done_all();
  //update_concentration_timer.done_all();
  compute_updates_timer.done_all();
  //accumulate_concentrations_timer.done_all();
  add_new_actives_timer.done_all();
  set_active_points_timer.done_all();
  sample_timer.done_all();
  sample_write_timer.done_all();
  log_timer.done_all();

  chrono::duration<double> t_elapsed = NOW() - start_t;
  SLOG("Finished ", _options->num_timesteps, " time steps in ", setprecision(4), fixed,
       t_elapsed.count(), " s (", (double)t_elapsed.count() / _options->num_timesteps,
       " s per step)\n");
  //printing the time distribution summary
  if (!rank_me()) {
    auto pct = [](int64_t a, int64_t t){ return (t > 0) ? (100.0 * double(a) / double(t)) : 0.0; };
    SLOG("\nGrazer time distribution (percent of grazer-steps over ", _sim_stats.grazer_steps_total, " (total time x grazer count) steps):\n",
         "  Coral with Algae Time: ", std::fixed, std::setprecision(2), pct(_sim_stats.grazer_steps_on_coral_w_algae, _sim_stats.grazer_steps_total), "%\n",
         "  Coral no Algae Time: ", std::fixed, std::setprecision(2), pct(_sim_stats.grazer_steps_on_coral_no_algae, _sim_stats.grazer_steps_total), "%\n",
         "  Sand with Algae Time : ", std::fixed, std::setprecision(2), pct(_sim_stats.grazer_steps_on_sand_w_algae, _sim_stats.grazer_steps_total), "%\n",
         "  Sand no Algae Time : ", std::fixed, std::setprecision(2), pct(_sim_stats.grazer_steps_on_sand_no_algae, _sim_stats.grazer_steps_total), "%\n");
  }
}

int main(int argc, char **argv) {

// Simulate random walk with Von Mises-distributed turning angles
  
  //gen.seed(static_cast<unsigned int>(314));
  gen.seed(static_cast<unsigned int>(std::time(nullptr))); // Use the current time to generate the seed
  
  upcxx::init();
  auto start_t = NOW();
  _options = make_shared<Options>();
  if (!_options->load(argc, argv)) return 0;
  
  
  ProgressBar::SHOW_PROGRESS = _options->show_progress;
  if (pin_thread(getpid(), local_team().rank_me()) == -1)
    WARN("Could not pin process ", getpid(), " to core ", rank_me());
  else
    SLOG_VERBOSE("Pinned processes, with process 0 (pid ", getpid(), ") pinned to core ",
                 local_team().rank_me(), "\n");
#ifdef BLOCK_PARTITION
  SLOG_VERBOSE("Using block partitioning\n");
#else
  SLOG_VERBOSE("Using linear partitioning\n");
#endif
#ifndef AVERAGE_SUBSAMPLE
  SLOG_VERBOSE("Not computing averages of subsamples\n");
#endif
  MemoryTrackerThread memory_tracker;
  memory_tracker.start();
  auto start_free_mem = get_free_mem();
  SLOG(KBLUE, "Starting with ", get_size_str(start_free_mem), " free on node 0", KNORM, "\n");
  Reef reef;
  SLOG(KBLUE, "Memory used on node 0 after initialization is  ",
       get_size_str(start_free_mem - get_free_mem()), KNORM, "\n");
  run_sim(reef);
  memory_tracker.stop();
  chrono::duration<double> t_elapsed = NOW() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
       " for SimForager version ", SIMFORAGER_VERSION, "\n");

  barrier();  // Ensure all ranks are done
  finalize_grazer_logs();
  if (rank_me() == 0) {
    finalize_video_writer();
  }
  barrier();  // Give other ranks time if needed

  upcxx::finalize();  // Only now shutdown UPCXX
  
  return 0;
}
