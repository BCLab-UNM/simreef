#pragma once

#include <math.h>
#include <stdarg.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/flat_aggr_store.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"
#include "options.hpp"

extern shared_ptr<Options> _options;

using upcxx::rank_me;
using upcxx::rank_n;

using std::array;
using std::list;
using std::pair;
using std::shared_ptr;
using std::to_string;
using std::vector;

enum class ViewObject { ALGAE, FISH, SUBSTRATE, CHEMOKINE };

inline string view_object_str(ViewObject view_object) {
  switch (view_object) {
    case ViewObject::ALGAE: return "algae";
    case ViewObject::FISH: return "fishreef";
    case ViewObject::SUBSTRATE: return "substrate";
    case ViewObject::CHEMOKINE: return "chemokine";
    default: DIE("Unknown view object");
  }
  return "";
}

struct GridBlocks {
  int64_t block_size, num_x, num_y, num_z, size_x, size_y, size_z;
};

inline GridBlocks _grid_blocks;

struct GridCoords;
inline shared_ptr<GridCoords> _grid_size = nullptr;

struct GridCoords {
  int x, y, z;

  GridCoords() {}

  GridCoords(int64_t x, int64_t y, int64_t z)
      : x(x)
      , y(y)
      , z(z) {}

  // create a grid point from 1d
  GridCoords(int64_t i);

  // create a random grid point
  GridCoords(shared_ptr<Random> rnd_gen);

  void set_rnd(shared_ptr<Random> rnd_gen);

  bool operator==(const GridCoords &coords) {
    return x == coords.x && y == coords.y && z == coords.z;
  }

  bool operator!=(const GridCoords &coords) {
    return x != coords.x || y != coords.y || z != coords.z;
  }

  int64_t to_1d() const;

  static int64_t to_1d(int x, int y, int z);

  static std::tuple<int, int, int> to_3d(int64_t i);
  
  // convert linear coord system to block - needed to use with external ecosystem model data
  static int64_t linear_to_block(int64_t i);

  string str() const {
    return "(" + to_string(x) + ", " + to_string(y) + ", " + to_string(z) + ")";
  }
};

struct Fish {
  string id;
  int binding_period = -1;
  int reef_time_steps = -1;
  bool moved = true;
  int x, y, z = -1;
  // turning angle used for CRW
  double angle = 0.0;
  UPCXX_SERIALIZED_FIELDS(id, binding_period, reef_time_steps, moved, x, y, z, angle);

  Fish(const string &id);

  Fish();
};

enum class SubstrateStatus { HEALTHY = 0, INCUBATING = 1, EXPRESSING = 2, APOPTOTIC = 3, DEAD = 4, NO_FISH = 5, FISH = 6};
const string SubstrateStatusStr[] = {"HEALTHY", "INCUBATING", "EXPRESSING", "APOPTOTIC", "DEAD", "NO_FISH", "FISH"};
enum class SubstrateType { CORAL, ALGAE, SAND, NONE };

inline std::string to_string(SubstrateType t) {
  switch (t) {
  case SubstrateType::CORAL: return "CORAL";
  case SubstrateType::ALGAE: return "ALGAE";
  case SubstrateType::SAND:  return "SAND";
   case SubstrateType::NONE:  return "NONE";
  }
  return "UNKNOWN";
}

class Substrate {
  int id;
  int incubation_time_steps = -1;
  int expressing_time_steps = -1;
  int apoptotic_time_steps = -1;

 public:
  SubstrateStatus status = SubstrateStatus::NO_FISH;
  SubstrateType type = SubstrateType::CORAL;
  bool infectable = true;

  Substrate(int id);

  string str();

  std::string str() const {
    std::ostringstream oss;
    oss 
      << "Substrate { "
      << "id="   << id
      << ", type="      << to_string(type)
      << "}";
      return oss.str();
  }

  void infect();
  bool transition_to_expressing();
  bool apoptosis_death();
  bool infection_death();
  bool is_active();
  double get_binding_prob();
  bool was_expressing();
};

// 3D size 4*3+26*8+16+8+16=372
struct GridPoint {
  GridCoords coords;
  // empty space is nullptr
  Substrate *substrate = nullptr;
  Fish *fish = nullptr;
  // starts off empty and if calculated because this grid point becomes active, it is saved
  vector<int64_t> *neighbors = nullptr;
  float chemokine = 0, nb_chemokine = 0;
  float floating_algaes = 0, nb_floating_algaes = 0;

  string str() const;

  bool is_active();
};

struct SampleData {
  double fishes = 0;
  bool has_substrate = false;
  SubstrateStatus substrate_status = SubstrateStatus::HEALTHY;
  SubstrateType substrate_type = SubstrateType::NONE;
  float floating_algaes = 0;
  float chemokine = 0;
};

inline int64_t get_num_grid_points() {
  return (int64_t)_grid_size->x * (int64_t)_grid_size->y * (int64_t)_grid_size->z;
}

class Reef {
 private:
  using grid_points_t = upcxx::dist_object<vector<GridPoint>>;
  grid_points_t grid_points;
  vector<GridPoint>::iterator grid_point_iter;

  // keeps track of all grid points that need to be updated
  using new_active_grid_points_t = upcxx::dist_object<HASH_TABLE<GridPoint *, bool>>;
  new_active_grid_points_t new_active_grid_points;

  HASH_TABLE<GridPoint *, bool> active_grid_points;
  HASH_TABLE<GridPoint *, bool>::iterator active_grid_point_iter;

  int64_t num_circulating_fishes;
  upcxx::dist_object<int64_t> fishes_generated;
  std::vector<SubstrateType> ecosystem_cells;

  // this is static for ease of use in rpcs
  static GridPoint *get_local_grid_point(grid_points_t &grid_points, int64_t grid_i);

  SubstrateType getSubstrateFromColor(uint8_t);

  std::vector<std::pair<int, SubstrateType>> load_bmp_cells();

  int load_bmp_file();

  int load_data_file(const string &fname, int num_grid_points, SubstrateType substrate_type);

  vector<int> get_model_dims(const string &fname);

 public:
  Reef();

  ~Reef() {}

  int64_t get_num_local_grid_points();

  intrank_t get_rank_for_grid_point(int64_t grid_i);

  vector<int64_t> *get_neighbors(GridCoords c);

  bool set_initial_infection(int64_t grid_i);

  void accumulate_chemokines(HASH_TABLE<int64_t, float> &chemokines_to_update,
                             IntermittentTimer &timer);

  void accumulate_floating_algaes(HASH_TABLE<int64_t, float> &floating_algaes_to_update, IntermittentTimer &timer);

  float get_chemokine(int64_t grid_i);

  bool fishes_in_neighborhood(GridPoint *grid_point);

  int64_t get_num_circulating_fishes();

  void change_num_circulating_fishes(int num);

  bool try_add_new_reef_fish(int64_t grid_i);

  bool try_add_reef_fish(int64_t grid_i, Fish &fish);

  SubstrateStatus try_bind_fish(int64_t grid_i);

  GridPoint *get_first_local_grid_point();
  GridPoint *get_next_local_grid_point();

  GridPoint *get_first_active_grid_point();
  GridPoint *get_next_active_grid_point();

  void set_active(GridPoint *grid_point);
  void erase_active(GridPoint *grid_point);

  void add_new_actives(IntermittentTimer &timer);

  size_t get_num_actives();

  SampleData get_grid_point_sample_data(int64_t grid_i);

  const std::vector<SubstrateType>& get_ecosystem_cells() const {
    return ecosystem_cells;
  }
  
  // int64_t get_random_airway_substrate_location();

#ifdef DEBUG
  void check_actives(int time_step);
#endif
};
