
// reef.hpp — Core simulation grid structures for SimForager
// ----------------------------------------------------------
// Defines GridCoords, GridPoint, Substrate, Fish, and the main Reef class.
// These are the primary data structures used to represent the spatial layout,
// biological states, and agent-based entities in the simulation.
//
// Grid is partitioned across processes using either a flat grid or a block-based
// decomposition (when BLOCK_PARTITION is defined). All 3D coordinates can be
// mapped to linear indices and vice versa, enabling efficient parallel access.

#ifndef REEF_HPP
#define REEF_HPP

#include <memory>
#include <vector>
#include <string>
#include <set>
#include <unordered_map>


// Restored required typedefs and extern variables
#include <unordered_map>
#include <set>
#include "bytell_hash_map.hpp"
#include "upcxx/upcxx.hpp"
#include "upcxx_utils.hpp"

template <typename K, typename V>
using HASH_TABLE = ankerl::unordered_dense::map<K, V>;

using grid_points_t = std::unordered_map<int64_t, GridPoint *>;
using new_active_grid_points_t = std::set<GridPoint *>;

extern std::shared_ptr<struct GridSize> _grid_size;
extern struct GridBlocks _grid_blocks;
extern std::shared_ptr<class Random> _rnd_gen;

#include "options.hpp"
#include "utils.hpp"

using std::string;
using std::shared_ptr;
using std::vector;

// Enum for substrate status, tracking progression of infection
enum class SubstrateStatus {
    HEALTHY,
    INCUBATING,
    EXPRESSING,
    APOPTOTIC,
    DEAD
};

// Enum for substrate type (e.g. alveoli, airway, coral etc.)
enum class SubstrateType {
    NONE,
    ALVEOLI,
    AIRWAY,
    CORAL
};

// Forward declaration
struct GridPoint;

// Main grid coordinate class
struct GridCoords {
    int64_t x, y, z;

    GridCoords() : x(0), y(0), z(0) {}
    GridCoords(int64_t x_, int64_t y_, int64_t z_) : x(x_), y(y_), z(z_) {}
    GridCoords(int64_t linear_index); // converts from 1D to 3D (defined in reef.cpp)
    GridCoords(shared_ptr<class Random> rnd_gen);  // randomly initialized

    static int64_t to_1d(int x, int y, int z); // 3D to 1D conversion
    int64_t to_1d() const;                     // uses internal x,y,z
    static int64_t linear_to_block(int64_t);   // map linear index if using blocks
    string str() const {
        return "(" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + ")";
    }

    void set_rnd(shared_ptr<class Random> rnd_gen);
};

// Stores the state of a biological substrate
struct Substrate {
    int id;
    SubstrateStatus status = SubstrateStatus::HEALTHY;
    SubstrateType type = SubstrateType::NONE;
    int incubation_time_steps = 0;
    int expressing_time_steps = 0;
    int apoptotic_time_steps = 0;
    bool infectable = false;

    Substrate(int id);

    string str();
    void infect();
    bool transition_to_expressing();
    bool was_expressing();  // returns true if expressing before apoptosis
    bool apoptosis_death(); // progresses to DEAD
    bool infection_death(); // progresses to DEAD
    bool is_active();       // not healthy or dead
    double get_binding_prob();
};

// Stores fish (mobile agents) that swim through the grid
struct Fish {
    string id;
    double angle = 0.0;       // current heading angle
    int reef_time_steps = 0;  // life duration
    int x = 0, y = 0, z = 0;  // grid location
    bool moved = false;

    Fish();
    Fish(const string &id);
};

// Each point in the grid stores local substrate, fish, and concentrations
struct GridPoint {
    GridCoords coords;
    Substrate *substrate = nullptr;
    Fish *fish = nullptr;
    float chemokine = 0;
    float nb_chemokine = 0;
    float floating_algaes = 0;
    float nb_floating_algaes = 0;
    vector<int64_t> *neighbors = nullptr;

    string str() const;
    bool is_active(); // whether this point has any active biological agents
};

// Stores sampling data for VTK output
struct SampleData {
    bool has_substrate = false;
    SubstrateStatus substrate_status;
    int fishes = 0;
    float chemokine = 0;
    float floating_algaes = 0;
};

// Holds block dimensions for block-decomposed grid
struct GridBlocks {
    int size_x, size_y, size_z;
    int num_x, num_y, num_z;
    int64_t block_size;
};

// Main simulation class
class Reef {
public:
    Reef();

    int64_t get_num_local_grid_points();
    int get_num_circulating_fishes() const { return num_circulating_fishes; }
    void change_num_circulating_fishes(int diff) { num_circulating_fishes += diff; }

    SampleData get_grid_point_sample_data(int64_t grid_i);
    static vector<int64_t> *get_neighbors(GridCoords c);
    GridPoint *get_first_active_grid_point();
    GridPoint *get_next_active_grid_point();
    void erase_active(GridPoint *grid_point);
    void add_new_actives(IntermittentTimer &timer);

    bool set_initial_infection(int64_t grid_i);
    float get_chemokine(int64_t grid_i);

    void accumulate_chemokines(HASH_TABLE<int64_t, float> &chemokines_to_update,
                               IntermittentTimer &timer);
    void accumulate_floating_algaes(HASH_TABLE<int64_t, float> &floating_algaes_to_update,
                                    IntermittentTimer &timer);

private:
    int load_data_file(const string &fname, int num_grid_points, SubstrateType substrate_type);
    intrank_t get_rank_for_grid_point(int64_t grid_i);
    static GridPoint *get_local_grid_point(grid_points_t &grid_points, int64_t grid_i);

    grid_points_t grid_points;
    new_active_grid_points_t new_active_grid_points;
    vector<SubstrateType> ecosystem_cells;
    int num_circulating_fishes;
    vector<int64_t> fishes_generated;
};

#endif
