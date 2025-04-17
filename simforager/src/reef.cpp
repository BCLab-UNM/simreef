
// reef.cpp — Simulation Grid Implementation for SimForager
// --------------------------------------------------------
// Implements the core grid data structures and behavior including:
// - Coordinate conversion and partitioning (block/grid-based)
// - Substrate state transitions (e.g., infection, death)
// - Fish behavior and interaction
// - Initialization of grid points based on ecosystem type
// - RPC-safe accessors and parallelized updates

#include "reef.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
using upcxx::make_view;
using upcxx::progress;
using upcxx::rpc;
using upcxx::view;

// Converts a 1D index into x, y, z coordinates. Handles block-based partitions if enabled.
GridCoords::GridCoords(int64_t i) {
#ifdef BLOCK_PARTITION
  // If BLOCK_PARTITION is defined, decompose index into block and intra-block coordinates
  int64_t blocknum = i / _grid_blocks.block_size;
  int64_t block_z = blocknum / (_grid_blocks.num_x * _grid_blocks.num_y);
  blocknum -= block_z * _grid_blocks.num_x * _grid_blocks.num_y;
  int64_t block_y = blocknum / _grid_blocks.num_x;
  int64_t block_x = blocknum % _grid_blocks.num_x;
  int64_t in_block_id = i % _grid_blocks.block_size;
  block_x *= _grid_blocks.size_x;
  block_y *= _grid_blocks.size_y;
  block_z *= _grid_blocks.size_z;
  int64_t dz = in_block_id / (_grid_blocks.size_x * _grid_blocks.size_y);
  in_block_id -= (dz * _grid_blocks.size_x * _grid_blocks.size_y);
  int64_t dy = in_block_id / _grid_blocks.size_x;
  int64_t dx = in_block_id % _grid_blocks.size_x;
  x = block_x + dx;
  y = block_y + dy;
  z = block_z + dz;
#else
  // Simple 3D layout: convert from row-major linear index
  z = i / (_grid_size->x * _grid_size->y);
  i = i % (_grid_size->x * _grid_size->y);
  y = i / _grid_size->x;
  x = i % _grid_size->x;
#endif
}

// Constructs random coordinates using a shared random generator
GridCoords::GridCoords(shared_ptr<Random> rnd_gen) {
  x = rnd_gen->get(0, _grid_size->x);
  y = rnd_gen->get(0, _grid_size->y);
  z = rnd_gen->get(0, _grid_size->z);
}

// Converts 3D grid coordinates to a 1D linear index, considering block layout if enabled
int64_t GridCoords::to_1d(int x, int y, int z) {
  if (x >= _grid_size->x || y >= _grid_size->y || z >= _grid_size->z)
    DIE("Grid point is out of range: ", x, " ", y, " ", z, " max size ", _grid_size->str());
#ifdef BLOCK_PARTITION
  int64_t block_x = x / _grid_blocks.size_x;
  int64_t block_y = y / _grid_blocks.size_y;
  int64_t block_z = z / _grid_blocks.size_z;
  int64_t block_id =
      block_x + block_y * _grid_blocks.num_x + block_z * _grid_blocks.num_x * _grid_blocks.num_y;
  int64_t in_block_x = x % _grid_blocks.size_x;
  int64_t in_block_y = y % _grid_blocks.size_y;
  int64_t in_block_z = z % _grid_blocks.size_z;
  int64_t in_block_id = in_block_x + in_block_y * _grid_blocks.size_x +
                        in_block_z * _grid_blocks.size_x * _grid_blocks.size_y;
  return in_block_id + block_id * _grid_blocks.block_size;
#else
  return (int64_t)x + (int64_t)y * _grid_size->x + (int64_t)z * _grid_size->x * _grid_size->y;
#endif
}

// Overload for internal member use
int64_t GridCoords::to_1d() const {
  return GridCoords::to_1d(x, y, z);
}

// Used when input ecosystem files list 1D indices assuming non-blocked layout
int64_t GridCoords::linear_to_block(int64_t i) {
  int z = i / (_grid_size->x * _grid_size->y);
  i = i % (_grid_size->x * _grid_size->y);
  int y = i / _grid_size->x;
  int x = i % _grid_size->x;
  return GridCoords::to_1d(x, y, z);
}

// Assigns random coordinates using a generator (same as constructor)
void GridCoords::set_rnd(shared_ptr<Random> rnd_gen) {
  x = rnd_gen->get(0, _grid_size->x);
  y = rnd_gen->get(0, _grid_size->y);
  z = rnd_gen->get(0, _grid_size->z);
}

// Constructs a substrate and samples its state durations
Substrate::Substrate(int id)
    : id(id) {
  incubation_time_steps = _rnd_gen->get_poisson(_options->incubation_period);
  expressing_time_steps = _rnd_gen->get_poisson(_options->expressing_period);
  apoptotic_time_steps = _rnd_gen->get_poisson(_options->apoptosis_period);
  DBG("init substrate ", str(), "
");
}

// Format a string representation of this substrate (for debug)
string Substrate::str() {
  ostringstream oss;
  oss << id << " " << SubstrateStatusStr[(int)status] << " " << incubation_time_steps << " "
      << expressing_time_steps << " " << apoptotic_time_steps;
  return oss.str();
}

// Mark a substrate as newly infected
void Substrate::infect() {
  assert(status == SubstrateStatus::HEALTHY);
  assert(infectable);
  status = SubstrateStatus::INCUBATING;
}

// Progress incubation → expressing
bool Substrate::transition_to_expressing() {
  assert(status == SubstrateStatus::INCUBATING);
  incubation_time_steps--;
  if (incubation_time_steps > 0) return false;
  status = SubstrateStatus::EXPRESSING;
  return true;
}

// Used to determine whether the substrate was expressing before it became apoptotic
bool Substrate::was_expressing() {
  assert(status == SubstrateStatus::APOPTOTIC);
  return (incubation_time_steps == 0);
}

// Progress apoptosis toward death
bool Substrate::apoptosis_death() {
  assert(status == SubstrateStatus::APOPTOTIC);
  apoptotic_time_steps--;
  if (apoptotic_time_steps > 0) return false;
  status = SubstrateStatus::DEAD;
  return true;
}

// Progress expressing substrate toward death
bool Substrate::infection_death() {
  expressing_time_steps--;
  if (expressing_time_steps > 0) return false;
  status = SubstrateStatus::DEAD;
  return true;
}

// Is this substrate currently in an active state?
bool Substrate::is_active() {
  return (status != SubstrateStatus::HEALTHY && status != SubstrateStatus::DEAD);
}

// Calculates the binding probability for a fish, based on infection state
double Substrate::get_binding_prob() {
  if (status == SubstrateStatus::EXPRESSING || status == SubstrateStatus::APOPTOTIC)
    return _options->max_binding_prob;
  double scaling = 1.0 - (double)incubation_time_steps / _options->incubation_period;
  if (scaling < 0) scaling = 0;
  double prob = _options->max_binding_prob * scaling;
  return min(prob, _options->max_binding_prob);
}

// More functionality defined in this file, but this covers core data structure setup
