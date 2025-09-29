#include "reef.hpp"
#include <filesystem>

using namespace std;
using upcxx::make_view;
using upcxx::progress;
using upcxx::rpc;
using upcxx::view;

// Random number generator
boost::random::mt19937 vonmises_gen;

std::tuple<int, int, int> GridCoords::to_3d(int64_t i) {
#ifdef BLOCK_PARTITION
  int64_t blocknum = i / _grid_blocks.block_size;

  int64_t block_z = blocknum / (_grid_blocks.num_x * _grid_blocks.num_y);
  blocknum -= block_z * _grid_blocks.num_x * _grid_blocks.num_y;

  int64_t block_y = blocknum / _grid_blocks.num_x;
  int64_t block_x = blocknum % _grid_blocks.num_x;

  int64_t in_block_id = i % _grid_blocks.block_size;
  int64_t dz = in_block_id / (_grid_blocks.size_x * _grid_blocks.size_y);
  in_block_id -= dz * _grid_blocks.size_x * _grid_blocks.size_y;

  int64_t dy = in_block_id / _grid_blocks.size_x;
  int64_t dx = in_block_id % _grid_blocks.size_x;

  int x = block_x * _grid_blocks.size_x + dx;
  int y = block_y * _grid_blocks.size_y + dy;
  int z = block_z * _grid_blocks.size_z + dz;

  return {x, y, z};
#else
  int z = i / (_grid_size->x * _grid_size->y);
  i = i % (_grid_size->x * _grid_size->y);
  int y = i / _grid_size->x;
  int x = i % _grid_size->x;
  return {x, y, z};
#endif
}

GridCoords::GridCoords(int64_t i) {
#ifdef BLOCK_PARTITION
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
  z = i / (_grid_size->x * _grid_size->y);
  i = i % (_grid_size->x * _grid_size->y);
  y = i / _grid_size->x;
  x = i % _grid_size->x;
#endif
}

GridCoords::GridCoords(shared_ptr<Random> rnd_gen) {
  x = rnd_gen->get(0, _grid_size->x);
  y = rnd_gen->get(0, _grid_size->y);
  z = rnd_gen->get(0, _grid_size->z);
}

int64_t GridCoords::to_1d(int x, int y, int z) {
  if (x >= _grid_size->x || y >= _grid_size->y || z >= _grid_size->z) {
    DIE("Grid point is out of range: ", x, " ", y, " ", z, " max size ", _grid_size->str());
  }
  
  if (x >= _grid_size->x) x = _grid_size->x-1;
  if (y >= _grid_size->y) y = _grid_size->y-1;
  if (z >= _grid_size->z) z = _grid_size->z-1;
  
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

int64_t GridCoords::linear_to_block(int64_t i) {
  int z = i / (_grid_size->x * _grid_size->y);
  i = i % (_grid_size->x * _grid_size->y);
  int y = i / _grid_size->x;
  int x = i % _grid_size->x;
  return GridCoords::to_1d(x, y, z);
}

int64_t GridCoords::to_1d() const { return GridCoords::to_1d(x, y, z); }

void GridCoords::set_rnd(shared_ptr<Random> rnd_gen) {
  x = rnd_gen->get(0, _grid_size->x);
  y = rnd_gen->get(0, _grid_size->y);
  z = rnd_gen->get(0, _grid_size->z);
}

Fish::Fish(const string &id)
    : id(id) {
  reef_time_steps = _rnd_gen->get_poisson(_options->fish_reef_period);
  DBG("init fish ", id, " ", reef_time_steps, "\n");
}

Fish::Fish() { reef_time_steps = _rnd_gen->get_poisson(_options->fish_reef_period); }

Substrate::Substrate(int id)
    : id(id) {
  incubation_time_steps = _rnd_gen->get_poisson(_options->incubation_period);
  expressing_time_steps = _rnd_gen->get_poisson(_options->expressing_period);
  apoptotic_time_steps = _rnd_gen->get_poisson(_options->apoptosis_period);
  DBG("init substrate ", str(), "\n");
}

string Substrate::str() {
  ostringstream oss;
  oss << id << " " << SubstrateStatusStr[(int)status] << " " << incubation_time_steps << " "
      << expressing_time_steps << " " << apoptotic_time_steps;
  return oss.str();
}

void Substrate::infect() {
  return; // We dont care about this anymore
   assert(status == SubstrateStatus::HEALTHY);
  assert(infectable);
  status = SubstrateStatus::INCUBATING;
}

bool Substrate::transition_to_expressing() {
  assert(status == SubstrateStatus::INCUBATING);
  incubation_time_steps--;
  if (incubation_time_steps > 0) return false;
  status = SubstrateStatus::EXPRESSING;
  return true;
}

bool Substrate::was_expressing() {
  // this is used to determine if the substrate was expressing before apoptosis was induced
  assert(status == SubstrateStatus::APOPTOTIC);
  return (incubation_time_steps == 0);
}

bool Substrate::apoptosis_death() {
  assert(status == SubstrateStatus::APOPTOTIC);
  apoptotic_time_steps--;
  if (apoptotic_time_steps > 0) return false;
  status = SubstrateStatus::DEAD;
  return true;
}

bool Substrate::infection_death() {
  expressing_time_steps--;
  if (expressing_time_steps > 0) return false;
  status = SubstrateStatus::DEAD;
  return true;
}

bool Substrate::is_active() {
  return status != SubstrateStatus::NO_FISH;
}

double Substrate::get_binding_prob() {
  // binding prob is linearly scaled from 0 to 1 for incubating cells over the course of the
  // incubation period, but is always 1 for expressing cells
  if (status == SubstrateStatus::EXPRESSING || status == SubstrateStatus::APOPTOTIC)
    return _options->max_binding_prob;
  double scaling = 1.0 - (double)incubation_time_steps / _options->incubation_period;
  if (scaling < 0) scaling = 0;
  double prob = _options->max_binding_prob * scaling;
  return min(prob, _options->max_binding_prob);
}

string GridPoint::str() const {
  ostringstream oss;
  oss << "xyz " << coords.str() << ", substrate " << (substrate ? substrate->str() : "none") << ", v "
      << floating_algaes << ", c " << chemokine;
  return oss.str();
}

bool GridPoint::is_active() {
  // it could be incubating but without anything else set
  return ((substrate && substrate->is_active()) || floating_algaes > 0 || chemokine > 0 || fish);
}

static int get_cube_block_dim(int64_t num_grid_points) {
  int block_dim = 1;
  for (int d = 1; d < _options->max_block_dim; d++) {
    // only allow dims that divide each main dimension perfectly
    if (remainder(_grid_size->x, d) || remainder(_grid_size->y, d) || remainder(_grid_size->z, d)) {
      DBG("dim ", d, " does not divide all main dimensions cleanly\n");
      continue;
    }
    size_t cube = (size_t)pow((double)d, 3.0);
    size_t num_cubes = num_grid_points / cube;
    DBG("cube size ", cube, " num cubes ", num_cubes, "\n");
    if (num_cubes < rank_n() * MIN_BLOCKS_PER_PROC) {
      DBG("not enough cubes ", num_cubes, " < ", rank_n() * MIN_BLOCKS_PER_PROC, "\n");
      break;
    }
    // there is a remainder - this is not a perfect division
    if (remainder(num_grid_points, cube)) {
      DBG("there is a remainder - don't use\n");
      continue;
    }
    // skip sizes that distribute the blocks in columns
    if (d > 1 && (_grid_size->x % (d * rank_n()) == 0 || _grid_size->y % (d * rank_n()) == 0)) {
      DBG("dim ", d, " gives perfect division of all blocks into x axis - skip\n");
      continue;
    }
    DBG("selected dim ", d, "\n");
    block_dim = d;
  }
  return block_dim;
}

static int get_square_block_dim(int64_t num_grid_points) {
  int block_dim = 1;
  for (int d = 1; d < _options->max_block_dim; d++) {
    // only allow dims that divide each main dimension perfectly
    if (remainder(_grid_size->x, d) || remainder(_grid_size->y, d)) {
      DBG("dim ", d, " does not divide all main dimensions cleanly\n");
      continue;
    }
    size_t square = (size_t)pow((double)d, 2.0);
    size_t num_squares = num_grid_points / square;
    DBG("square size ", square, " num squares ", num_squares, "\n");
    if (num_squares < rank_n() * MIN_BLOCKS_PER_PROC) {
      DBG("not enough squares ", num_squares, " < ", rank_n() * MIN_BLOCKS_PER_PROC, "\n");
      break;
    }
    // there is a remainder - this is not a perfect division
    if (remainder(num_grid_points, square)) {
      DBG("there is a remainder - don't use\n");
      continue;
    }
    // skip sizes that distribute the blocks in columns
    if (d > 1 && _grid_size->x % (d * rank_n()) == 0) {
      DBG("dim ", d, " gives perfect division of all blocks into x axis - skip\n");
      continue;
    }
    DBG("selected dim ", d, "\n");
    block_dim = d;
  }
  return block_dim;
}

Reef::Reef()
    : grid_points({})
    , new_active_grid_points({})
    , num_circulating_fishes(0)
    , fishes_generated({0}) {
  auto remainder = [](int64_t numerator, int64_t denominator) -> bool {
    return ((double)numerator / denominator - (numerator / denominator) != 0);
  };
  BarrierTimer timer(__FILEFUNC__, false, true);

  _grid_size = make_shared<GridCoords>(
      GridCoords(_options->dimensions[0], _options->dimensions[1], _options->dimensions[2]));
  int64_t num_grid_points = get_num_grid_points();
  // find the biggest cube that perfectly divides the grid and gives enough
  // data for at least two cubes per rank (for load
  // balance)
  // This is a trade-off: the more data is blocked, the better the locality,
  // but load balance could be a problem if not all
  // ranks get the same number of cubes. Also, having bigger cubes could lead
  // to load imbalance if all of the computation is
  // happening within a cube.
  int64_t block_dim = (_grid_size->z > 1 ? get_cube_block_dim(num_grid_points) :
                                           get_square_block_dim(num_grid_points));
  if (block_dim == 1)
    SWARN("Using a block size of 1: this will result in a lot of "
          "communication. You should change the dimensions.");

  _grid_blocks.block_size =
      (_grid_size->z > 1 ? block_dim * block_dim * block_dim : block_dim * block_dim);
  _grid_blocks.num_x = _grid_size->x / block_dim;
  _grid_blocks.num_y = _grid_size->y / block_dim;
  _grid_blocks.num_z = (_grid_size->z > 1 ? _grid_size->z / block_dim : 1);
  _grid_blocks.size_x = block_dim;
  _grid_blocks.size_y = block_dim;
  _grid_blocks.size_z = (_grid_size->z > 1 ? block_dim : 1);

  int64_t num_blocks = num_grid_points / _grid_blocks.block_size;

  int64_t blocks_per_rank = ceil((double)num_blocks / rank_n());

  bool threeD = _grid_size->z > 1;
  SLOG("Dividing ", num_grid_points, " grid points into ", num_blocks,
       (threeD ? " blocks" : " squares"), " of size ", _grid_blocks.block_size, " (", block_dim,
       "^", (threeD ? 3 : 2), "), with ", blocks_per_rank, " per process\n");
  double sz_grid_point = sizeof(GridPoint) + (double)sizeof(Substrate);
  auto mem_reqd = sz_grid_point * blocks_per_rank * _grid_blocks.block_size;
  SLOG("Total initial memory required per process is at least ", get_size_str(mem_reqd),
       " with each grid point requiring on average ", sz_grid_point, " bytes\n");
  int64_t num_ecosystem_cells = 0;
  if (!_options->ecosystem_model_dir.empty()) {
    ecosystem_cells.resize(num_grid_points, SubstrateType::NONE);
    Timer t_load_ecosystem_model("load ecosystem model");
    t_load_ecosystem_model.start();
    num_ecosystem_cells += load_data_file(_options->ecosystem_model_dir + "/bronchiole.dat", num_grid_points,
                                     SubstrateType::CORAL_WITH_ALGAE);
    t_load_ecosystem_model.stop();
    SLOG("Lung model loaded ", num_ecosystem_cells, " substrate in ", fixed, setprecision(2),
         t_load_ecosystem_model.get_elapsed(), " s\n");
  }
  // Load BMP file to substrate types
  // ecosystem_cells.resize(num_grid_points, SubstrateType::NONE);
  Timer t_bmp("load BMP substrate map");
  t_bmp.start();
  num_ecosystem_cells += load_bmp_file(); // += hmm?
  t_bmp.stop();
  SLOG("BMP substrate map loaded in ", t_bmp.get_elapsed(), " s \n");
  // FIXME: these blocks need to be stride distributed to better load balance
  grid_points->reserve(blocks_per_rank * _grid_blocks.block_size);
  for (int64_t i = 0; i < blocks_per_rank; i++) {
    int64_t start_id = (i * rank_n() + rank_me()) * _grid_blocks.block_size;
    if (start_id >= num_grid_points) break;
    for (auto id = start_id; id < start_id + _grid_blocks.block_size; id++) {
      assert(id < num_grid_points);
      GridCoords coords(id);
      if (num_ecosystem_cells) {
        if (ecosystem_cells[id] != SubstrateType::NONE) {
          Substrate *substrate = new Substrate(id);
          substrate->type = ecosystem_cells[id];
          substrate->infectable = true;
          //seeding the initial algea count
          GridPoint gp{coords, substrate};

          if (substrate->type == SubstrateType::CORAL_WITH_ALGAE || substrate->type == SubstrateType::SAND_WITH_ALGAE ) {
              gp.algae_on_substrate = static_cast<float>(_options->algae_init_count);
          } else {
              gp.algae_on_substrate = 0.0;
            } 

          //grid_points->emplace_back(GridPoint({coords, substrate}));
          grid_points->emplace_back(std::move(gp));

        } else {  // Add empty space == air
          //grid_points->emplace_back(GridPoint({coords, nullptr}));
          GridPoint gp{coords, nullptr};
          gp.algae_on_substrate = 0.0;  // nothing attached in empty space
          grid_points->emplace_back(std::move(gp));
        }
      } else {
        Substrate *substrate = new Substrate(id);
        substrate->type = SubstrateType::CORAL_WITH_ALGAE;
        // substrate->status = static_cast<SubstrateStatus>(rank_me() % 4);
        substrate->infectable = true;
        grid_points->emplace_back(GridPoint({coords, substrate}));
      }
#ifdef DEBUG
      DBG("adding grid point ", id, " at ", coords.str(), "\n");
      auto id_1d = coords.to_1d();
      if (id_1d != id) DIE("id ", id, " is not same as returned by to_1d ", id_1d);
      auto nbs = get_neighbors(coords);
      ostringstream oss;
      for (auto nb_grid_i : *nbs) {
        oss << GridCoords(nb_grid_i).str() << " ";
      }
      DBG("nbs: ", oss.str(), "\n");
#endif
    }
  }
  barrier();
}

SubstrateType Reef::getSubstrateFromColor(uint8_t value) {
  switch (value) {
  case 1: return SubstrateType::NONE;
  case 2: return SubstrateType::CORAL_NO_ALGAE;
  case 3: return SubstrateType::SAND_WITH_ALGAE;
  case 4: return SubstrateType::CORAL_WITH_ALGAE;
  case 5: return SubstrateType::SAND_NO_ALGAE;
  default: SDIE("Reef:getSubstrateFromColor() unknown type.");
  }
  
  return SubstrateType::NONE;
}

std::vector<std::pair<int, SubstrateType>> Reef::load_bmp_cells() {
  auto pixels = readBMPColorMap(_options->substrate_bitmap_path);
  int height = (int)pixels.size();
  int width  = height > 0 ? (int)pixels[0].size() : 0;
  int depth  = _grid_size->z;
  std::vector<std::pair<int, SubstrateType>> cells;
  cells.reserve((size_t)height * width * depth);

  // BMP rows are bottom-up
  for (int row = 0; row < height; ++row) {
	  for (int col = 0; col < width; ++col) {
		  for (int z = 0; z < depth; ++z)
		  {
			  int id = GridCoords::to_1d(row, col, z);
			  id = GridCoords::linear_to_block(id);
			  cells.emplace_back(id, SubstrateType::NONE);
			  }
		  }
	  }
for (int row = 0; row < height; ++row) {
	for (int col = 0; col < width; ++col) {
		if (row == col) {
			int id = GridCoords::to_1d(row, col, 0);
                          id = GridCoords::linear_to_block(id);
                          cells.emplace_back(id, SubstrateType::CORAL_WITH_ALGAE);
		}
	}
}

  /*
  for (int row = 0; row < height; ++row) {
      int flipped_y = height - 1 - row;
      for (int col = 0; col < width; ++col) {
          uint8_t color = pixels[row][col];
          if (color == 0) continue;  // skip
          SubstrateType st = getSubstrateFromColor(color);
	  int id = GridCoords::to_1d(col, flipped_y,
          // replicate this 2D map across all z layers
	  
          for (int z = 0; z < depth; ++z) {
              int id = GridCoords::to_1d(col, flipped_y, z);
	  
#ifdef BLOCK_PARTITION
              id = GridCoords::linear_to_block(id);
#endif
	      cells.emplace_back(id, st);
          }
      }
  }
  */
  return cells;
}

int Reef::load_bmp_file() {
    SLOG("Current directory:", std::filesystem::current_path(), "\n");

    std::vector<std::vector<uint8_t>> substrate_array = readBMPColorMap(_options->substrate_bitmap_path);
    debugColorMapData(_options->substrate_bitmap_path, substrate_array);

    int height = substrate_array.size();
    int width = substrate_array[0].size();
    int num_ecosystem_cells = 0;

    ecosystem_cells.resize(width * height, SubstrateType::NONE);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            uint8_t code = substrate_array[y][x];
            SubstrateType type;

            /*switch (code) {
	    case 1: type = SubstrateType::NONE; break;
	    case 2: type = SubstrateType::SAND; break;
	    case 3: type = SubstrateType::ALGAE; break;
	    case 4: type = SubstrateType::CORAL; break;
	    default:
	      DIE("Unexpected substrate code ", code, " at (", x, ",", y, ")");
            }
  */
      switch (code) {
        case 1: type = SubstrateType::NONE;            break; // black
        case 2: type = SubstrateType::CORAL_NO_ALGAE;  break; // blue
        case 3: type = SubstrateType::SAND_WITH_ALGAE; break; // green
        case 4: type = SubstrateType::CORAL_WITH_ALGAE;break; // red
        case 5: type = SubstrateType::SAND_NO_ALGAE;   break; // yellow (NEW)
        default: DIE("Unexpected substrate code ", code, " at (", x, ",", y, ")");
      }

            int id = GridCoords::to_1d(x, y, 0);  // z=0 for 2D substrate
#ifdef BLOCK_PARTITION
            id = GridCoords::linear_to_block(id);
#endif
            ecosystem_cells[id] = type;
            num_ecosystem_cells++;
        }
    }

    return num_ecosystem_cells;
}

int Reef::load_data_file(const string &fname, int num_grid_points, SubstrateType substrate_type) {
  ifstream f(fname, ios::in | ios::binary);
  if (!f) SDIE("Couldn't open file ", fname);
  f.seekg(0, ios::end);
  auto fsize = f.tellg();
  auto num_ids = fsize / sizeof(int);
  if (num_ids > num_grid_points + 3) DIE("Too many ids in ", fname, " max is ", num_grid_points);
  f.clear();
  f.seekg(0, ios::beg);
  vector<int> id_buf(num_ids);
  if (!f.read(reinterpret_cast<char *>(&(id_buf[0])), fsize))
    DIE("Couldn't read all bytes in ", fname);
  int num_ecosystem_cells = 0;
  // skip first three wwhich are dimensions
  for (int i = 3; i < id_buf.size(); i++) {
    auto id = id_buf[i];
#ifdef BLOCK_PARTITION
    id = GridCoords::linear_to_block(id);
#endif
    ecosystem_cells[id] = substrate_type;
    num_ecosystem_cells++;
  }
  f.close();
  return num_ecosystem_cells;
}

intrank_t Reef::get_rank_for_grid_point(int64_t grid_i) {
  int64_t block_i = grid_i / _grid_blocks.block_size;
  return block_i % rank_n();
}

GridPoint *Reef::get_local_grid_point(grid_points_t &grid_points, int64_t grid_i) {
  int64_t block_i = grid_i / _grid_blocks.block_size / rank_n();
  int64_t i = grid_i % _grid_blocks.block_size + block_i * _grid_blocks.block_size;
  assert(i < grid_points->size());
  GridPoint *grid_point = &(*grid_points)[i];
  if (grid_point->coords.to_1d() != grid_i)
    DIE("mismatched coords to grid i ", grid_point->coords.to_1d(), " != ", grid_i);
  return grid_point;
}

SampleData Reef::get_grid_point_sample_data(int64_t grid_i) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, int64_t grid_i) {
               GridPoint *grid_point = Reef::get_local_grid_point(grid_points, grid_i);
               SampleData sample;
               if (grid_point->fish) {
		 sample.has_fish = true;
		 sample.fishes = 1;
		 sample.fish_alert = grid_point->fish->alert;
		 sample.fish_type = grid_point->fish->type;
		 if (sample.fish_type != FishType::NONE){
		   //SLOG("SampleData::get_grid_point(): fish->type ", to_string(grid_point->fish->type),"\n");
		   //SLOG("SampleData::get_grid_point(): sample.fish_type ", to_string(sample.fish_type),"\n");
		   //SLOG("SampleData::get_grid_point(): sample.fish_alert ", to_string(sample.fish_alert),"\n");
		 }
	       }
               if (grid_point->substrate) {
                 sample.has_substrate = true;
                 sample.substrate_status = grid_point->substrate->status;
                 sample.substrate_type = grid_point->substrate->type;
               }
               sample.floating_algaes = grid_point->floating_algaes;
               sample.chemokine = grid_point->chemokine;
               return sample;
             },
             grid_points, grid_i)
      .wait();
}

// Original simcov get neighbor function
vector<int64_t> *Reef::get_neighbors(GridCoords c) {
  GridPoint *grid_point = Reef::get_local_grid_point(grid_points, c.to_1d());
  if (!grid_point->neighbors) {
    grid_point->neighbors = new vector<int64_t>;
    int newx, newy, newz;
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          newx = c.x + i;
          newy = c.y + j;
          newz = c.z + k;
          if ((newx >= 0 && newx < _grid_size->x) && (newy >= 0 && newy < _grid_size->y) &&
              (newz >= 0 && newz < _grid_size->z)) {
            if (newx != c.x || newy != c.y || newz != c.z) {
              grid_point->neighbors->push_back(GridCoords::to_1d(newx, newy, newz));
            }
          }
        }
      }
    }
  }
  return grid_point->neighbors;
}

/**
 * @brief Get all grid points within a given radius of a location.
 *
 * Returns a list of 1D indices for all grid points within `radius` of the
 * input coordinates `c`, including the centre itself. The neighbourhood shape
 * depends on `metric`:
 *
 *   - Chebyshev: max(|dx|, |dy|, |dz|) ≤ radius   → square (2D) or cube (3D)
 *   - Manhattan: |dx| + |dy| + |dz| ≤ radius      → diamond (2D) or octahedron (3D)
 *   - Euclidean: dx² + dy² + dz² ≤ radius²        → disc (2D) or sphere (3D)
 *
 * Grid boundaries are enforced (no wrap-around in this implementation).
 *
 * Performance:
 *   Time:  O(radius³) in 3D (or O(radius²) in 2D)
 *   Space: O(neighbour count) (result vector is returned by value, RVO applies)
 *
 * Assumptions:
 *   - `_grid_size->x/y/z` define grid dimensions
 *   - `GridCoords::to_1d()` converts (x,y,z) to a 1D index
 *
 * @param c       Centre coordinates
 * @param radius  Non-negative integer radius
 * @param metric  Radius metric type
 * @return vector<int64_t> containing all valid neighbour indices
 */
vector<int64_t>
Reef::get_neighbors(GridCoords c, int radius, RadiusMetric metric) const {
    vector<int64_t> result;

    // Special case: radius ≤ 0 → only the centre cell
    if (radius <= 0) {
        result.push_back(GridCoords::to_1d(c.x, c.y, c.z));
        return result;
    }

    // Estimate max size for efficiency
    int span = 2 * radius + 1;
    bool is3D = (_grid_size->z > 1);
    size_t worst_case = is3D
        ? static_cast<size_t>(span) * span * span
        : static_cast<size_t>(span) * span;
    result.reserve(worst_case);

    // Helper: check if (x,y,z) is inside grid bounds
    auto in_bounds = [&](int x, int y, int z) {
        return (x >= 0 && x < _grid_size->x) &&
               (y >= 0 && y < _grid_size->y) &&
               (z >= 0 && z < _grid_size->z);
    };

    // Helper: check if offset (dx,dy,dz) satisfies radius metric
    auto within_radius = [&](int dx, int dy, int dz) {
        switch (metric) {
            case RadiusMetric::Chebyshev:
                return max({abs(dx), abs(dy), abs(dz)}) <= radius;
            case RadiusMetric::Manhattan:
                return (abs(dx) + abs(dy) + abs(dz)) <= radius;
            case RadiusMetric::Euclidean:
                return (dx*dx + dy*dy + dz*dz) <= radius * radius;
        }
        return false; // Should never happen
    };

    int zmin = is3D ? -radius : 0;
    int zmax = is3D ?  radius : 0;

    // Iterate over cubic bounding box centred at (c.x,c.y,c.z)
    for (int dz = zmin; dz <= zmax; ++dz) {
        for (int dy = -radius; dy <= radius; ++dy) {
            for (int dx = -radius; dx <= radius; ++dx) {
                if (!within_radius(dx, dy, dz)) continue;

                int newx = c.x + dx;
                int newy = c.y + dy;
                int newz = c.z + dz;

                if (in_bounds(newx, newy, newz)) {
                    result.push_back(GridCoords::to_1d(newx, newy, newz));
                }
            }
        }
    }

    return result; // RVO ensures no extra copy
}

int64_t Reef::get_num_local_grid_points() { return grid_points->size(); }

/*
int64_t Reef::get_random_airway_substrate_location() {
  std::set<int>::iterator it = airway.begin();
  std::advance(it, airway.size() / 2);
  return *it;
}
*/

bool Reef::set_initial_infection(int64_t grid_i) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                int64_t grid_i) {
               GridPoint *grid_point = Reef::get_local_grid_point(grid_points, grid_i);
               DBG("set infected for grid point ", grid_point, " ", grid_point->str(), "\n");
               if (!grid_point->substrate) return false;
               grid_point->floating_algaes = _options->initial_infection;
               new_active_grid_points->insert({grid_point, true});
               return true;
             },
             grid_points, new_active_grid_points, grid_i)
      .wait();
}

void Reef::accumulate_chemokines(HASH_TABLE<int64_t, float> &chemokines_to_update,
                                   IntermittentTimer &timer) {
  timer.start();
  // accumulate updates for each target rank
  HASH_TABLE<intrank_t, vector<pair<int64_t, float>>> target_rank_updates;
  for (auto &[coords_1d, chemokines] : chemokines_to_update) {
    progress();
    target_rank_updates[get_rank_for_grid_point(coords_1d)].push_back({coords_1d, chemokines});
  }
  future<> fut_chain = make_future<>();
  // dispatch all updates to each target rank in turn
  for (auto &[target_rank, update_vector] : target_rank_updates) {
    progress();
    auto fut = rpc(
        target_rank,
        [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
           view<pair<int64_t, float>> update_vector) {
          for (auto &[grid_i, chemokine] : update_vector) {
            GridPoint *grid_point = Reef::get_local_grid_point(grid_points, grid_i);
            new_active_grid_points->insert({grid_point, true});
            // just accumulate the concentrations. We will adjust them to be the average
            // of all neighbors later
            grid_point->nb_chemokine += chemokine;
          }
        },
        grid_points, new_active_grid_points, make_view(update_vector));
    fut_chain = when_all(fut_chain, fut);
  }
  fut_chain.wait();
  timer.stop();
}

void Reef::accumulate_floating_algaes(HASH_TABLE<int64_t, float> &floating_algaes_to_update,
                                IntermittentTimer &timer) {
  timer.start();
  // accumulate updates for each target rank
  HASH_TABLE<intrank_t, vector<pair<int64_t, float>>> target_rank_updates;
  for (auto &[coords_1d, floating_algaes] : floating_algaes_to_update) {
    progress();
    target_rank_updates[get_rank_for_grid_point(coords_1d)].push_back({coords_1d, floating_algaes});
  }
  future<> fut_chain = make_future<>();
  // dispatch all updates to each target rank in turn
  for (auto &[target_rank, update_vector] : target_rank_updates) {
    progress();
    auto fut = rpc(
        target_rank,
        [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
           view<pair<int64_t, float>> update_vector) {
          for (auto &[grid_i, floating_algaes] : update_vector) {
            GridPoint *grid_point = Reef::get_local_grid_point(grid_points, grid_i);
            new_active_grid_points->insert({grid_point, true});
            grid_point->nb_floating_algaes += floating_algaes;
          }
        },
        grid_points, new_active_grid_points, make_view(update_vector));
    fut_chain = when_all(fut_chain, fut);
  }
  fut_chain.wait();
  timer.stop();
}

float Reef::get_chemokine(int64_t grid_i) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, int64_t grid_i) {
               GridPoint *grid_point = Reef::get_local_grid_point(grid_points, grid_i);
               return grid_point->chemokine;
             },
             grid_points, grid_i)
      .wait();
}

int64_t Reef::get_num_circulating_fishes() { return num_circulating_fishes; }

void Reef::change_num_circulating_fishes(int num) {
  num_circulating_fishes += num;
  if (num_circulating_fishes < 0) num_circulating_fishes = 0;
}

bool Reef::try_add_new_reef_fish(int64_t grid_i) {
  auto res = rpc(
                 get_rank_for_grid_point(grid_i),
                 [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                    int64_t grid_i, dist_object<int64_t> &fishes_generated) {
                   GridPoint *grid_point = Reef::get_local_grid_point(grid_points, grid_i);
                   // grid point is already occupied by a fish, don't add
                   if (grid_point->fish) return false;
                   //if (grid_point->chemokine < _options->min_chemokine) return false;
                   new_active_grid_points->insert({grid_point, true});
                   string fish_id = to_string(rank_me()) + "-" + to_string(*fishes_generated);
                   (*fishes_generated)++;
                   grid_point->fish = new Fish(fish_id);
                   grid_point->fish->angle = sample_vonmises(0.0, 0.0, vonmises_gen);
                   grid_point->fish->moved = true;
                   // Set current fish coords to grid point coords
                   grid_point->fish->x = grid_point->coords.x;
                   grid_point->fish->y = grid_point->coords.y;
                   grid_point->fish->z = grid_point->coords.z;
		   grid_point->substrate->status = SubstrateStatus::FISH;

		   grid_point->fish->type = FishType::GRAZER;
		   
		   // Some fraction of the fish will be predators
		   if (_options->predator_ratio > 0) 
		     if (!(rand() % _options->predator_ratio)) // 1/ratio_predators chance
		       grid_point->fish->type = FishType::PREDATOR;
		   		   
                   return true;
                 },
                 grid_points, new_active_grid_points, grid_i, fishes_generated)
                 .wait();
  //if (res) num_circulating_fishes--;
  //assert(num_circulating_fishes >= 0);
  //return res;
  return 0;
}

bool Reef::try_add_reef_fish(int64_t grid_i, Fish &fish) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                int64_t grid_i, Fish fish) {
               GridPoint *grid_point = Reef::get_local_grid_point(grid_points, grid_i);
               // grid point is already occupied by a fish, don't add
               if (grid_point->fish) return false;
               new_active_grid_points->insert({grid_point, true});
               fish.moved = true;
               grid_point->fish = new Fish(fish);
               // Set current fish coords to grid point coords
               grid_point->fish->x = grid_point->coords.x;
               grid_point->fish->y = grid_point->coords.y;
               grid_point->fish->z = grid_point->coords.z;
	       grid_point->fish->type = fish.type;
	       grid_point->fish->alert = fish.alert;
	       
	       grid_point->substrate->status = SubstrateStatus::FISH;
               return true;
             },
             grid_points, new_active_grid_points, grid_i, fish)
      .wait();
}

SubstrateStatus Reef::try_bind_fish(int64_t grid_i) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points_t,
                int64_t grid_i) {
               GridPoint *grid_point = Reef::get_local_grid_point(grid_points, grid_i);
               if (!grid_point->substrate) return SubstrateStatus::DEAD;
               if (grid_point->substrate->status == SubstrateStatus::HEALTHY ||
                   grid_point->substrate->status == SubstrateStatus::DEAD)
                 return grid_point->substrate->status;

               // if (grid_point->substrate->status == SubstrateStatus::DEAD) return
               // SubstrateStatus::DEAD;

               double binding_prob = grid_point->substrate->get_binding_prob();
               if (_rnd_gen->trial_success(binding_prob)) {
                 auto prev_status = grid_point->substrate->status;
                 grid_point->substrate->status = SubstrateStatus::APOPTOTIC;
                 return prev_status;
               }
               return SubstrateStatus::DEAD;
             },
             grid_points, new_active_grid_points, grid_i)
      .wait();
}

GridPoint *Reef::get_first_local_grid_point() {
  grid_point_iter = grid_points->begin();
  if (grid_point_iter == grid_points->end()) return nullptr;
  auto grid_point = &(*grid_point_iter);
  ++grid_point_iter;
  return grid_point;
}

GridPoint *Reef::get_next_local_grid_point() {
  if (grid_point_iter == grid_points->end()) return nullptr;
  auto grid_point = &(*grid_point_iter);
  ++grid_point_iter;
  return grid_point;
}

GridPoint *Reef::get_first_active_grid_point() {
  active_grid_point_iter = active_grid_points.begin();
  if (active_grid_point_iter == active_grid_points.end()) return nullptr;
  auto grid_point = active_grid_point_iter->first;
  ++active_grid_point_iter;
  return grid_point;
}

GridPoint *Reef::get_next_active_grid_point() {
  if (active_grid_point_iter == active_grid_points.end()) return nullptr;
  auto grid_point = active_grid_point_iter->first;
  ++active_grid_point_iter;
  return grid_point;
}

void Reef::set_active(GridPoint *grid_point) {
  new_active_grid_points->insert({grid_point, true});
}

void Reef::erase_active(GridPoint *grid_point) { active_grid_points.erase(grid_point); }

void Reef::add_new_actives(IntermittentTimer &timer) {
  timer.start();
  DBG("add ", new_active_grid_points->size(), " new active grid points\n");
  for (auto elem : *new_active_grid_points) {
    // DBG("inserting from new active ", elem.first, " ", elem.first->str(), "\n");
    active_grid_points.insert(elem);
  }
  new_active_grid_points->clear();
  timer.stop();
}

size_t Reef::get_num_actives() { return active_grid_points.size(); }

// Return true if a fish of the specified type is within the specified radius
bool Reef::detect_neighbour_fish(const GridPoint* center,
                                 FishType fish_type,
                                 int radius,
                                 RadiusMetric metric)
{
    if (!center) return false;

    auto neighbour_indices = get_neighbors(center->coords, radius, metric);

    // Iterate over the neighbourhood looking for a matching fish type
    for (auto idx : neighbour_indices) {
        GridPoint* gp = Reef::get_local_grid_point(grid_points, idx);
        if (!gp) continue;

	if (gp->fish && gp->fish->type == fish_type) {
	  return true;  // match found
        }
	
    }

    // No matching fish of the required type found
    return false;
}

int Reef::count_neighbour_substrate(const GridPoint* center,
                                   SubstrateType type,
                                   int radius,
                                   RadiusMetric metric)
{
    if (!center) return 0;

    // Neighbour indices (includes center by design)
    auto ids = get_neighbors(center->coords, radius, metric);

    int count = 0;
    for (auto idx : ids) {
        GridPoint* gp = Reef::get_local_grid_point(grid_points, idx);
        if (gp && gp->substrate && gp->substrate->type == type) {
            ++count;
        }
    }
    return count;
}


#ifdef DEBUG
void Reef::check_actives(int time_step) {
  for (int64_t i = 0; i < grid_points->size(); i++) {
    GridPoint *grid_point = &(*grid_points)[i];
    if (time_step == 0)
      DBG("coords ", grid_point->coords.str(), " to 1d ", grid_point->coords.to_1d(), "\n");
    bool found_active = (active_grid_points.find(grid_point) != active_grid_points.end());
    if (grid_point->is_active() && !found_active)
      DIE(time_step, ": active grid point ", grid_point->str(), " not found in active list");
    if (!grid_point->is_active() && found_active)
      DIE(time_step, ": inactive grid point ", grid_point->str(), " found in active list");
  }
}
#endif
