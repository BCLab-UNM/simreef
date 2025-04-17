
# SimForager

SimForager is a high-performance distributed simulation framework designed to model biological or ecological processes on a 2D or 3D spatial grid. It uses UPC++ for scalable parallelism and features biologically inspired agents ("fishes") that interact with substrates in an evolving environment.

Originally inspired by SimCov, SimForager supports simulation of infection-like dynamics, chemical signaling (chemokines), and diffusive interactions across a discretised ecosystem grid, including block-wise partitioning for parallel decomposition.

---

## 🧠 Key Concepts

- **Grid-based environment**: Each grid point may contain a substrate (e.g., alveoli, coral), floating particles (algae), and/or fishes.
- **Agents (Fishes)**: Mobile units with stochastic movement, potentially governed by gradient following or Von Mises angular drift.
- **Substrates**: Cells with life cycles—healthy, incubating, expressing, apoptotic, or dead—modulated by contact and diffusion.
- **Diffusing fields**: Chemokines and algae spread via neighborhood-based diffusion and decay.
- **Parallel runtime**: Uses UPC++ to distribute simulation logic across processes, including RPC for remote grid updates.

---

## 🗂 File Overview

### `main.cpp`

- Entry point for the simulation.
- Initializes options, reef grid, and agent state.
- Runs the simulation loop with time-stepped updates.
- Handles sampling, diffusion, fish movement, infection, and logging.
- Writes VTK output for visualisation.

### `reef.hpp` / `reef.cpp`

Defines and implements:

- `GridCoords`: Converts between 3D and 1D indexing.
- `GridPoint`: Stores substrate, fish, chemokine, algae at a single cell.
- `Substrate`: Models lifecycle state and transitions.
- `Fish`: Represents mobile agents with turning angles and persistence.
- `Reef`: Central class orchestrating the grid, partitioning, and simulation updates.

### `utils.hpp` / `utils.cpp`

- `readBMPColorMap`: Parses `.bmp` files into grid maps.
- `debugColorMapData`: Validates image-based maps.
- `pin_thread`: Sets thread-core affinity.
- `dump_single_file`: Atomic, rank-safe parallel file writer.

### `vonmises.hpp` / `vonmises.cpp`

- Implements the Von Mises distribution for simulating angular persistence in movement.
- `sample_vonmises`: Rejection sampling of angles.
- `polar_to_cartesian`: Converts angles to movement deltas.

### `options.hpp`

- Contains the `Options` struct with all tunable simulation parameters:
  - Grid size, durations, diffusion rates, output paths, fish behaviour
- Shared across the simulation as `_options`.

### `version.h`

- Declares string constants for the build version, date, and branch.
- Defines macros for CUDA context awareness.

### `CMakeLists.txt`

- CMake build configuration file.
- Uses UPC++ as a required dependency.
- Builds `simforager` binary from source files.

---

## 📦 Third-Party Headers

### `bytell_hash_map.hpp`

- High-performance, memory-compact hash map (likely from ankerl::unordered_dense).
- Used for efficient per-grid-point lookup and indexing.

### `CLI11.hpp`

- Command-line argument parser.
- Parses flags, options, and arguments to populate the simulation configuration.

---

## 🚀 Building and Running

### Prerequisites

- **UPC++**: Must be installed and sourced (`upcxx-meta install`).
- CMake ≥ 3.10
- C++14-compatible compiler

### Building

```bash
mkdir build && cd build
cmake ..
make -j
```

### Running

```bash
upcxx-run -n 4 ./simforager --dimensions 100 100 1 --output_dir results/
```

---

## 🧪 Output and Visualisation

- Statistics are logged to `simforager.stats`
- Spatial samples are output as `.vtk` binary files, compatible with Paraview or VisIt
- Sampling includes substrate states, fish locations, algae, and chemokine concentrations

---

## 📝 License & Acknowledgements

- Code originally by Steven Hofmeyr (LBNL) and extended by Matthew Fricke (UNM)
- Built using UPC++ from Lawrence Berkeley National Laboratory
- Includes headers from [ankerl::unordered_dense](https://github.com/martinus/unordered_dense) and [CLI11](https://github.com/CLIUtils/CLI11)


---

## 🔄 Origin and Licensing

SimForager is a **forked adaptation** of:

> **SimCov**  
> Copyright (c) 2021, The Regents of the University of California,  
> through Lawrence Berkeley National Laboratory (subject to receipt of  
> any required approvals from the U.S. Dept. of Energy), Arizona State  
> University and University of New Mexico. All rights reserved.

> NOTICE. This Software was developed under funding from the U.S. Department of Energy.  
> The U.S. Government retains a paid-up, nonexclusive, irrevocable, worldwide license  
> to reproduce, distribute, prepare derivative works, and perform/display publicly.

For licensing inquiries, contact Berkeley Lab's Intellectual Property Office at  
📧 [IPO@lbl.gov](mailto:IPO@lbl.gov).

---

## 🛠️ Installation and Building

SimForager requires:

- [UPC++](https://bitbucket.org/berkeleylab/upcxx/wiki/Home)
- CMake
- A C++14 compatible compiler

This repository includes a submodule, so use:

```bash
git clone --recurse-submodules git@github.com:cswritlarge/simforager.git
```

Then build using the included script:

```bash
./build.sh Release   # or ./build.sh Debug
```

The resulting binary will be installed in:

```
<simforager-repo-directory>/install/bin
```

---

## ▶️ Running SimForager

Use `upcxx-run` to launch the simulation:

```bash
upcxx-run -n <number_processes> -N <number_nodes> -- simforager
```

To see available parameters:

```bash
simforager -h
```

### Output

- A run creates a detailed log file `simforager.log`
- A config summary is saved to `simforager.config`
- Simulation snapshots are stored in `samples/` (viewable in Paraview)

### Using Config Files

Runs can also be launched with a config file:

```bash
upcxx-run -n <number_processes> -N <number_nodes> -- simforager --config <config_file>
```

Example config format:

```ini
; Dimensions: x y z
  dim = 100,100,100

; Number of timesteps
  timesteps = 14000
```

Any CLI arguments override values in the config file.

### Debug Mode

If compiled with debugging, a `per_thread/` directory is created containing output from each rank.

---

## 📊 Visualising Results in Paraview

To generate a Paraview state file:

```bash
scripts/generate_paraview_state.py --data <output_dir>/samples --stats <output_dir>/simcov.stats -o paraview-state
```

Then open it with:

```bash
paraview paraview-state.pvsm
```

> **Note**: Use `pvpython`, not standard Python, for the script.

---

## 📚 Publication and Version Info

SimCov commit `a7aa9132` (Aug 30, 2021) was used for:

- Preprint: [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.05.19.444569v3)
- Published: [PLOS Comp Bio](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009735)

⚠️ Config folder numbers correspond to the preprint figures and differ from the published version.

