simforager is a forked adaptation of:

simcov Copyright (c) 2021, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy), Arizona State
University and University of New Mexico. All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

# SimForager #

This is a model for simulating large amounts of [currently homogeneous] foragers.

## Installing and building

It requires [UPC++](https://bitbucket.org/berkeleylab/upcxx/wiki/Home), C++ and cmake.

This repo contains a submodule, so to install, it's best to run

`git clone --recurse-submodules git@github.com:cswritlarge/simforager.git`

to fully initialize the submodules.

Once downloaded, build it with

`./build.sh Release`

or

`./build.sh Debug`

The executable will be installed in

`<simforager-repo-directory>/install/bin`

## Running

To run, execute

`upcxx-run -n <number_processes> -N <number_nodes> -- simforager`

To see the parameters available, run with `-h`.

It will create an output directory, which will contain a detailed log file (`simforager.log`). It will also create a file containing all the configuration parameters (`simforager.config`).

A run can also be executed with a config file as

`upcxx-run -n <number_processes> -N <number_nodes> -- simforager --config <config_file>`

The config file consists of a list of all the command line options as key-value pairs, with semi-colons denoting comments.
For example:

```
; Dimensions: x y z
  dim =                    100,100,100

; Number of timesteps
  timesteps =                   14000

```

Any options not specified in the config file will be set to the defaults (seen when run with `-h`).

When running with a config file, any parameters passed on the command line will override parameters read from the config file.


If running in debug mode, a subdirectory
called `per_thread` will appear, with one directory per process that contains debug information produced by that proc

Once the run has completed, the outputs can be viewed in paraview by opening the samples subdirectory of the output directory. To help with viewing, a python script can generate a paraview state file:

```
scripts/generate_paraview_state.py --data <output_dir>/samples --stats <output_dir>/simcov.stats -o paraview-state
```

Note that the `pvpython` wrapper is needed (not just plain python), and is used in the script.

The above command will generate a file called `paraview-state.pvsm` and this can be loaded into paraview, e.g.:

```
paraview paraview-state.pvsm
```

## Versions

The version of simcov at commit id a7aa913288bd8139cd6fef5410811c2fad097236 on Aug 30, 2021 was used with pre-print publication https://www.biorxiv.org/content/10.1101/2021.05.19.444569v3 and the published version https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009735.

Note that the config folders in this version of the source code are numbered according to the figure order corresponding to the preprint, but different from the final published version.
