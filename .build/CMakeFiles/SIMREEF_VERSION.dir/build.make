# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/spack/opt/spack/linux-rocky8-skylake_avx512/gcc-8.5.0/cmake-3.11.4-qkyjqkhc6gp4zhxt7rhhf6wy4vf76ocj/bin/cmake

# The command to remove a file.
RM = /opt/spack/opt/spack/linux-rocky8-skylake_avx512/gcc-8.5.0/cmake-3.11.4-qkyjqkhc6gp4zhxt7rhhf6wy4vf76ocj/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /users/mfricke/simreef

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /users/mfricke/simreef/.build

# Include any dependencies generated for this target.
include CMakeFiles/SIMREEF_VERSION.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SIMREEF_VERSION.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SIMREEF_VERSION.dir/flags.make

CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.o: CMakeFiles/SIMREEF_VERSION.dir/flags.make
CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.o: makeVersionFile/version.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/mfricke/simreef/.build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.o"
	/opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/openmpi-4.1.3-j6zbgs4rx7w7mb4imwl6fqk2wxvglehb/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.o -c /users/mfricke/simreef/.build/makeVersionFile/version.cpp

CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.i"
	/opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/openmpi-4.1.3-j6zbgs4rx7w7mb4imwl6fqk2wxvglehb/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/mfricke/simreef/.build/makeVersionFile/version.cpp > CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.i

CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.s"
	/opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/openmpi-4.1.3-j6zbgs4rx7w7mb4imwl6fqk2wxvglehb/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/mfricke/simreef/.build/makeVersionFile/version.cpp -o CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.s

SIMREEF_VERSION: CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.o
SIMREEF_VERSION: CMakeFiles/SIMREEF_VERSION.dir/build.make

.PHONY : SIMREEF_VERSION

# Rule to build all files generated by this target.
CMakeFiles/SIMREEF_VERSION.dir/build: SIMREEF_VERSION

.PHONY : CMakeFiles/SIMREEF_VERSION.dir/build

CMakeFiles/SIMREEF_VERSION.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SIMREEF_VERSION.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SIMREEF_VERSION.dir/clean

CMakeFiles/SIMREEF_VERSION.dir/depend:
	cd /users/mfricke/simreef/.build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /users/mfricke/simreef /users/mfricke/simreef /users/mfricke/simreef/.build /users/mfricke/simreef/.build /users/mfricke/simreef/.build/CMakeFiles/SIMREEF_VERSION.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SIMREEF_VERSION.dir/depend

