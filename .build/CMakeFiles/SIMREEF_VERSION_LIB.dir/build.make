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
include CMakeFiles/SIMREEF_VERSION_LIB.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SIMREEF_VERSION_LIB.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SIMREEF_VERSION_LIB.dir/flags.make

# Object files for target SIMREEF_VERSION_LIB
SIMREEF_VERSION_LIB_OBJECTS =

# External object files for target SIMREEF_VERSION_LIB
SIMREEF_VERSION_LIB_EXTERNAL_OBJECTS = \
"/users/mfricke/simreef/.build/CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.o"

libSIMREEF_VERSION_LIB.a: CMakeFiles/SIMREEF_VERSION.dir/makeVersionFile/version.cpp.o
libSIMREEF_VERSION_LIB.a: CMakeFiles/SIMREEF_VERSION_LIB.dir/build.make
libSIMREEF_VERSION_LIB.a: CMakeFiles/SIMREEF_VERSION_LIB.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/users/mfricke/simreef/.build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX static library libSIMREEF_VERSION_LIB.a"
	$(CMAKE_COMMAND) -P CMakeFiles/SIMREEF_VERSION_LIB.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SIMREEF_VERSION_LIB.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SIMREEF_VERSION_LIB.dir/build: libSIMREEF_VERSION_LIB.a

.PHONY : CMakeFiles/SIMREEF_VERSION_LIB.dir/build

CMakeFiles/SIMREEF_VERSION_LIB.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SIMREEF_VERSION_LIB.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SIMREEF_VERSION_LIB.dir/clean

CMakeFiles/SIMREEF_VERSION_LIB.dir/depend:
	cd /users/mfricke/simreef/.build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /users/mfricke/simreef /users/mfricke/simreef /users/mfricke/simreef/.build /users/mfricke/simreef/.build /users/mfricke/simreef/.build/CMakeFiles/SIMREEF_VERSION_LIB.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SIMREEF_VERSION_LIB.dir/depend
