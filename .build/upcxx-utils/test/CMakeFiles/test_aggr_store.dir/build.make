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
include upcxx-utils/test/CMakeFiles/test_aggr_store.dir/depend.make

# Include the progress variables for this target.
include upcxx-utils/test/CMakeFiles/test_aggr_store.dir/progress.make

# Include the compile flags for this target's objects.
include upcxx-utils/test/CMakeFiles/test_aggr_store.dir/flags.make

upcxx-utils/test/CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.o: upcxx-utils/test/CMakeFiles/test_aggr_store.dir/flags.make
upcxx-utils/test/CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.o: ../upcxx-utils/test/test_aggr_store.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/mfricke/simreef/.build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object upcxx-utils/test/CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.o"
	cd /users/mfricke/simreef/.build/upcxx-utils/test && /opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/upcxx-2020.10.0-6eh2prmiaolqfqinq4wjbb5by6z2phw6/bin/upcxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.o -c /users/mfricke/simreef/upcxx-utils/test/test_aggr_store.cpp

upcxx-utils/test/CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.i"
	cd /users/mfricke/simreef/.build/upcxx-utils/test && /opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/upcxx-2020.10.0-6eh2prmiaolqfqinq4wjbb5by6z2phw6/bin/upcxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/mfricke/simreef/upcxx-utils/test/test_aggr_store.cpp > CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.i

upcxx-utils/test/CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.s"
	cd /users/mfricke/simreef/.build/upcxx-utils/test && /opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/upcxx-2020.10.0-6eh2prmiaolqfqinq4wjbb5by6z2phw6/bin/upcxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/mfricke/simreef/upcxx-utils/test/test_aggr_store.cpp -o CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.s

# Object files for target test_aggr_store
test_aggr_store_OBJECTS = \
"CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.o"

# External object files for target test_aggr_store
test_aggr_store_EXTERNAL_OBJECTS =

upcxx-utils/test/test_aggr_store: upcxx-utils/test/CMakeFiles/test_aggr_store.dir/test_aggr_store.cpp.o
upcxx-utils/test/test_aggr_store: upcxx-utils/test/CMakeFiles/test_aggr_store.dir/build.make
upcxx-utils/test/test_aggr_store: upcxx-utils/src/libUPCXX_UTILS.a
upcxx-utils/test/test_aggr_store: upcxx-utils/test/CMakeFiles/test_aggr_store.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/users/mfricke/simreef/.build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_aggr_store"
	cd /users/mfricke/simreef/.build/upcxx-utils/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_aggr_store.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
upcxx-utils/test/CMakeFiles/test_aggr_store.dir/build: upcxx-utils/test/test_aggr_store

.PHONY : upcxx-utils/test/CMakeFiles/test_aggr_store.dir/build

upcxx-utils/test/CMakeFiles/test_aggr_store.dir/clean:
	cd /users/mfricke/simreef/.build/upcxx-utils/test && $(CMAKE_COMMAND) -P CMakeFiles/test_aggr_store.dir/cmake_clean.cmake
.PHONY : upcxx-utils/test/CMakeFiles/test_aggr_store.dir/clean

upcxx-utils/test/CMakeFiles/test_aggr_store.dir/depend:
	cd /users/mfricke/simreef/.build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /users/mfricke/simreef /users/mfricke/simreef/upcxx-utils/test /users/mfricke/simreef/.build /users/mfricke/simreef/.build/upcxx-utils/test /users/mfricke/simreef/.build/upcxx-utils/test/CMakeFiles/test_aggr_store.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : upcxx-utils/test/CMakeFiles/test_aggr_store.dir/depend

