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
include upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/depend.make

# Include the progress variables for this target.
include upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/progress.make

# Include the compile flags for this target's objects.
include upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/flags.make

upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.o: upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/flags.make
upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.o: ../upcxx-utils/src/mem_profile.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/mfricke/simreef/.build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.o"
	cd /users/mfricke/simreef/.build/upcxx-utils/src && /opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/upcxx-2020.10.0-6eh2prmiaolqfqinq4wjbb5by6z2phw6/bin/upcxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.o -c /users/mfricke/simreef/upcxx-utils/src/mem_profile.cpp

upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.i"
	cd /users/mfricke/simreef/.build/upcxx-utils/src && /opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/upcxx-2020.10.0-6eh2prmiaolqfqinq4wjbb5by6z2phw6/bin/upcxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/mfricke/simreef/upcxx-utils/src/mem_profile.cpp > CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.i

upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.s"
	cd /users/mfricke/simreef/.build/upcxx-utils/src && /opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/upcxx-2020.10.0-6eh2prmiaolqfqinq4wjbb5by6z2phw6/bin/upcxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/mfricke/simreef/upcxx-utils/src/mem_profile.cpp -o CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.s

upcxx_utils_mem_profile: upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/mem_profile.cpp.o
upcxx_utils_mem_profile: upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/build.make

.PHONY : upcxx_utils_mem_profile

# Rule to build all files generated by this target.
upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/build: upcxx_utils_mem_profile

.PHONY : upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/build

upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/clean:
	cd /users/mfricke/simreef/.build/upcxx-utils/src && $(CMAKE_COMMAND) -P CMakeFiles/upcxx_utils_mem_profile.dir/cmake_clean.cmake
.PHONY : upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/clean

upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/depend:
	cd /users/mfricke/simreef/.build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /users/mfricke/simreef /users/mfricke/simreef/upcxx-utils/src /users/mfricke/simreef/.build /users/mfricke/simreef/.build/upcxx-utils/src /users/mfricke/simreef/.build/upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : upcxx-utils/src/CMakeFiles/upcxx_utils_mem_profile.dir/depend
