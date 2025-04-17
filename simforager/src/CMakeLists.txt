
# CMakeLists.txt — Build Configuration for SimForager
# ----------------------------------------------------
# Defines build rules and dependencies for compiling the SimForager project with UPC++.

cmake_minimum_required(VERSION 3.10)
project(simforager CXX)

# Enable C++14 for compatibility with UPC++ and STL
set(CMAKE_CXX_STANDARD 14)

# Compiler warnings and debug symbols
add_compile_options(-Wall -Wextra -O2 -g)

# Find UPC++ installation (required)
find_package(UPCXX REQUIRED)

# Define the main executable target
add_executable(simforager
  main.cpp
  reef.cpp
  utils.cpp
  vonmises.cpp
)

# Link against UPC++ libraries
target_link_libraries(simforager PRIVATE UPCXX::upcxx)

# Optionally pass macro definitions (e.g. for GPU, debug flags)
# add_definitions(-DUSE_CUDA)
