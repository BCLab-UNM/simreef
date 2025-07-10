cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

message(STATUS "finding UPCXX for SimCov build")
find_package(UPCXX REQUIRED)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13 AND DEFINED UPCXX_LIBRARIES)
  message(STATUS "UPCXX_UTILS is using the UPCXX::upcxx library interface: ${UPCXX_LIBRARIES}")
else()
  find_program(UPCXX_EXEC upcxx)
  set(CMAKE_CXX_COMPILER ${UPCXX_EXEC})
  message(STATUS "UPCXX_UTILS is using upcxx directly as the UPCXX::upcxx library interface is not available in this low version of cmake: ${CMAKE_VERSION}")
  unset(UPCXX_LIBRARIES)
endif()

message("Building ${CMAKE_BUILD_TYPE} version")

find_package(Threads REQUIRED)
find_package(OpenCV REQUIRED)
message(STATUS "OpenCV_INCLUDE_DIRS = ${OpenCV_INCLUDE_DIRS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-function -Wno-catch-value -pthread")

# Collect compile-time definitions
get_directory_property(dir_defs COMPILE_DEFINITIONS)
set(all_defs "")
message("Definitions:")
foreach(d ${dir_defs})
    message("  " ${d})
    set(all_defs "${all_defs} ${d}")
endforeach()
add_definitions(-DALL_DEFNS=\"${all_defs}\")

# Main target
add_executable(simforager main.cpp reef.cpp utils.cpp vonmises.cpp)

# OpenCV include (accounting for opencv4 layout if necessary)
target_include_directories(simforager PRIVATE
    ${OpenCV_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}/opencv4
)

# Link libraries using plain-style linking
target_link_libraries(simforager
    Threads::Threads
    ${UPCXX_LIBRARIES}
    ${UPCXX_UTILS_LIBRARIES}
    SIMFORAGER_VERSION_LIB
    ${OpenCV_LIBS}
)

# Install the binary
install(TARGETS simforager DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
