# Install script for directory: /users/mfricke/simreef/upcxx-utils

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/users/mfricke/simreef/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/UPCXX_UTILS.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/UPCXX_UTILS.cmake"
         "/users/mfricke/simreef/.build/upcxx-utils/CMakeFiles/Export/cmake/UPCXX_UTILS.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/UPCXX_UTILS-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/UPCXX_UTILS.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "/users/mfricke/simreef/.build/upcxx-utils/CMakeFiles/Export/cmake/UPCXX_UTILS.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "/users/mfricke/simreef/.build/upcxx-utils/CMakeFiles/Export/cmake/UPCXX_UTILS-release.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "/users/mfricke/simreef/upcxx-utils/cmake/UPCXX_UTILSConfig.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES
    "/users/mfricke/simreef/.build/upcxx-utils/makeVersionFile/UPCXX_UTILSConfigVersion.cmake"
    "/users/mfricke/simreef/.build/upcxx-utils/makeVersionFile/UPCXX_UTILS_VERSION"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/upcxx_utils" TYPE FILE FILES
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/colors.h"
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/flat_aggr_store.hpp"
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/log.hpp"
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/mem_profile.hpp"
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/progress_bar.hpp"
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/split_rank.hpp"
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/timers.hpp"
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/two_tier_aggr_store.hpp"
    "/users/mfricke/simreef/upcxx-utils/include/upcxx_utils/version.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE FILE FILES
    "/users/mfricke/simreef/upcxx-utils/LICENSE.txt"
    "/users/mfricke/simreef/upcxx-utils/LEGAL.txt"
    "/users/mfricke/simreef/upcxx-utils/README.md"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/users/mfricke/simreef/.build/upcxx-utils/src/cmake_install.cmake")
  include("/users/mfricke/simreef/.build/upcxx-utils/test/cmake_install.cmake")

endif()

