# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# C++-related tests for SUNDIALS CMake-based configuration.
# ---------------------------------------------------------------

# Enable the C++
enable_language(CXX)
set(CXX_FOUND TRUE)

# show some cache variables
mark_as_advanced(CLEAR
  CMAKE_CXX_COMPILER
  CMAKE_CXX_FLAGS)

# hide all build type specific flags
mark_as_advanced(FORCE
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_RELWITHDEBINFO)

# only show flags for the current build type
# these flags are appended to CMAKE_CXX_FLAGS
if(CMAKE_BUILD_TYPE)
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    message("Appending CXX debug flags")
    mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_DEBUG)
  elseif(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
    message("Appending CXX min size release flags")
    mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_MINSIZEREL)
  elseif(CMAKE_BUILD_TYPE MATCHES "Release")
    message("Appending CXX release flags")
    mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_RELEASE)
  elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    message("Appending CXX release with debug info flags")
    mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_RELWITHDEBINFO)
  endif()
endif()
