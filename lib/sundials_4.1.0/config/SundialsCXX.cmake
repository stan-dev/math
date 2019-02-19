# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU, Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
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

# Converts a CXX standard number to the flag needed to set it.
macro(CXX_STD2FLAG std flag_var)
  set(flag_var "-std=c++${std}")
endmacro()

# Sets the CMAKE_CXX_STANDARD variable to the ${std} if the compiler
# supports the flag. E.g. USE_CXX_STD(11) sets CMAKE_CXX_STANDARD=11.
# If CUDA is enabled, it adds the correct flag to CUDA_NVCC_FLAGS.
#
# Requires:
#   CMake > 3.1.3.
# Notes: 
#   If the compiler is not supprted by the CMake version in use, then
#   the flag will have to be added manually.
macro(USE_CXX_STD std)
  if(NOT (CMAKE_CXX_STANDARD EQUAL ${std}))
    include(CheckCXXCompilerFlag)
    CXX_STD2FLAG(${std} flag_var)
    CHECK_CXX_COMPILER_FLAG(${flag_var} COMPILER_SUPPORTS_STDFLAG)
    if(COMPILER_SUPPORTS_STDFLAG)
      set(CMAKE_CXX_STANDARD ${std})
      message(STATUS "Set CMAKE_CXX_STANDARD to ${std}")
      if(CUDA_ENABLE)
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std c++${std}")
        message(STATUS "Add -std c++${std} to CUDA_NVCC_FLAGS")
      endif()
    else()
      PRINT_WARNING("Could not set CMAKE_CXX_STANDARD to ${std}.")
    endif()
  endif()
endmacro()

