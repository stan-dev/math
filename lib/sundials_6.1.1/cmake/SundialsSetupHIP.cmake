# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Setup the HIP language and libraries.
# ---------------------------------------------------------------

if(NOT DEFINED ROCM_PATH)
  if(NOT DEFINED ENV{ROCM_PATH})
    set(ROCM_PATH "/opt/rocm/" CACHE PATH "Path to which ROCm has been installed")
  else()
    set(ROCM_PATH "$ENV{ROCM_PATH}" CACHE PATH "Path to which ROCm has been installed")
  endif()
endif()

if(NOT DEFINED HIP_PATH)
  if(NOT DEFINED ENV{HIP_PATH})
    set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
  else()
    set(HIP_PATH "$ENV{HIP_PATH}" CACHE PATH "Path to which HIP has been installed")
  endif()
endif()

if(NOT DEFINED HIP_PLATFORM)
  if(NOT DEFINED ENV{HIP_PLATFORM})
    set(HIP_PLATFORM "hcc" CACHE STRING "HIP platform (hcc, nvcc)")
  else()
    set(HIP_PLATFORM "$ENV{HIP_PLATFORM}" CACHE STRING "HIP platform (hcc, nvcc)")
  endif()
endif()

# Set CMAKE_PREFIX_PATH as the hip-config.cmake has some find_package calls
# which don't have the proper path set (not sure if this is a bug or
# intentional), so without this they will fail even if we provide the PATH
# option to find_package(HIP).
set(CMAKE_PREFIX_PATH "${ROCM_PATH};${HIP_PATH}")
find_package(HIP REQUIRED)

if("${HIP_COMPILER}" STREQUAL "hcc")
  print_error("Deprecated HCC compiler is not supported" "Please update ROCm")
endif()

message(STATUS "HIP version:      ${HIP_VERSION}")
message(STATUS "HIP platform:     ${HIP_PLATFORM}")
message(STATUS "HIP compiler:     ${HIP_COMPILER}")
message(STATUS "HIP linker:       ${CMAKE_CXX_LINK_EXECUTABLE}")
message(STATUS "AMD targets:      ${AMDGPU_TARGETS}")
