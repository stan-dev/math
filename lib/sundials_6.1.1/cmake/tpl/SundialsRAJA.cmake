# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup RAJA correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# TPL_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_RAJA_INCLUDED)
  set(SUNDIALS_RAJA_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

if((SUNDIALS_RAJA_BACKENDS MATCHES "CUDA") AND (NOT ENABLE_CUDA))
  message(FATAL_ERROR "RAJA with a CUDA backend requires ENABLE_CUDA = ON")
endif()

if((SUNDIALS_RAJA_BACKENDS MATCHES "HIP") AND (NOT ENABLE_HIP))
  message(FATAL_ERROR "RAJA with a HIP backend requires ENABLE_HIP = ON")
endif()

if((SUNDIALS_RAJA_BACKENDS MATCHES "SYCL") AND (NOT ENABLE_SYCL))
  message(FATAL_ERROR "RAJA with a SYCL backend requires ENABLE_SYCL = ON")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# find the library configuration file
find_file(RAJA_CONFIGHPP_PATH config.hpp
          HINTS "${RAJA_DIR}"
          PATH_SUFFIXES include include/RAJA
          NO_DEFAULT_PATH)
mark_as_advanced(FORCE RAJA_CONFIGHPP_PATH)

# Look for CMake configuration file in RAJA installation
find_package(RAJA CONFIG
             PATHS "${RAJA_DIR}" "${RAJA_DIR}/share/raja/cmake"
             NO_DEFAULT_PATH
             REQUIRED)

# determine the backends
foreach(_backend CUDA HIP OPENMP TARGET_OPENMP SYCL)
  file(STRINGS "${RAJA_CONFIGHPP_PATH}" _raja_has_backend REGEX "^#define RAJA_ENABLE_${_backend}\$")
  if(_raja_has_backend)
    set(RAJA_BACKENDS "${_backend};${RAJA_BACKENDS}")
  endif()
endforeach()

message(STATUS "RAJA Version:  ${RAJA_VERSION_MAJOR}.${RAJA_VERSION_MINOR}.${RAJA_VERSION_PATCHLEVEL}")
message(STATUS "RAJA Backends: ${RAJA_BACKENDS}")

# Check if RAJA uses the Threads::Threads target. Currently this target gets
# created by find_package(CUDA). The example CMake files will need to call
# find_package(Threads) to create this target.
set(RAJA_NEEDS_THREADS OFF)
set(_raja_target_list RAJA RAJA::cuda)
foreach(_raja_target ${_raja_target_list})
  if(TARGET ${_raja_target})
    get_target_property(_raja_target_interface_link_libraries
      ${_raja_target} INTERFACE_LINK_LIBRARIES)
    if(_raja_target_interface_link_libraries MATCHES "Threads::Threads")
      set(RAJA_NEEDS_THREADS ON)
      if(NOT TARGET Threads::Threads)
        find_package(Threads)
      endif()
      break()
    endif()
  endif()
endforeach()

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if((SUNDIALS_RAJA_BACKENDS MATCHES "CUDA") AND
   (NOT RAJA_BACKENDS MATCHES "CUDA"))
  print_error("Requested that SUNDIALS uses the CUDA RAJA backend, but RAJA was not built with the CUDA backend.")
endif()

if((SUNDIALS_RAJA_BACKENDS MATCHES "HIP") AND
   (NOT RAJA_BACKENDS MATCHES "HIP"))
  print_error("Requested that SUNDIALS uses the HIP RAJA backend, but RAJA was not built with the HIP backend.")
endif()

if(NOT ENABLE_OPENMP AND RAJA_BACKENDS MATCHES "OPENMP")
  print_error("RAJA was built with OpenMP, but OpenMP is not enabled. Set ENABLE_OPENMP to ON.")
endif()

if(NOT ENABLE_OPENMP_DEVICE AND RAJA_BACKENDS MATCHES "TARGET_OPENMP")
  print_error("RAJA was built with OpenMP device offloading, but OpenMP with device offloading is not enabled. Set ENABLE_OPENMP_DEVICE to ON.")
endif()

if((SUNDIALS_RAJA_BACKENDS MATCHES "SYCL") AND
    (NOT RAJA_BACKENDS MATCHES "SYCL"))
  print_error("Requested that SUNDIALS uses the SYCL RAJA backend, but RAJA was not built with the SYCL backend.")
endif()
