# ---------------------------------------------------------------------------
# Programmer: David J. Gardner, and Cody J. Balos @ LLNL
# ---------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------
# Locate OpenMP and test for OpenMP device offloading support
# ---------------------------------------------------------------------------

find_package(OpenMP)

# Check for OpenMP offloading support
if(OPENMP_FOUND AND OPENMP_DEVICE_ENABLE)

  if(SKIP_OPENMP_DEVICE_CHECK)

    # The user has asked for checks to be skipped, assume offloading is supported
    set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)

  else()

    # If CMake version is 3.9 or newer, the FindOpenMP module checks the
    # OpenMP version.
    if((CMAKE_VERSION VERSION_EQUAL 3.9) OR (CMAKE_VERSION VERSION_GREATER 3.9))

      message(STATUS "Checking whether OpenMP supports device offloading")

      if((OpenMP_C_VERSION VERSION_EQUAL 4.5) OR (OpenMP_C_VERSION VERSION_GREATER 4.5))
        message(STATUS "Checking whether OpenMP supports device offloading -- yes")
        set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
      else()
        message(STATUS "Checking whether OpenMP supports device offloading -- no")
        set(OPENMP_SUPPORTS_DEVICE_OFFLOADING FALSE)
      endif()

    else()

      # CMake OpenMP version check not available, assume offloading is supported
      set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
      print_warning("Unable to determine OpenMP offloading support."
        "OPENMP_DEVICE_ENABLE is ON but device offloading may not function.")

    endif()

  endif()

endif()
