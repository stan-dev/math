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
# Module to find and setup OpenMP correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# TPL_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
#
# Creates the variables:
#   OPENMP_FOUND - was OpenMP found
#   OPENMP45_FOUND - was OpenMP v4.5 or greater found
#   OPENMP_SUPPORTS_DEVICE_OFFLOADING - is device offloading supported
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_OPENMP_INCLUDED)
  set(SUNDIALS_OPENMP_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

set(OPENMP45_FOUND FALSE)
set(OPENMP_SUPPORTS_DEVICE_OFFLOADING FALSE)
find_package(OpenMP REQUIRED)

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# Work around a bug in setting OpenMP version variables in CMake >= 3.9. The
# OpenMP version information is not stored in cache variables and is not set
# on repeated calls to find OpenMP (i.e., when using ccmake). To ensure these
# variables exist store copies of the values.
set(OpenMP_C_VERSION "${OpenMP_C_VERSION}" CACHE INTERNAL "" FORCE)
set(OpenMP_CXX_VERSION "${OpenMP_CXX_VERSION}" CACHE INTERNAL "" FORCE)
set(OpenMP_Fortran_VERSION "${OpenMP_Fortran_VERSION}" CACHE INTERNAL "" FORCE)

# Check for OpenMP offloading support
if(OPENMP_FOUND AND (ENABLE_OPENMP_DEVICE OR SUPERLUDIST_OpenMP))

  if(OPENMP_DEVICE_WORKS)

    # The user has asked for checks to be skipped, assume offloading is supported
    set(OPENMP45_FOUND TRUE)
    set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
    print_warning("Skipping OpenMP device/version check." "SUNDIALS OpenMP functionality dependent on OpenMP 4.5+ is not guaranteed.")

  else()

    # Check the OpenMP version
    message(STATUS "Checking whether OpenMP supports device offloading")

    if((OpenMP_C_VERSION VERSION_EQUAL 4.5) OR (OpenMP_C_VERSION VERSION_GREATER 4.5))
      message(STATUS "Checking whether OpenMP supports device offloading -- yes")
      set(OPENMP45_FOUND TRUE)
      set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
    else()
      message(STATUS "Checking whether OpenMP supports device offloading -- no")
      set(OPENMP45_FOUND FALSE)
      set(OPENMP_SUPPORTS_DEVICE_OFFLOADING FALSE)
      print_error("The found OpenMP version does not support device offloading.")
    endif()

  endif()

endif()