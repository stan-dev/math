# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# Module to find and setup oneMKL correctly.
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

if(NOT DEFINED SUNDIALS_ONEMKL_INCLUDED)
  set(SUNDIALS_ONEMKL_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# oneMKL does not support extended precision
if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(FATAL_ERROR
    "oneMKL is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# oneMKL does not support 32-bit index sizes
if(SUNDIALS_INDEX_SIZE MATCHES "32")
  message(FATAL_ERROR
    "oneMKL is not compatible with ${SUNDIALS_INDEX_SIZE}-bit indices")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# Look for CMake configuration file in oneMKL installation
find_package(MKL CONFIG
             PATHS "${ONEMKL_DIR}" "${ONEMKL_DIR}/lib/cmake/mkl"
             NO_DEFAULT_PATH
             REQUIRED)

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(MKL_FOUND AND (NOT ONEMKL_WORKS))
  message(STATUS "Checking if oneMKL works... OK")
  set(ONEMKL_WORKS TRUE CACHE BOOL "oneMKL works with SUNDIALS as configured" FORCE)
else()
  message(STATUS "Skipped oneMKL tests, assuming oneMKL works with SUNDIALS.")
endif()
