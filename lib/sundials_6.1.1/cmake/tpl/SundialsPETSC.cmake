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
# Module to find and setup PETSC correctly.
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

if(NOT DEFINED SUNDIALS_PETSC_INCLUDED)
  set(SUNDIALS_PETSC_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# Using PETSc requires building with MPI enabled
if(ENABLE_PETSC AND NOT ENABLE_MPI)
  print_error("MPI is required for PETSc support. Set ENABLE_MPI to ON.")
endif()

if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  print_error("SUNDIALS is not compatible with PETSc when using ${SUNDIALS_PRECISION} precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(PETSC REQUIRED)

message(STATUS "PETSC_DIR:         ${PETSC_DIR}")
message(STATUS "PETSC_LIBRARIES:   ${PETSC_LIBRARIES_}")
message(STATUS "PETSC_INCLUDES:    ${PETSC_INCLUDES_}")
message(STATUS "PETSC_INDEX_SIZE:  ${PETSC_INDEX_SIZE}")
message(STATUS "PETSC_PRECISION:   ${PETSC_PRECISION}\n")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(PETSC_FOUND AND (NOT PETSC_WORKS))
  # No need for any compile tests because the FindPETSC module
  # does compile tests already.

  if(NOT ("${SUNDIALS_INDEX_SIZE}" MATCHES "${PETSC_INDEX_SIZE}"))
    string(CONCAT _err_msg_string
          "PETSc not functional due to index size mismatch:\n"
          "SUNDIALS_INDEX_SIZE=${SUNDIALS_INDEX_SIZE}, "
          "but PETSc was built with ${PETSC_INDEX_SIZE}-bit indices\n"
          "PETSC_DIR: ${PETSC_DIR}\n")
    print_error("${_err_msg_string}")
  endif()

  string(TOUPPER "${PETSC_PRECISION}" _petsc_precision)
  string(TOUPPER "${SUNDIALS_PRECISION}" _sundials_precision)
  if(NOT ("${_sundials_precision}" MATCHES "${_petsc_precision}"))
    string(CONCAT _err_msg_string
          "PETSc not functional due to real type precision mismatch:\n"
          "SUNDIALS_PRECISION=${_sundials_precision}, "
          "but PETSc was built with ${_petsc_precision} precision\n"
          "PETSC_DIR: ${PETSC_DIR}\n")
    print_error("${_err_msg_string}")
  endif()

  set(PETSC_WORKS TRUE CACHE BOOL "PETSC works with SUNDIALS as configured" FORCE)
elseif(PETSC_FOUND AND PETSC_WORKS)
  message(STATUS "Skipped PETSC tests, assuming PETSC works with SUNDIALS. Set PETSC_WORKS=FALSE to (re)run compatibility test.")
endif()
