# ---------------------------------------------------------------
# Programmer(s):  Eddy Banks, Cody J. Balos @ LLNL
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
# PETSc tests for SUNDIALS CMake-based configuration.
# ---------------------------------------------------------------

### This is only set if running GUI - simply return first time enabled
if(PETSC_DISABLED)
  set(PETSC_DISABLED FALSE CACHE INTERNAL "GUI - now enabled" FORCE)
  return()
endif()

if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  PRINT_ERROR("SUNDIALS is not compatible with PETSc when using ${SUNDIALS_PRECISION} precision")
endif()

message(STATUS "Checking for PETSc support... ")

# --- Find PETSc and test it --- #
find_package(PETSC REQUIRED)

# If we have the PETSC libraries, check that index size and precision match.
if(PETSC_FOUND)
  if(NOT ("${SUNDIALS_INDEX_SIZE}" MATCHES "${PETSC_INDEX_SIZE}"))
    string(CONCAT _err_msg_string
          "PETSc not functional due to index size mismatch:\n"
          "SUNDIALS_INDEX_SIZE=${SUNDIALS_INDEX_SIZE}, "
          "but PETSc was built with ${PETSC_INDEX_SIZE}-bit indices\n"
          "PETSC_DIR: ${PETSC_DIR}\n")
    PRINT_ERROR("${_err_msg_string}")
  endif()

  string(TOUPPER "${PETSC_PRECISION}" _petsc_precision)
  string(TOUPPER "${SUNDIALS_PRECISION}" _sundials_precision)
  if(NOT ("${_sundials_precision}" MATCHES "${_petsc_precision}"))
    string(CONCAT _err_msg_string
          "PETSc not functional due to real type precision mismatch:\n"
          "SUNDIALS_PRECISION=${_sundials_precision}, "
          "but PETSc was built with ${_petsc_precision} precision\n"
          "PETSC_DIR: ${PETSC_DIR}\n")
    PRINT_ERROR("${_err_msg_string}")
  endif()

  message(STATUS "Checking for PETSc support... OK")
  # sundials_config.h symbols
  set(SUNDIALS_PETSC TRUE)
else()
  message(STATUS "Checking for PETSc support... FAILED")
  # sundials_config.h symbols
  set(SUNDIALS_PETSC FALSE)
endif()
