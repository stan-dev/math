# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
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
# CMakeLists.txt file for the superlu-dist SUNLinearSolver library
# -----------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLUDIST\n\")")

# Use MPI
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Add the library
sundials_add_library(sundials_sunlinsolsuperludist
  SOURCES
    sunlinsol_superludist.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_superludist.h
  INCLUDE_SUBDIR
    sunlinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC
    sundials_sunmatrixslunrloc
    SUNDIALS::SUPERLUDIST
    $<IF:$<BOOL:${SUPERLUDIST_OpenMP}>,OpenMP::OpenMP_C,>
  COMPILE_OPTIONS
    PRIVATE ${_compile_options}
  OUTPUT_NAME
    sundials_sunlinsolsuperludist
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_VERSION}
)

message(STATUS "Added SUNLINSOL_SUPERLUDIST module")
