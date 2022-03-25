# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the LapackBand SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKBAND\n\")")

# Add the library
sundials_add_library(sundials_sunlinsollapackband
  SOURCES
    sunlinsol_lapackband.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_lapackband.h
  INCLUDE_SUBDIR
    sunlinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC sundials_sunmatrixband ${LAPACK_LIBRARIES}
  OUTPUT_NAME
    sundials_sunlinsollapackband
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_VERSION}
)

message(STATUS "Added SUNLINSOL_LAPACKBAND module")
