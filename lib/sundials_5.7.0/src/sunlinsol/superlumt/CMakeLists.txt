# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU, Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the superlumt SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLUMT\n\")")

# Include OpenMP flags if SuperLU_MT is using OpenMP,
# otherwise we use threads.
if(SUPERLUMT_THREAD_TYPE STREQUAL "OPENMP")
  set(_threads OpenMP::OpenMP_C)
else()
  set(_threads ${CMAKE_THREAD_LIBS_INIT})
endif()

# Add the library
sundials_add_library(sundials_sunlinsolsuperlumt
  SOURCES
    sunlinsol_superlumt.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_superlumt.h
  INCLUDE_SUBDIR
    sunlinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC
    sundials_sunmatrixsparse
    SUNDIALS::SUPERLUMT ${_threads}
  OUTPUT_NAME
    sundials_sunlinsolsuperlumt
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_VERSION}
)

message(STATUS "Added SUNLINSOL_SUPERLUMT module")

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fsunlinsolsuperlumt
    SOURCES
      fsunlinsol_superlumt.c
    LINK_LIBRARIES
      PUBLIC sundials_sunlinsolsuperlumt
      PRIVATE sundials_fnvecserial sundials_fsunmatrixsparse
    OUTPUT_NAME
      sundials_fsunlinsolsuperlumt
    VERSION
      ${sunlinsollib_VERSION}
    SOVERSION
      ${sunlinsollib_SOVERSION}
  )
  message(STATUS "Added SUNLINSOL_SUPERLUMT F77 interface")
endif()
