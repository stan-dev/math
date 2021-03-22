# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the LapackDense SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKDENSE\n\")")

# Add the library
sundials_add_library(sundials_sunlinsollapackdense
  SOURCES
    sunlinsol_lapackdense.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_lapackdense.h
  INCLUDE_SUBDIR
    sunlinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC sundials_sunmatrixdense ${LAPACK_LIBRARIES}
  OUTPUT_NAME
    sundials_sunlinsollapackdense
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_VERSION}
)

message(STATUS "Added SUNLINSOL_LAPACKDENSE module")

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fsunlinsollapackdense
    SOURCES
    fsunlinsol_lapackdense.c
    LINK_LIBRARIES
      PRIVATE sundials_fnvecserial sundials_fsunmatrixdense sundials_sunlinsollapackdense
    OUTPUT_NAME
      sundials_fsunlinsollapackdense
    VERSION
      ${sunlinsollib_VERSION}
    SOVERSION
      ${sunlinsollib_SOVERSION}
  )
  message(STATUS "Added SUNLINSOL_LAPACKDENSE F77 interface")
endif()
