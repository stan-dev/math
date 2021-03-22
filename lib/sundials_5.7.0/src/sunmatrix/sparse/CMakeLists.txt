# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU,
#                and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the sparse SUNMatrix library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_SPARSE\n\")")

# Add the sunmatrix_sparse library
sundials_add_library(sundials_sunmatrixsparse
  SOURCES
    sunmatrix_sparse.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_sparse.h
  INCLUDE_SUBDIR
    sunmatrix
  OBJECT_LIBRARIES
    sundials_generic_obj
  OUTPUT_NAME
    sundials_sunmatrixsparse
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_SPARSE module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fsunmatrixsparse
    SOURCES
      fsunmatrix_sparse.c
    LINK_LIBRARIES
      PRIVATE sundials_fnvecserial sundials_sunmatrixsparse
    OUTPUT_NAME
      sundials_fsunmatrixsparse
    VERSION
      ${sunmatrixlib_VERSION}
    SOVERSION
      ${sunmatrixlib_SOVERSION}
  )
  message(STATUS "Added SUNMATRIX_SPARSE F77 interface")
endif()
