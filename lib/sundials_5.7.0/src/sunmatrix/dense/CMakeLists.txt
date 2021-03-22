# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the dense SUNMatrix library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_DENSE\n\")")

# Add the sunmatrix_dense library
sundials_add_library(sundials_sunmatrixdense
  SOURCES
    sunmatrix_dense.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_dense.h
  INCLUDE_SUBDIR
    sunmatrix
  OBJECT_LIBRARIES
    sundials_generic_obj
  OUTPUT_NAME
    sundials_sunmatrixdense
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_DENSE module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fsunmatrixdense
    SOURCES
      fsunmatrix_dense.c
    LINK_LIBRARIES
      PRIVATE sundials_fnvecserial sundials_sunmatrixdense
    OUTPUT_NAME
      sundials_fsunmatrixdense
    VERSION
      ${sunmatrixlib_VERSION}
    SOVERSION
      ${sunmatrixlib_SOVERSION}
  )
  message(STATUS "Added SUNMATRIX_DENSE F77 interface")
endif()
