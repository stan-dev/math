# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_ONEMKLDENSE\n\")")

# Create the library
sundials_add_library(sundials_sunmatrixonemkldense
  SOURCES
    sunmatrix_onemkldense.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_onemkldense.h
  INCLUDE_SUBDIR
    sunmatrix
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_sunmemsycl_obj
  COMPILE_FEATURES
    PUBLIC cxx_std_17
  INCLUDE_DIRECTORIES
    PUBLIC ${MKL_INCLUDE_DIR}
  LINK_LIBRARIES
    PUBLIC MKL::MKL_DPCPP
  OUTPUT_NAME
    sundials_sunmatrixonemkldense
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_ONEMKLDENSE module")
