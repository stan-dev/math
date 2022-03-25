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
# CMakeLists.txt file for the SYCL NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_SYCL\n\")")

# Create the library
sundials_add_library(sundials_nvecsycl
  SOURCES
    nvector_sycl.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_sycl.h
  INCLUDE_SUBDIR
    nvector
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_sunmemsycl_obj
  COMPILE_FEATURES
    PUBLIC cxx_std_17
  OUTPUT_NAME
    sundials_nvecsycl
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_SYCL module")
