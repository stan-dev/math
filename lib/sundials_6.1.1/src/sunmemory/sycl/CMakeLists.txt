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

# Create library
sundials_add_library(sundials_sunmemsycl
  SOURCES
    sundials_sycl_memory.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmemory/sunmemory_sycl.h
  INCLUDE_SUBDIR
    sunmemory
  OBJECT_LIB_ONLY
)
