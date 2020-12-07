# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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

# Create a library out of the generic sundials modules
sundials_add_library(sundials_sunmemcuda
  SOURCES
    sundials_cuda_memory.cu
  OBJECT_LIB_ONLY
)

# Install the header files
install(FILES ${SUNDIALS_SOURCE_DIR}/include/sunmemory/sunmemory_cuda.h
  DESTINATION include/sunmemory)