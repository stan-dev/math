# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# matrix level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# Always add SUNDIALS provided matrix modules
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(sparse)

if(SUPERLUDIST_ENABLE AND SUPERLUDIST_FOUND)
  add_subdirectory(slunrloc)
endif(SUPERLUDIST_ENABLE AND SUPERLUDIST_FOUND)

if(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cusparse)
  endif()
endif(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
