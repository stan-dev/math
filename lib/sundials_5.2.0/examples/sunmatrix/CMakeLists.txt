# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# examples/sunmatrix level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Always add the serial sunmatrix dense/band/sparse examples
add_subdirectory(dense)
add_subdirectory(band)
add_subdirectory(sparse)

if(MPI_ENABLE AND MPI_C_FOUND AND MPI_CXX_FOUND AND
    SUNDIALS_SUPERLUDIST)
  add_subdirectory(slunrloc)
endif()

if(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cusparse)
  endif()
endif(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
