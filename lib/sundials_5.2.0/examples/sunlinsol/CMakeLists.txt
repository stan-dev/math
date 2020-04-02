# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
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
# examples/sunlinsol level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Always add the serial sunlinearsolver dense/band/spils examples
add_subdirectory(band)
add_subdirectory(dense)

if(KLU_ENABLE AND KLU_FOUND)
  add_subdirectory(klu)
endif(KLU_ENABLE AND KLU_FOUND)

if(SUPERLUMT_ENABLE AND SUPERLUMT_FOUND)
  add_subdirectory(superlumt)
endif(SUPERLUMT_ENABLE AND SUPERLUMT_FOUND)

if(MPI_ENABLE AND MPI_C_FOUND AND MPI_CXX_FOUND AND SUNDIALS_SUPERLUDIST)
  add_subdirectory(superludist)
endif()

if(LAPACK_ENABLE AND LAPACK_FOUND)
  add_subdirectory(lapackband)
  add_subdirectory(lapackdense)
endif(LAPACK_ENABLE AND LAPACK_FOUND)

if(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cusolversp)
  endif()
endif(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)

# Always add serial sunlinearsolver iterative examples
add_subdirectory(spgmr/serial)
add_subdirectory(spfgmr/serial)
add_subdirectory(spbcgs/serial)
add_subdirectory(sptfqmr/serial)
add_subdirectory(pcg/serial)

if(MPI_ENABLE AND MPI_C_FOUND)
  add_subdirectory(spgmr/parallel)
  add_subdirectory(spfgmr/parallel)
  add_subdirectory(spbcgs/parallel)
  add_subdirectory(sptfqmr/parallel)
  #ADD_SUBDIRECTORY(pcg/parallel)
endif()
