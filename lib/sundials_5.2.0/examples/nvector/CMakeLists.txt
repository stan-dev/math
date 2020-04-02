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
# examples/nvector level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Always add the nvector serial examples
add_subdirectory(serial)

if(MPI_ENABLE AND MPI_C_FOUND)
  add_subdirectory(parallel)
  add_subdirectory(mpimanyvector)
  add_subdirectory(mpiplusx)
endif()
add_subdirectory(manyvector)

if(HYPRE_ENABLE AND HYPRE_FOUND)
  add_subdirectory(parhyp)
endif()

if(PTHREAD_ENABLE AND PTHREADS_FOUND)
  add_subdirectory(pthreads)
endif()

if(OPENMP_ENABLE AND OPENMP_FOUND)
  add_subdirectory(C_openmp)
endif()

if(OPENMP_DEVICE_ENABLE AND OPENMP_FOUND AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
  add_subdirectory(openmpdev)
endif()

if(PETSC_ENABLE AND PETSC_FOUND)
  add_subdirectory(petsc)
endif()

if(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
  add_subdirectory(cuda)
  if(MPI_ENABLE AND MPI_C_FOUND)
    add_subdirectory(mpicuda)
  endif()
endif()

if(RAJA_ENABLE AND RAJA_FOUND)
  add_subdirectory(raja)
  if(MPI_ENABLE AND MPI_C_FOUND)
    add_subdirectory(mpiraja)
  endif()
endif()

if(Trilinos_ENABLE AND Trilinos_FOUND)
  add_subdirectory(trilinos)
endif()
