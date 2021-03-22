# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# examples/nvector level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Always add the nvector serial examples
add_subdirectory(serial)

if(ENABLE_MPI AND MPI_C_FOUND)
  add_subdirectory(parallel)
  add_subdirectory(mpimanyvector)
  add_subdirectory(mpiplusx)
endif()
add_subdirectory(manyvector)

if(BUILD_NVECTOR_PARHYP)
  add_subdirectory(parhyp)
endif()

if(BUILD_NVECTOR_PTHREADS)
  add_subdirectory(pthreads)
endif()

if(BUILD_NVECTOR_OPENMP)
  add_subdirectory(C_openmp)
endif()

if(BUILD_NVECTOR_OPENMPDEV)
  add_subdirectory(openmpdev)
endif()

if(BUILD_NVECTOR_PETSC)
  add_subdirectory(petsc)
endif()

if(BUILD_NVECTOR_CUDA)
  add_subdirectory(cuda)
  if(ENABLE_MPI AND MPI_C_FOUND)
    add_subdirectory(mpicuda)
  endif()
endif()

if(BUILD_NVECTOR_HIP)
  add_subdirectory(hip)
endif()

if(BUILD_NVECTOR_SYCL)
  add_subdirectory(sycl)
endif()

if(BUILD_NVECTOR_RAJA)
  add_subdirectory(raja)
  if(ENABLE_MPI AND MPI_C_FOUND)
    add_subdirectory(mpiraja)
  endif()
endif()

if(BUILD_NVECTOR_TRILINOS)
  add_subdirectory(trilinos)
endif()
