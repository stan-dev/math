# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# examples/nvector level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} caliper)
endif()

# Always add the nvector serial examples
add_subdirectory(serial)

# Build the nvector test utilities
add_library(test_nvector_obj OBJECT test_nvector.c)
target_link_libraries(test_nvector_obj PRIVATE sundials_nvecserial)

if(ENABLE_MPI AND MPI_C_FOUND)
  add_subdirectory(parallel)
  add_subdirectory(mpimanyvector)
  add_subdirectory(mpiplusx)
  # Build the mpi nvector test utilities
  add_library(test_nvectormpi_obj OBJECT test_mpinvector.c)
  target_link_libraries(test_nvectormpi_obj PRIVATE MPI::MPI_C sundials_nvecparallel)
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
