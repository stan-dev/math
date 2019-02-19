# ---------------------------------------------------------------
# Programmer:  David Gardner, Slaven Peles, and Cody Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# examples/nvector level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# Always add the nvector serial examples
ADD_SUBDIRECTORY(serial)

IF(MPI_ENABLE AND MPI_C_FOUND)
  ADD_SUBDIRECTORY(parallel)
ENDIF()

IF(HYPRE_ENABLE AND HYPRE_FOUND)
  ADD_SUBDIRECTORY(parhyp)
ENDIF()

IF(PTHREADS_ENABLE AND PTHREADS_FOUND)
  ADD_SUBDIRECTORY(pthreads)
ENDIF()

IF(OPENMP_ENABLE AND OPENMP_FOUND)
  ADD_SUBDIRECTORY(C_openmp)
ENDIF()

IF(OPENMP_DEVICE_ENABLE AND OPENMP_FOUND AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
  ADD_SUBDIRECTORY(openmpdev)
ENDIF()

IF(PETSC_ENABLE AND PETSC_FOUND)
  ADD_SUBDIRECTORY(petsc)
ENDIF()

IF(CUDA_ENABLE AND CUDA_FOUND)
  ADD_SUBDIRECTORY(cuda)
  IF(MPI_ENABLE AND MPI_C_FOUND)
    ADD_SUBDIRECTORY(mpicuda)
  ENDIF()
ENDIF()

IF(RAJA_ENABLE AND RAJA_FOUND)
  ADD_SUBDIRECTORY(raja)
  IF(MPI_ENABLE AND MPI_C_FOUND)
    ADD_SUBDIRECTORY(mpiraja)
  ENDIF()
ENDIF()

IF(Trilinos_ENABLE AND Trilinos_FUNCTIONAL)
  ADD_SUBDIRECTORY(trilinos)
ENDIF()
