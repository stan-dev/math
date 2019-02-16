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

IF(MPI_C_FOUND)
  ADD_SUBDIRECTORY(parallel)
ENDIF()

IF(HYPRE_FOUND)
  ADD_SUBDIRECTORY(parhyp)
ENDIF()

IF(PTHREADS_FOUND)
  ADD_SUBDIRECTORY(pthreads)
ENDIF()

IF(OPENMP_FOUND)
  IF(OPENMP_ENABLE)
    ADD_SUBDIRECTORY(C_openmp)
  ENDIF()
  IF(OPENMP_DEVICE_ENABLE)
    ADD_SUBDIRECTORY(openmpdev)
  ENDIF()
ENDIF()

IF(PETSC_FOUND)
  ADD_SUBDIRECTORY(petsc)
ENDIF()

IF(CUDA_FOUND)
  ADD_SUBDIRECTORY(cuda)
  IF(MPI_C_FOUND)
    ADD_SUBDIRECTORY(mpicuda)
  ENDIF()
ENDIF()

IF(RAJA_FOUND)
  ADD_SUBDIRECTORY(raja)
  IF(MPI_C_FOUND)
    ADD_SUBDIRECTORY(mpiraja)
  ENDIF()
ENDIF()
