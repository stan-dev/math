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
# examples/ida level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# C examples
IF(EXAMPLES_ENABLE_C)
  ADD_SUBDIRECTORY(serial)
  IF(OPENMP_ENABLE AND OPENMP_FOUND)
    ADD_SUBDIRECTORY(C_openmp)
  ENDIF()
  IF(MPI_C_FOUND)
    ADD_SUBDIRECTORY(parallel)
  ENDIF()
  IF(PETSC_FOUND)
    ADD_SUBDIRECTORY(petsc)
  ENDIF()
ENDIF()

# Fortran examples
IF(EXAMPLES_ENABLE_F77 AND F77_FOUND)
  ADD_SUBDIRECTORY(fcmix_serial)
  IF(OPENMP_ENABLE AND OPENMP_FOUND)
    ADD_SUBDIRECTORY(fcmix_openmp)
  ENDIF()
  IF(PTHREADS_FOUND)
    ADD_SUBDIRECTORY(fcmix_pthreads)
  ENDIF()
  IF(MPI_Fortran_FOUND)
    ADD_SUBDIRECTORY(fcmix_parallel)
  ENDIF()
ENDIF()

# cuda examples
IF(EXAMPLES_ENABLE_CUDA)
  IF(CUDA_ENABLE AND CUDA_FOUND)
    ADD_SUBDIRECTORY(cuda)
    IF(MPI_ENABLE AND MPI_C_FOUND)
      ADD_SUBDIRECTORY(mpicuda)
    ENDIF()
  ENDIF()
ENDIF(EXAMPLES_ENABLE_CUDA)

# raja examples
IF(EXAMPLES_ENABLE_RAJA)
  IF(RAJA_ENABLE AND RAJA_FOUND)
    ADD_SUBDIRECTORY(raja)
    IF(MPI_ENABLE AND MPI_C_FOUND)
      ADD_SUBDIRECTORY(mpiraja)
    ENDIF()
  ENDIF()
ENDIF(EXAMPLES_ENABLE_RAJA)

