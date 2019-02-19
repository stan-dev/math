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
  IF(MPI_ENABLE AND MPI_C_FOUND)
    ADD_SUBDIRECTORY(parallel)
  ENDIF()
  IF(PETSC_ENABLE AND PETSC_FOUND)
    ADD_SUBDIRECTORY(petsc)
  ENDIF()
ENDIF()

# C++ examples
if(EXAMPLES_ENABLE_CXX)
  if(Trilinos_ENABLE AND Trilinos_FUNCTIONAL)
    add_subdirectory(trilinos)
  endif()
endif()

# Fortran examples
IF(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F77 AND F77_FOUND)
  ADD_SUBDIRECTORY(fcmix_serial)
  IF(OPENMP_ENABLE AND OPENMP_FOUND)
    ADD_SUBDIRECTORY(fcmix_openmp)
  ENDIF()
  IF(PTHREADS_ENABLE AND PTHREADS_FOUND)
    ADD_SUBDIRECTORY(fcmix_pthreads)
  ENDIF()
  IF(MPI_ENABLE AND MPI_Fortran_FOUND)
    ADD_SUBDIRECTORY(fcmix_parallel)
  ENDIF()
ENDIF()

# cuda examples
if(EXAMPLES_ENABLE_CUDA)

  if(CUDA_ENABLE AND CUDA_FOUND)
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

endif(EXAMPLES_ENABLE_CUDA)

