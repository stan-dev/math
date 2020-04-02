# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Shelby Lockhart @ LLNL
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
# examples/cvode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C examples
if(EXAMPLES_ENABLE_C)
  add_subdirectory(serial)
  if(OPENMP_ENABLE AND OPENMP_FOUND)
    add_subdirectory(C_openmp)
  endif()
  if(OPENMP_DEVICE_ENABLE AND OPENMP_FOUND AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    add_subdirectory(C_openmpdev)
  endif()
  if(MPI_ENABLE AND MPI_C_FOUND)
    add_subdirectory(parallel)
    add_subdirectory(C_mpimanyvector)
  endif()
  if(HYPRE_ENABLE AND HYPRE_FOUND)
    add_subdirectory(parhyp)
  endif()
  if(PETSC_ENABLE AND PETSC_FOUND)
    add_subdirectory(petsc)
  endif()
endif()

# CXX examples
if(EXAMPLES_ENABLE_CXX)
  if(MPI_ENABLE AND MPI_CXX_FOUND AND SUNDIALS_SUPERLUDIST)
    add_subdirectory(superludist)
  endif()
endif()

# Fortran examples
if(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F77 AND F77_FOUND)
  add_subdirectory(fcmix_serial)
  if(MPI_ENABLE AND MPI_Fortran_FOUND)
    add_subdirectory(fcmix_parallel)
  endif()
endif()

if(F2003_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F2003 AND F2003_FOUND)
  add_subdirectory(F2003_serial)
endif()

# CUDA examples
if(EXAMPLES_ENABLE_CUDA)

  if(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
    add_subdirectory(cuda)
  endif()

  if(RAJA_ENABLE AND RAJA_FOUND)
    add_subdirectory(raja)
  endif()

endif(EXAMPLES_ENABLE_CUDA)
