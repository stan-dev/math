# ---------------------------------------------------------------
# Programmer(s): David Gardner, Slaven Peles, and
#                Shelby Lockhart @ LLNL
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
# examples/cvode level CMakeLists.txt for SUNDIALS (for CMake)
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
  endif()
  if(HYPRE_ENABLE AND HYPRE_FOUND)
    add_subdirectory(parhyp)
  endif()
endif()

# Fortran examples
if(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F77 AND F77_FOUND)
  add_subdirectory(fcmix_serial)
  if(MPI_ENABLE AND MPI_Fortran_FOUND)
    add_subdirectory(fcmix_parallel)
  endif()
endif()

if(F2003_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F90 AND F90_FOUND)
  add_subdirectory(F90_serial)
endif()

# cuda examples
if(EXAMPLES_ENABLE_CUDA)

  if(CUDA_ENABLE AND CUDA_FOUND)
    add_subdirectory(cuda)
  endif()

  if(RAJA_ENABLE AND RAJA_FOUND)
    add_subdirectory(raja)
  endif()

endif(EXAMPLES_ENABLE_CUDA)

