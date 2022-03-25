# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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

#
# Deprecated Fortran options
#

if(DEFINED F2003_INTERFACE_ENABLE)
  print_warning("The CMake option F2003_INTERFACE_ENABLE is deprecated"
                "Use BUILD_FORTRAN_MODULE_INTERFACE instead"
                MODE DEPRECATION)
  set(BUILD_FORTRAN_MODULE_INTERFACE ${F2003_INTERFACE_ENABLE} CACHE BOOL "Enable Fortran 2003 module interfaces")
endif()

unset(F2003_INTERFACE_ENABLE CACHE)

#
# Deprecated TPL options
#

if(DEFINED MPI_ENABLE)
  print_warning("The CMake option MPI_ENABLE is deprecated" "Use ENABLE_MPI instead"
                MODE DEPRECATION)
  set(ENABLE_MPI ${MPI_ENABLE} CACHE BOOL "Enable MPI support" FORCE)
  unset(MPI_ENABLE CACHE)
endif()

if(DEFINED OPENMP_ENABLE)
  print_warning("The CMake option OPENMP_ENABLE is deprecated" "Use ENABLE_OPENMP instead"
                MODE DEPRECATION)
  set(ENABLE_OPENMP ${OPENMP_ENABLE} CACHE BOOL "Enable OpenMP support" FORCE)
  unset(OPENMP_ENABLE CACHE)
endif()

if(DEFINED OPENMP_DEVICE_ENABLE)
  print_warning("The CMake option OPENMP_DEVICE_ENABLE is deprecated"
                "Use ENABLE_OPENMP_DEVICE instead"
                MODE DEPRECATION)
  set(ENABLE_OPENMP_DEVICE ${OPENMP_DEVICE_ENABLE} CACHE BOOL
      "Enable OpenMP device offloading support" FORCE)
  unset(OPENMP_DEVICE_ENABLE CACHE)
endif()

if(DEFINED SKIP_OPENMP_DEVICE_CHECK)
  print_warning("The CMake option SKIP_OPENMP_DEVICE_CHECK is deprecated"
                "Use OPENMP_DEVICE_WORKS instead"
                MODE DEPRECATION)
  set(OPENMP_DEVICE_WORKS ${SKIP_OPENMP_DEVICE_CHECK} CACHE BOOL
     "Skip the compiler check for OpenMP device offloading" FORCE)
  unset(SKIP_OPENMP_DEVICE_CHECK CACHE)
endif()

if(DEFINED PTHREAD_ENABLE)
  print_warning("The CMake option PTHREAD_ENABLE is deprecated" "Use ENABLE_PTHREAD instead"
                MODE DEPRECATION)
  set(ENABLE_PTHREAD ${PTHREAD_ENABLE} CACHE BOOL "Enable Pthreads support" FORCE)
  unset(PTHREAD_ENABLE CACHE)
endif()

if(DEFINED CUDA_ENABLE)
  print_warning("The CMake option CUDA_ENABLE is deprecated" "Use ENABLE_CUDA instead"
                MODE DEPRECATION)
  set(ENABLE_CUDA ${CUDA_ENABLE} CACHE BOOL "Enable CUDA support" FORCE)
  unset(CUDA_ENABLE CACHE)
endif()

if(DEFINED LAPACK_ENABLE)
  print_warning("The CMake option LAPACK_ENABLE is deprecated" "Use ENABLE_LAPACK instead"
                MODE DEPRECATION)
  set(ENABLE_LAPACK ${LAPACK_ENABLE} CACHE BOOL "Enable LAPACK support" FORCE)
  unset(LAPACK_ENABLE CACHE)
endif()

if(DEFINED SUPERLUDIST_ENABLE)
  print_warning("The CMake option SUPERLUDIST_ENABLE is deprecated"
                "Use ENABLE_SUPERLUDIST instead"
                MODE DEPRECATION)
  set(ENABLE_SUPERLUDIST ${SUPERLUDIST_ENABLE} CACHE BOOL "Enable SuperLU_DIST support" FORCE)
  unset(SUPERLUDIST_ENABLE CACHE)
endif()

if(DEFINED SUPERLUMT_ENABLE)
  print_warning("The CMake option SUPERLUMT_ENABLE is deprecated" "Use ENABLE_SUPERLUMT instead"
                MODE DEPRECATION)
  set(ENABLE_SUPERLUMT ${SUPERLUMT_ENABLE} CACHE BOOL "Enable SuperLU_MT support" FORCE)
  unset(SUPERLUMT_ENABLE CACHE)
endif()

if(DEFINED KLU_ENABLE)
  print_warning("The CMake option KLU_ENABLE is deprecated" "Use ENABLE_KLU instead"
                MODE DEPRECATION)
  set(ENABLE_KLU ${KLU_ENABLE} CACHE BOOL "Enable KLU support" FORCE)
  unset(KLU_ENABLE CACHE)
endif()

if(DEFINED HYPRE_ENABLE)
  print_warning("The CMake option HYPRE_ENABLE is deprecated" "Use ENABLE_HYPRE instead"
                MODE DEPRECATION)
  set(ENABLE_HYPRE ${HYPRE_ENABLE} CACHE BOOL "Enable HYPRE support" FORCE)
  unset(HYPRE_ENABLE CACHE)
endif()

if(DEFINED PETSC_ENABLE)
  print_warning("The CMake option PETSC_ENABLE is deprecated" "Use ENABLE_PETSC instead"
                MODE DEPRECATION)
  set(ENABLE_PETSC ${PETSC_ENABLE} CACHE BOOL "Enable PETSC support" FORCE)
  unset(PETSC_ENABLE CACHE)
endif()

if(DEFINED Trilinos_ENABLE)
  print_warning("The CMake option Trilinos_ENABLE is deprecated" "Use ENABLE_TRILINOS instead"
                MODE DEPRECATION)
  set(ENABLE_TRILINOS ${Trilinos_ENABLE} CACHE BOOL "Enable Trilinos support" FORCE)
  unset(Trilinos_ENABLE CACHE)
endif()

if(DEFINED RAJA_ENABLE)
  print_warning("The CMake option RAJA_ENABLE is deprecated" "Use ENABLE_RAJA instead"
                MODE DEPRECATION)
  set(ENABLE_RAJA ${RAJA_ENABLE} CACHE BOOL "Enable RAJA support" FORCE)
  unset(RAJA_ENABLE CACHE)
endif()

#
# Deprecated CUDA_ARCH option
#

if(DEFINED CUDA_ARCH)
  print_warning("The CMake option CUDA_ARCH is deprecated" "Use CMAKE_CUDA_ARCHITECTURES instead"
                MODE DEPRECATION)
  # convert sm_** to just **
  string(REGEX MATCH "[0-9]+" arch_name "${CUDA_ARCH}")
  set(CMAKE_CUDA_ARCHITECTURES ${arch_name} CACHE STRING "CUDA Architectures" FORCE)
  unset(CUDA_ARCH)
endif()
