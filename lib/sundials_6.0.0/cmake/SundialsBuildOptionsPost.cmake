# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# SUNDIALS build options that are interepreted after all other
# CMake configuration.
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Option to use specialized fused kernels in the packages.
# Currently only available in CVODE.
# ---------------------------------------------------------------

sundials_option(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS BOOL "Build specialized fused GPU kernels" OFF
                DEPENDS_ON BUILD_CVODE
                DEPENDS_ON_THROW_ERROR
                SHOW_IF BUILD_CVODE)

# ---------------------------------------------------------------
# Options to enable/disable build for NVECTOR modules.
# ---------------------------------------------------------------

# required modules are in the build list, but cannot be disabled
set(BUILD_NVECTOR_SERIAL TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_SERIAL")

sundials_option(BUILD_NVECTOR_CUDA BOOL "Build the NVECTOR_CUDA module (requires CUDA)" ON
                DEPENDS_ON ENABLE_CUDA CMAKE_CUDA_COMPILER
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_CUDA")

sundials_option(BUILD_NVECTOR_HIP BOOL "Build the NVECTOR_HIP module (requires HIP)" ON
                DEPENDS_ON ENABLE_HIP
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_HIP")

sundials_option(BUILD_NVECTOR_SYCL BOOL "Build the NVECTOR_SYCL module (requires SYCL)" ON
                DEPENDS_ON ENABLE_SYCL
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_SYCL")

sundials_option(BUILD_NVECTOR_MANYVECTOR BOOL "Build the NVECTOR_MANYVECTOR module" ON
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_MANYVECTOR")

sundials_option(BUILD_NVECTOR_MPIMANYVECTOR BOOL "Build the NVECTOR_MPIMANYVECTOR module (requires MPI)" ON
                DEPENDS_ON ENABLE_MPI MPI_C_FOUND
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_MPIMANYVECTOR")

sundials_option(BUILD_NVECTOR_MPIPLUSX BOOL "Build the NVECTOR_MPIPLUSX module (requires MPI)" ON
                DEPENDS_ON ENABLE_MPI MPI_C_FOUND BUILD_NVECTOR_MPIMANYVECTOR
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_MPIPLUSX")

sundials_option(BUILD_NVECTOR_PARALLEL BOOL "Build the NVECTOR_PARALLEL module (requires MPI)" ON
                DEPENDS_ON ENABLE_MPI MPI_C_FOUND
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_PARALLEL")

sundials_option(BUILD_NVECTOR_OPENMP BOOL "Build the NVECTOR_OPENMP module" ON
                DEPENDS_ON ENABLE_OPENMP
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_OPENMP")

sundials_option(BUILD_NVECTOR_OPENMPDEV BOOL "Build the NVECTOR_OPENMPDEV module" ON
                DEPENDS_ON ENABLE_OPENMP_DEVICE OPENMP_SUPPORTS_DEVICE_OFFLOADING
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_OPENMPDEV")

sundials_option(BUILD_NVECTOR_PARHYP BOOL "Build the NVECTOR_PARHYP module (requires hypre)" ON
                DEPENDS_ON ENABLE_HYPRE HYPRE_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_PARHYP")

sundials_option(BUILD_NVECTOR_PETSC BOOL "Build the NVECTOR_PETSC module (requires PETSc)" ON
                DEPENDS_ON ENABLE_PETSC PETSC_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_PETSC")

sundials_option(BUILD_NVECTOR_PTHREADS BOOL "Build the NVECTOR_PTHREADS module" ON
                DEPENDS_ON ENABLE_PTHREAD
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_PTHREADS")

sundials_option(BUILD_NVECTOR_RAJA BOOL "Build the NVECTOR_RAJA module (requires RAJA)" ON
                DEPENDS_ON ENABLE_RAJA
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_RAJA")

sundials_option(BUILD_NVECTOR_TRILINOS BOOL "Build the NVECTOR_TRILINOS module (requires Trilinos)" ON
                DEPENDS_ON ENABLE_TRILINOS Trilinos_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_NVECTOR_TRILINOS")


# ---------------------------------------------------------------
# Options to enable/disable build for SUNMATRIX modules.
# ---------------------------------------------------------------

# required modules are in the build list, but cannot be disabled
set(BUILD_SUNMATRIX_BAND TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNMATRIX_BAND")
set(BUILD_SUNMATRIX_DENSE TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNMATRIX_DENSE")
set(BUILD_SUNMATRIX_SPARSE TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNMATRIX_SPARSE")

set(_COMPATIBLE_INDEX_SIZE FALSE)
if(SUNDIALS_INDEX_SIZE MATCHES "32")
  set(_COMPATIBLE_INDEX_SIZE TRUE)
endif()
sundials_option(BUILD_SUNMATRIX_CUSPARSE BOOL "Build the SUNMATRIX_CUSPARSE module (requires CUDA and 32-bit indexing)" ON
                DEPENDS_ON ENABLE_CUDA CMAKE_CUDA_COMPILER _COMPATIBLE_INDEX_SIZE BUILD_NVECTOR_CUDA
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNMATRIX_CUSPARSE")

sundials_option(BUILD_SUNMATRIX_MAGMADENSE BOOL "Build the SUNMATRIX_MAGMADENSE module (requires MAGMA)" ON
                DEPENDS_ON ENABLE_MAGMA MAGMA_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNMATRIX_MAGMADENSE")

sundials_option(BUILD_SUNMATRIX_ONEMKLDENSE BOOL "Build the SUNMATRIX_ONEMKLDENSE module (requires oneMKL)" ON
                DEPENDS_ON ENABLE_ONEMKL ONEMKL_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNMATRIX_ONEMKLDENSE")

sundials_option(BUILD_SUNMATRIX_SLUNRLOC BOOL "Build the SUNMATRIX_SLUNRLOC module (requires SuperLU_DIST)" ON
                DEPENDS_ON ENABLE_SUPERLUDIST SUPERLUDIST_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNMATRIX_SLUNRLOC")

# ---------------------------------------------------------------
# Options to enable/disable build for SUNLINSOL modules.
# ---------------------------------------------------------------

# required modules are in the build list, but cannot be disabled
set(BUILD_SUNLINSOL_BAND TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_BAND")
set(BUILD_SUNLINSOL_DENSE TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_DENSE")
set(BUILD_SUNLINSOL_PCG TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_PCG")
set(BUILD_SUNLINSOL_SPBCGS TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_SPBCGS")
set(BUILD_SUNLINSOL_SPFGMR TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_SPFGMR")
set(BUILD_SUNLINSOL_SPGMR TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_SPGMR")
set(BUILD_SUNLINSOL_SPTFQMR TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_SPTFQMR")

sundials_option(BUILD_SUNLINSOL_CUSOLVERSP BOOL "Build the SUNLINSOL_CUSOLVERSP module (requires CUDA and 32-bit indexing)" ON
                DEPENDS_ON ENABLE_CUDA CMAKE_CUDA_COMPILER BUILD_NVECTOR_CUDA BUILD_SUNMATRIX_CUSPARSE
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_CUSOLVERSP")

sundials_option(BUILD_SUNLINSOL_KLU BOOL "Build the SUNLINSOL_KLU module (requires KLU)" ON
                DEPENDS_ON ENABLE_KLU KLU_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_KLU")

sundials_option(BUILD_SUNLINSOL_LAPACKBAND BOOL "Build the SUNLINSOL_LAPACKBAND module (requires LAPACK)" ON
                DEPENDS_ON ENABLE_LAPACK LAPACK_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_LAPACKBAND")

sundials_option(BUILD_SUNLINSOL_LAPACKDENSE BOOL "Build the SUNLINSOL_LAPACKDENSE module (requires LAPACK)" ON
                DEPENDS_ON ENABLE_LAPACK LAPACK_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_LAPACKDENSE")

sundials_option(BUILD_SUNLINSOL_MAGMADENSE BOOL "Build the SUNLINSOL_MAGMADENSE module (requires MAGMA)" ON
                DEPENDS_ON ENABLE_MAGMA MAGMA_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_MAGMADENSE")

sundials_option(BUILD_SUNLINSOL_ONEMKLDENSE BOOL "Build the SUNLINSOL_ONEMKLDENSE module (requires oneMKL)" ON
                DEPENDS_ON ENABLE_ONEMKL ONEMKL_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_ONEMKLDENSE")

sundials_option(BUILD_SUNLINSOL_SUPERLUDIST BOOL "Build the SUNLINSOL_SUPERLUDIST module (requires SUPERLUDIST)" ON
                DEPENDS_ON ENABLE_SUPERLUDIST SUPERLUDIST_WORKS BUILD_SUNMATRIX_SLUNRLOC
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_SUPERLUDIST")

sundials_option(BUILD_SUNLINSOL_SUPERLUMT BOOL "Build the SUNLINSOL_SUPERLUMT module (requires SUPERLUMT)" ON
                DEPENDS_ON ENABLE_SUPERLUMT SUPERLUMT_WORKS
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNLINSOL_SUPERLUMT")


# ---------------------------------------------------------------
# Options to enable/disable build for SUNNONLINSOL modules.
# ---------------------------------------------------------------

# required modules are in the build list, but cannot be disabled
set(BUILD_SUNNONLINSOL_NEWTON TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNNONLINSOL_NEWTON")
set(BUILD_SUNNONLINSOL_FIXEDPOINT TRUE)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNNONLINSOL_FIXEDPOINT")

sundials_option(BUILD_SUNNONLINSOL_PETSCSNES BOOL "Build the SUNNONLINSOL_PETSCSNES module (requires PETSc)" ON
                DEPENDS_ON ENABLE_PETSC PETSC_FOUND
                ADVANCED)
list(APPEND SUNDIALS_BUILD_LIST "BUILD_SUNNONLINSOL_PETSCSNES")
