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
# SUNDIALS options for third-party libraries
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Enable MPI support?
# ---------------------------------------------------------------
sundials_option(ENABLE_MPI BOOL "Enable MPI support" OFF)

# ---------------------------------------------------------------
# Enable OpenMP support?
# ---------------------------------------------------------------
sundials_option(ENABLE_OPENMP BOOL "Enable OpenMP support" OFF)

# ---------------------------------------------------------------
# Enable OpenMP target offloading support?
# ---------------------------------------------------------------
sundials_option(ENABLE_OPENMP_DEVICE BOOL
                "Enable OpenMP device offloading support" OFF)

# Advanced option to skip OpenMP device offloading support check.
# This is needed for a specific compiler that doesn't correctly
# report its OpenMP spec date (with CMake >= 3.9).
sundials_option(OPENMP_DEVICE_WORKS BOOL
                "Skip the OpenMP device offloading support check" OFF
                ADVANCED)

# ---------------------------------------------------------------
# Enable Pthread support?
# ---------------------------------------------------------------
sundials_option(ENABLE_PTHREAD BOOL "Enable Pthreads support" OFF)

# -------------------------------------------------------------
# Enable CUDA support?
# -------------------------------------------------------------
sundials_option(ENABLE_CUDA BOOL "Enable CUDA support" OFF)

# CMake 3.18 adds this option.
sundials_option(CMAKE_CUDA_ARCHITECTURES STRING "Target CUDA architecture" "70"
                SHOW_IF ENABLE_CUDA)

# -------------------------------------------------------------
# Enable HIP support?
# -------------------------------------------------------------
sundials_option(ENABLE_HIP BOOL "Enable HIP support" OFF)

# -------------------------------------------------------------
# Enable SYCL support?
# -------------------------------------------------------------
sundials_option(ENABLE_SYCL BOOL "Enable SYCL support" OFF)

# ---------------------------------------------------------------
# Enable LAPACK support?
# ---------------------------------------------------------------
sundials_option(ENABLE_LAPACK BOOL "Enable Lapack support" OFF)

sundials_option(LAPACK_LIBRARIES STRING "Lapack and Blas libraries" "${LAPACK_LIBRARIES}"
                SHOW_IF ENABLE_LAPACK)

sundials_option(LAPACK_WORKS BOOL "Set to ON to force CMake to accept a given LAPACK configuration" OFF
                SHOW_IF ENABLE_LAPACK
                ADVANCED)

# ---------------------------------------------------------------
# Enable MAGMA support?
# ---------------------------------------------------------------
sundials_option(ENABLE_MAGMA BOOL "Enable MAGMA support" OFF)

sundials_option(MAGMA_DIR PATH "Path to the root of a MAGMA installation" "${MAGMA_DIR}"
                SHOW_IF ENABLE_MAGMA)

sundials_option(SUNDIALS_MAGMA_BACKENDS STRING "Which MAGMA backend under the SUNDIALS MAGMA interfaces (CUDA, HIP)" "CUDA"
                OPTIONS "CUDA;HIP"
                SHOW_IF ENABLE_MAGMA)

sundials_option(MAGMA_WORKS BOOL "Set to ON to force CMake to accept a given MAGMA configuration" OFF
                SHOW_IF ENABLE_MAGMA
                ADVANCED)

# ---------------------------------------------------------------
# Enable SuperLU_DIST support?
# ---------------------------------------------------------------
sundials_option(ENABLE_SUPERLUDIST BOOL "Enable SuperLU_DIST support" OFF)

sundials_option(SUPERLUDIST_INCLUDE_DIR PATH "SuperLU_DIST include directory" "${SUPERLUDIST_INCLUDE_DIR}"
                SHOW_IF ENABLE_SUPERLUDIST)

sundials_option(SUPERLUDIST_LIBRARY_DIR PATH "SuperLU_DIST library directory" "${SUPERLUDIST_LIBRARY_DIR}"
                SHOW_IF ENABLE_SUPERLUDIST)

sundials_option(SUPERLUDIST_LIBRARIES STRING "Semi-colon separated list of additional libraries needed for SuperLU_DIST." "${SUPERLUDIST_LIBRARIES}"
                SHOW_IF ENABLE_SUPERLUDIST)

sundials_option(SUPERLUDIST_OpenMP BOOL "Enable SUNDIALS support for SuperLU_DIST OpenMP on-node parallelism" OFF
                SHOW_IF ENABLE_SUPERLUDIST)

sundials_option(SUPERLUDIST_WORKS BOOL "Set to ON to force CMake to accept a given SUPERLUDIST configuration" OFF
                SHOW_IF ENABLE_SUPERLUDIST
                ADVANCED)

# ---------------------------------------------------------------
# Enable SuperLU_MT support?
# ---------------------------------------------------------------
sundials_option(ENABLE_SUPERLUMT BOOL "Enable SuperLU_MT support" OFF)

sundials_option(SUPERLUMT_INCLUDE_DIR PATH "SuperLU_MT include directory" "${SUPERLUMT_INCLUDE_DIR}"
                SHOW_IF ENABLE_SUPERLUMT)

sundials_option(SUPERLUMT_LIBRARY_DIR PATH "SuperLU_MT library directory" "${SUPERLUMT_LIBRARY_DIR}"
                SHOW_IF ENABLE_SUPERLUMT)

sundials_option(SUPERLUMT_LIBRARIES STRING "Semi-colon separated list of additional libraries needed for SuperLU_MT." "${SUPERLUMT_LIBRARIES}"
                SHOW_IF ENABLE_SUPERLUMT)

sundials_option(SUPERLUMT_THREAD_TYPE STRING "SuperLU_MT threading type: OPENMP or PTHREAD" "PTHREAD"
                SHOW_IF ENABLE_SUPERLUMT)

sundials_option(SUPERLUMT_WORKS BOOL "Set to ON to force CMake to accept a given SUPERLUMT configuration" OFF
                SHOW_IF ENABLE_SUPERLUMT
                ADVANCED)

# ---------------------------------------------------------------
# Enable KLU support?
# ---------------------------------------------------------------
sundials_option(ENABLE_KLU BOOL "Enable KLU support" OFF)

sundials_option(KLU_INCLUDE_DIR PATH "KLU include directory" "${KLU_INCLUDE_DIR}"
                SHOW_IF ENABLE_KLU)

sundials_option(KLU_LIBRARY_DIR PATH "KLU library directory" "${KLU_LIBRARY_DIR}"
                SHOW_IF ENABLE_KLU)

sundials_option(KLU_WORKS BOOL "Set to ON to force CMake to accept a given KLU configuration" OFF
                SHOW_IF ENABLE_KLU
                ADVANCED)

# ---------------------------------------------------------------
# Enable hypre support?
# ---------------------------------------------------------------
sundials_option(ENABLE_HYPRE BOOL "Enable hypre support" OFF)

sundials_option(HYPRE_INCLUDE_DIR PATH "HYPRE include directory" "${HYPRE_INCLUDE_DIR}"
                SHOW_IF ENABLE_HYPRE)

sundials_option(HYPRE_LIBRARY_DIR PATH "HYPRE library directory" "${HYPRE_LIBRARY_DIR}"
                SHOW_IF ENABLE_HYPRE)

sundials_option(HYPRE_WORKS BOOL "Set to ON to force CMake to accept a given hypre configuration" OFF
                SHOW_IF ENABLE_HYPRE
                ADVANCED)

# ---------------------------------------------------------------
# Enable PETSc support?
# ---------------------------------------------------------------

sundials_option(ENABLE_PETSC BOOL "Enable PETSc support" OFF)

sundials_option(PETSC_DIR PATH "Path to the root of a PETSc installation" "${PETSC_DIR}"
                SHOW_IF ENABLE_PETSC)

sundials_option(PETSC_ARCH STRING "PETSc architecture (optional)" "${PETSC_ARCH}"
                SHOW_IF ENABLE_PETSC)

sundials_option(PETSC_LIBRARIES STRING "Semi-colon separated list of PETSc link libraries" "${PETSC_LIBRARIES}"
                SHOW_IF ENABLE_PETSC
                ADVANCED)

sundials_option(PETSC_INCLUDES STRING "Semi-colon separated list of PETSc include directories" "${PETSC_INCLUDES}"
                SHOW_IF ENABLE_PETSC
                ADVANCED)

sundials_option(PETSC_WORKS BOOL "Set to ON to force CMake to accept a given PETSc configuration" OFF
                SHOW_IF ENABLE_PETSC
                ADVANCED)

# -------------------------------------------------------------
# Enable RAJA support?
# -------------------------------------------------------------
sundials_option(ENABLE_RAJA BOOL "Enable RAJA support" OFF)

sundials_option(RAJA_DIR PATH "Path to root of RAJA installation" "${RAJA_DIR}"
                SHOW_IF ENABLE_RAJA)

sundials_option(SUNDIALS_RAJA_BACKENDS STRING "Which RAJA backend under the SUNDIALS RAJA interfaces (CUDA, HIP)" "CUDA"
                OPTIONS "CUDA;HIP"
                SHOW_IF ENABLE_RAJA)

# ---------------------------------------------------------------
# Enable Trilinos support?
# ---------------------------------------------------------------
sundials_option(ENABLE_TRILINOS BOOL "Enable Trilinos support" OFF)

sundials_option(Trilinos_DIR PATH "Path to root of Trilinos installation" "${Trilinos_DIR}"
                SHOW_IF ENABLE_TRILINOS)

sundials_option(Trilinos_INTERFACE_CXX_COMPILER STRING
                "C++ compiler for Trilinos interface" "${Trilinos_CXX_COMPILER}"
                SHOW_IF ENABLE_TRILINOS
                ADVANCED)

sundials_option(Trilinos_INTERFACE_C_COMPILER STRING
                "C compiler for Trilinos interface" "${Trilinos_C_COMPILER}"
                SHOW_IF ENABLE_TRILINOS
                ADVANCED)

sundials_option(Trilinos_INTERFACE_CXX_COMPILER_FLAGS STRING
                "C++ compiler flags for Trilinos interface" "${Trilinos_CXX_COMPILER_FLAGS}"
                SHOW_IF ENABLE_TRILINOS
                ADVANCED)

sundials_option(Trilinos_INTERFACE_C_COMPILER_FLAGS STRING
                "C compiler flags for Trilinos interface" "${Trilinos_C_COMPILER_FLAGS}"
                SHOW_IF ENABLE_TRILINOS
                ADVANCED)

sundials_option(Trilinos_INTERFACE_MPIEXEC STRING
                "MPI executable for Trilinos interface" "${Trilinos_MPI_EXEC}"
                SHOW_IF ENABLE_TRILINOS
                ADVANCED)

sundials_option(Trilinos_WORKS BOOL "Set to ON to force CMake to accept a given Trilinos configuration" OFF
                SHOW_IF ENABLE_TRILINOS
                ADVANCED)

# ---------------------------------------------------------------
# Enable XBraid support?
# ---------------------------------------------------------------

sundials_option(ENABLE_XBRAID BOOL "Enable XBraid support" OFF)

sundials_option(XBRAID_DIR PATH "Path to the root of an XBraid installation" "${XBRAID_DIR}"
                DEPENDS_ON ENABLE_XBRAID)

sundials_option(XBRAID_LIBRARIES STRING "Semi-colon separated list of XBraid link libraries" "${XBRAID_LIBRARIES}"
                DEPENDS_ON ENABLE_XBRAID
                ADVANCED)

sundials_option(XBRAID_INCLUDES STRING "Semi-colon separated list of XBraid include directories" "${XBRAID_INCLUDES}"
                DEPENDS_ON ENABLE_XBRAID
                ADVANCED)

sundials_option(XBRAID_WORKS BOOL "Set to ON to force CMake to accept a given XBraid configuration" OFF
                DEPENDS_ON ENABLE_XBRAID
                ADVANCED)
