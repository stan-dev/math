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
# CMakeLists.txt file for PETSc nvector examples
# ---------------------------------------------------------------

message(STATUS "Added PETSc NVECTOR benchmark")

sundials_add_nvector_benchmark(nvector_petsc_benchmark
  SOURCES test_nvector_performance_petsc.c
  SUNDIALS_TARGETS sundials_nvecpetsc
  LINK_LIBRARIES MPI::MPI_C PUBLIC SUNDIALS::PETSC_VEC
  INSTALL_SUBDIR benchmarks/nvector/petsc
  )
