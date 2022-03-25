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

message(STATUS "Added Parallel NVECTOR benchmark")

sundials_add_nvector_benchmark(nvector_parallel_benchmark
  SOURCES test_nvector_performance_parallel.c
  SUNDIALS_TARGETS sundials_nvecparallel
  LINK_LIBRARIES MPI::MPI_CXX
  INSTALL_SUBDIR benchmarks/nvector/parallel
  )
