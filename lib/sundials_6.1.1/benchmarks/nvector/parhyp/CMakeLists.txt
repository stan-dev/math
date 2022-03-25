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
# CMakeLists.txt file for hypre nvector examples
# ---------------------------------------------------------------

message(STATUS "Added hypre NVECTOR benchmark")

sundials_add_nvector_benchmark(nvector_parhyp_benchmark
  SOURCES test_nvector_performance_parhyp.c
  SUNDIALS_TARGETS sundials_nvecparhyp
  LINK_LIBRARIES MPI::MPI_C SUNDIALS::HYPRE
  INSTALL_SUBDIR benchmarks/nvector/parhyp
  )
