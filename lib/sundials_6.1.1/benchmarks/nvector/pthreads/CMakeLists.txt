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
# CMakeLists.txt file for PThread nvector examples
# ---------------------------------------------------------------

message(STATUS "Added PThreads NVECTOR benchmark")

sundials_add_nvector_benchmark(nvector_pthreads_benchmark
  SOURCES test_nvector_performance_pthreads.c
  SUNDIALS_TARGETS sundials_nvecpthreads
  LINK_LIBRARIES Threads::Threads
  INSTALL_SUBDIR benchmarks/nvector/pthreads
  )
