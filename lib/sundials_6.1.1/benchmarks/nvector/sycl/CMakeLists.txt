# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

message(STATUS "Added SYCL NVECTOR benchmark")

set(BENCHMARKS_DIR ${PROJECT_SOURCE_DIR}/benchmarks)

sundials_add_nvector_benchmark(test_nvector_performance_sycl
  SOURCES test_nvector_performance_sycl.cpp
  SUNDIALS_TARGETS sundials_nvecsycl
  )

install(TARGETS ${example_target}
  DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/nvector/raja")
