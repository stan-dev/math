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

message(STATUS "Added MPIPlusX NVECTOR benchmark")

sundials_add_nvector_benchmark(nvector_mpiplusx_benchmark
  SOURCES test_nvector_performance_mpiplusx.c
  SUNDIALS_TARGETS sundials_nvecserial sundials_nvecmpiplusx
  LINK_LIBRARIES MPI::MPI_CXX
  INSTALL_SUBDIR benchmarks/nvector/mpiplusx
  )
