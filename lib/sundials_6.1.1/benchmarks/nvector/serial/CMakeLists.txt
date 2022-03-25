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
# CMakeLists.txt file for serial nvector examples
# ---------------------------------------------------------------

message(STATUS "Added Serial NVECTOR benchmark")

sundials_add_nvector_benchmark(nvector_serial_benchmark
  SOURCES test_nvector_performance_serial.c
  SUNDIALS_TARGETS sundials_nvecserial
  INSTALL_SUBDIR benchmarks/nvector/serial
  )
