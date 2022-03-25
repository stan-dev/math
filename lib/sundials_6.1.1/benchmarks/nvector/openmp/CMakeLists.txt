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
# CMakeLists.txt file for OpenMP nvector examples
# ---------------------------------------------------------------

message(STATUS "Added OpenMP NVECTOR benchmark")

# Set-up linker flags and link libraries
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

sundials_add_nvector_benchmark(nvector_openmp_benchmark
  SOURCES test_nvector_performance_openmp.c
  SUNDIALS_TARGETS sundials_nvecopenmp
  INSTALL_SUBDIR benchmarks/nvector/openmp
  )
