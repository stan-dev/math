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

message(STATUS "Added OpenMPDEV NVECTOR benchmark")

set(NVECS_LIB sundials_nvecopenmpdev)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXTRA_LINK_LIBS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

add_executable(test_nvector_performance_openmpdev
  test_nvector_performance_openmpdev.c
  ../test_nvector_performance.c
  ../../../src/sundials/sundials_nvector.c)

# folder to organize targets in an IDE
set_target_properties(test_nvector_performance_openmp PROPERTIES FOLDER "Benchmarks")

target_include_directories(test_nvector_performance_openmpdev PRIVATE ..)

# libraries to link against
target_link_libraries(test_nvector_performance_openmpdev ${SUNDIALS_LIBS})

install(TARGETS test_nvector_performance_openmpdev
  DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/nvector/openmpdev")
