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

message(STATUS "Added RAJA NVECTOR benchmark")

set(example test_nvector_performance_raja)

set(BENCHMARKS_DIR ${PROJECT_SOURCE_DIR}/benchmarks)

foreach(backend ${SUNDIALS_RAJA_BACKENDS})

  if(backend MATCHES "CUDA")
    set(_lang CUDA)
    set(_defines USE_CUDA)
    set(example_target "${example}_cuda")
    set_source_files_properties(${example}.cpp PROPERTIES LANGUAGE CUDA)
  elseif(backend MATCHES "HIP")
    set(_lang CXX)
    set(_defines USE_HIP)
    set(example_target "${example}_hip")
  elseif(backend MATCHES "SYCL")
    set(_lang CXX)
    set(_defines USE_SYCL)
    set(example_target "${example}_sycl")
  else()
    continue()
  endif()

  sundials_add_nvector_benchmark(${example_target}
    SOURCES test_nvector_performance_raja.cpp
    SUNDIALS_TARGETS sundials_nvecraja
    )

  target_compile_definitions(${example_target}
    PRIVATE ${_defines})

  if(backend MATCHES "CUDA")
    set_target_properties(${example_target} PROPERTIES LINKER_LANGUAGE CXX)
  endif()

  if(RAJA_BACKENDS MATCHES "OPENMP")
    target_link_libraries(${example_target} PRIVATE OpenMP::OpenMP_CXX)
  endif()

  install(TARGETS ${example_target}
    DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/nvector/raja")

endforeach()
