# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------

# list of tests
set(tests )

if(BUILD_ARKODE)
  list(APPEND tests "arkode\;BENCHMARK_ODE")
endif()

if(BUILD_CVODE)
  list(APPEND tests "cvode\;BENCHMARK_ODE")
endif()

if(BUILD_IDA)
  list(APPEND tests "ida\;BENCHMARK_DAE")
endif()

# create executables
foreach(test_tuple ${tests})

  # parse the test tuple
  list(GET test_tuple 0 package)
  list(GET test_tuple 1 problem_type)

  set(sources
    ${benchmark_prefix}/main_${package}.cpp
    ${shared_sources}
    buffers.cpp
    diffusion.cpp
    solution.cpp
    utils.cpp)

  # set the target name
  set(target ${package}_diffusion_2D_mpi)

  # create executable
  add_executable(${target} ${sources})

  set_target_properties(${target} PROPERTIES FOLDER "Benchmarks")

  target_compile_definitions(${target} PRIVATE ${problem_type})

  target_include_directories(${target} PRIVATE ${benchmark_prefix})

  target_link_libraries(${target}
    PRIVATE
    sundials_${package}
    sundials_nvecparallel
    MPI::MPI_CXX)

  install(TARGETS ${target}
    DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/diffusion_2D")

endforeach()
