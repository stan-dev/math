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
# CMake macro for adding executables.
# ---------------------------------------------------------------

macro(sundials_add_nvector_benchmark NAME)

  set(options )
  set(singleValueArgs )
  set(multiValueArgs SOURCES SUNDIALS_TARGETS LINK_LIBRARIES
    INSTALL_SUBDIR)

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  set(BENCHMARKS_DIR ${PROJECT_SOURCE_DIR}/benchmarks)

  add_executable(${NAME}
    ${BENCHMARKS_DIR}/nvector/test_nvector_performance.c
    ${arg_SOURCES})

  set_target_properties(${NAME} PROPERTIES FOLDER "Benchmarks")

  target_include_directories(${NAME} PRIVATE
    ${BENCHMARKS_DIR}/nvector)

  target_link_libraries(${NAME} PRIVATE
    ${arg_SUNDIALS_TARGETS} ${arg_LINK_LIBRARIES} -lm)


  if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
    target_include_directories(${NAME} PRIVATE ${caliper_INCLUDE_DIR})
    target_link_libraries(${NAME} PRIVATE caliper)
  endif()

  install(TARGETS ${NAME}
    DESTINATION "${CMAKE_INSTALL_BINDIR}/${arg_INSTALL_SUBDIR}")

endmacro(sundials_add_nvector_benchmark)
