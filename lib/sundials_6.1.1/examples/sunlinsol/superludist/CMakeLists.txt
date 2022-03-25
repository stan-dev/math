# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# -----------------------------------------------------------------------------
# CMakeLists.txt file for the superlu-dist SUNLinearSolver examples
# -----------------------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;type" where the
# 'develop' for examples excluded from 'make test' in releases

# Examples using the SuperLU_MT linear solver
set(sunlinsol_superludist_examples
  "test_sunlinsol_superludist\;100 1 1 0\;1\;"
  "test_sunlinsol_superludist\;400 1 2 0\;2\;"
  "test_sunlinsol_superludist\;900 3 1 0\;3\;"
  "test_sunlinsol_superludist\;1000 2 2 0\;4\;"
)

# Dependencies for nvector examples
set(sunlinsol_superludist_dependencies
  test_sunlinsol
  sundials_nvector
  sundials_matrix
  sundials_linearsolver
)

if(MPI_C_COMPILER AND MPI_CXX_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

# Add source directory to include directories
include_directories(. ..)

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_superludist_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.cpp ../test_sunlinsol.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example}
      PRIVATE
        sundials_sunlinsolsuperludist
        sundials_sunmatrixslunrloc
        sundials_sunmatrixdense
        sundials_nvecparallel
        sundials_nvecserial)

  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_nodes}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_nodes})
    endif()
  else()
    if("${number_of_nodes}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_nodes}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_nodes}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cpp
      ../test_sunlinsol.h
      ../test_sunlinsol.c
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_linearsolver.c
      ../../../src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superludist)
  endif()

endforeach(example_tuple ${sunlinsol_superludist_examples})


if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superludist)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "SUNMatrix_SLUNRloc")
  string(CONCAT SUNDIALS_LIBS
    "-lsundials_nvecparallel "
    "-lsundials_sunmatrixslunrloc "
    "-lsundials_sunlinsolsuperludist")

  examples2string(sunlinsol_superludist_examples EXAMPLES)
  examples2string(sunlinsol_superludist_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_superludist_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunlinsol/superludist/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/superludist/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superludist
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_superludist_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/superludist/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/superludist/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superludist
      RENAME Makefile
      )
  endif()

endif()
