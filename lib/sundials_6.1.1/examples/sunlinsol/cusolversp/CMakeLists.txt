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
# CMakeLists.txt file for cuSolverSp examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS cuda nvector
set(sunlinsol_cusolversp_examples
  "test_sunlinsol_cusolversp_batchqr\;100 1 0\;\;\;"
  "test_sunlinsol_cusolversp_batchqr\;10 10 0\;\;\;"
  )

# Dependencies for nvector examples
set(sunlinsol_examples_dependencies
  test_sunlinsol
  sundials_nvector
  )

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(SUNDIALS_LIBS sundials_nveccuda
                  sundials_nvecserial
                  sundials_sunmatrixdense
                  sundials_sunmatrixsparse
                  sundials_sunlinsolcusolversp
                  ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_cusolversp_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.cu
                  ../test_sunlinsol.c
                  ../../../src/sundials/sundials_linearsolver.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cu
      ../test_sunlinsol.c
      ../test_sunlinsol.h
      ../../../src/sundials/sundials_nvector.c
      ../../../src/sundials/sundials_linearsolver.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/cusolversp)
  endif()

endforeach(example_tuple ${sunlinsol_cusolversp_examples})


if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/cusolversp)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER_LIB "")
  set(NVECTOR_LIB "sundials_nvecserial sundials_nveccuda")
  set(SUNMAT_LIB "sundials_sunmatrixdense sundials_sunmatrixsparse sundials_sunmatrixcusparse")
  set(SUNLS_LIB "sundials_sunlinsolcusolversp")

  examples2string(sunlinsol_cusolversp_examples EXAMPLES_CUSOLVER)
  examples2string(sunlinsol_examples_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_CUDA_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunlinsol/cusolversp/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/cusolversp/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/cusolversp
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_CUDA_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/cusolversp/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/cusolversp/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/cusolversp
      RENAME Makefile
      )
  endif()

endif()
