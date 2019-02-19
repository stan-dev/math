# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for CVODES parallel examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
SET(CVODES_examples
  "cvsAdvDiff_ASAp_non_p\;\;2\;2\;develop"
  "cvsAdvDiff_FSA_non_p\;-sensi stg t\;2\;2\;develop"
  "cvsAdvDiff_FSA_non_p\;-sensi sim t\;2\;2\;develop"
  "cvsAdvDiff_non_p\;\;2\;2\;develop"
  "cvsAtmDisp_ASAi_kry_bbd_p\;\;2\;8\;develop"
  "cvsDiurnal_FSA_kry_p\;-sensi stg t\;2\;4\;develop"
  "cvsDiurnal_FSA_kry_p\;-sensi sim t\;2\;4\;develop"
  "cvsDiurnal_kry_bbd_p\;\;2\;4\;develop"
  "cvsDiurnal_kry_p\;\;2\;4\;develop"
  )

IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
ELSE()
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${CVODES_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})


# Add the build and install targets for each example
FOREACH(example_tuple ${CVODES_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_args)
  LIST(GET example_tuple 2 number_of_nodes)
  LIST(GET example_tuple 3 number_of_tasks)
  LIST(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  IF(NOT TARGET ${example})
    # example source files
    ADD_EXECUTABLE(${example} ${example}.c)

    # folder to organize targets in an IDE
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

    IF(NOT MPI_C_COMPILER)
      TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
    ENDIF()
  ENDIF()

  # check if example args are provided and set the test name
  IF("${example_args}" STREQUAL "")
    SET(test_name ${example})
  ELSE()
    STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  ENDIF()

  # add example to regression tests
  SUNDIALS_ADD_TEST(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  FILE(GLOB example_out ${example}.out*)

  # install example source and .out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/parallel)
  ENDIF()

ENDFOREACH(example_tuple ${CVODES_examples})


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/parallel)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CVODES")
  SET(SOLVER_LIB "sundials_cvodes")

  EXAMPLES2STRING(CVODES_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/parallel
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/parallel
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
