# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# CMakeLists.txt file for sunnonlinsol fixed point examples
# ------------------------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Example programs
SET(examples
  "test_sunnonlinsol_fixedpoint\;\;"
)

# Dependencies
SET(dependencies
  sundials_nvector
  sundials_nonlinearsolver
)

SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nonlinearsolver.c
)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(SUNDIALS_LIBS sundials_nvecserial_static)
  LIST(APPEND SUNDIALS_LIBS sundials_sunnonlinsolfixedpoint_static)
ELSE()
  SET(SUNDIALS_LIBS sundials_nvecserial_shared)
  LIST(APPEND SUNDIALS_LIBS sundials_sunnonlinsolfixedpoint_shared)
ENDIF()

# Set-up linker flags and link libraries
LIST(APPEND SUNDIALS_LIBS ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
FOREACH(example_tuple ${examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_args)
  LIST(GET example_tuple 2 example_type)

  # This is used to get around DLL linkage issue since we are
  # manually including shared sources, which are normally in
  # a library that is included. If this is not set build system
  # thinks the sources are externally linked.
  IF(WIN32)
    ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)
  ENDIF()

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  IF(NOT TARGET ${example})
    # example source files
    ADD_EXECUTABLE(${example}
      ${shared_SOURCES} ${example}.c)

    # folder to organize targets in an IDE
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
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
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${shared_SOURCES} ${example}.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunnonlinsol/fixedpoint)
  ENDIF()

ENDFOREACH(example_tuple ${examples})


IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunnonlinsol/fixedpoint)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER_LIB "sundials_sunnonlinsolfixedpoint")
  SET(LIBS "${LIBS} -lsundials_sunmatrixdense -lsundials_sunlinsoldense")

  # Set the link directory for the dense sunmatrix and linear solver library 
  # The generated CMakeLists.txt does not use find_library() locate it
  SET(EXTRA_LIBS_DIR "${libdir}")

  EXAMPLES2STRING(examples EXAMPLES)
  EXAMPLES2STRING(dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunnonlinsol/fixedpoint/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/sunnonlinsol/fixedpoint/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunnonlinsol/fixedpoint
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunnonlinsol/fixedpoint/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/sunnonlinsol/fixedpoint/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunnonlinsol/fixedpoint
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
