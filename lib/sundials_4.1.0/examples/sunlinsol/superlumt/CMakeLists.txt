# ---------------------------------------------------------------
# Programmer(s): Daniel Reynolds @ SMU
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
# CMakeLists.txt file for sunlinsol SuperLU_MT examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using the SuperLU_MT linear solver
SET(sunlinsol_superlumt_examples
  "test_sunlinsol_superlumt\;300 0 1 0\;"
  "test_sunlinsol_superlumt\;300 1 1 0\;"
  "test_sunlinsol_superlumt\;1000 0 3 0\;develop"
  "test_sunlinsol_superlumt\;1000 1 3 0\;develop"
)

# Dependencies for nvector examples
SET(sunlinsol_superlumt_dependencies
  test_sunlinsol
  sundials_nvector
  sundials_matrix
  sundials_linearsolver
)

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ..)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(SUNMATS_LIB sundials_sunmatrixdense_static)
  LIST(APPEND SUNMATS_LIB sundials_sunmatrixsparse_static)
  SET(SUNLINSOLS_LIB sundials_sunlinsolsuperlumt_static)
ELSE()
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(SUNMATS_LIB sundials_sunmatrixdense_shared)
  LIST(APPEND SUNMATS_LIB sundials_sunmatrixsparse_shared)
  SET(SUNLINSOLS_LIB sundials_sunlinsolsuperlumt_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${NVECS_LIB} ${SUNMATS_LIB} ${SUNLINSOLS_LIB} ${EXTRA_LINK_LIBS})

LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})

IF(LAPACK_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

IF(BLAS_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${BLAS_LIBRARIES})
ENDIF(BLAS_FOUND)


# Add the build and install targets for each example
FOREACH(example_tuple ${sunlinsol_superlumt_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_args)
  LIST(GET example_tuple 2 example_type)

  # This is used to get around DLL linkage issue since we are
  # manually including sundials_nvector.c here, which is normally in
  # a library that is included.  If this is not set build system
  # thinks nvector is externally linked.
  IF(WIN32)
    ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)
  ENDIF(WIN32)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  IF(NOT TARGET ${example})
    # example source files
    ADD_EXECUTABLE(${example} ${example}.c
      ../test_sunlinsol.c
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_linearsolver.c
      ../../../src/sundials/sundials_nvector.c)

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
    INSTALL(FILES ${example}.c
      ../test_sunlinsol.h
      ../test_sunlinsol.c
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_linearsolver.c
      ../../../src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superlumt)
  ENDIF()

ENDFOREACH(example_tuple ${sunlinsol_superlumt_examples})


IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superlumt)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER_LIB "sundials_sunlinsolsuperlumt")
  SET(LIBS "${LIBS} -lsundials_sunmatrixdense -lsundials_sunmatrixsparse")

  # Set the link directory for the dense and sparse sunmatrix libraries
  # The generated CMakeLists.txt does not use find_library() locate them
  SET(EXTRA_LIBS_DIR "${libdir}")

  EXAMPLES2STRING(sunlinsol_superlumt_examples EXAMPLES_SLUMT)
  EXAMPLES2STRING(sunlinsol_superlumt_dependencies EXAMPLES_DEPENDENCIES_SLUMT)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunlinsol/superlumt/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/superlumt/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superlumt
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/superlumt/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/superlumt/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superlumt
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
