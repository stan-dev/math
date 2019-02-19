# ---------------------------------------------------------------
# Programmer(s): Daniel Reynolds @ SMU
#                David J. Gardner, Cody J. Balos @ LLNL
# -----------------------------------------------------------------
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
# CMakeLists.txt file for banded sunmatrix examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS banded matrix
SET(sunmatrix_band_examples
  "test_sunmatrix_band\;10 2 3 0\;"
  "test_sunmatrix_band\;300 7 4 0\;"
  "test_sunmatrix_band\;1000 8 8 0\;"
  "test_sunmatrix_band\;5000 3 20 0\;"
)

# Dependencies for sunmatrix examples
SET(sunmatrix_band_dependencies
  test_sunmatrix
  sundials_nvector
  sundials_matrix
)

# If building F2003 tests
if (F2003_INTERFACE_ENABLE)
  set(sunmatrix_band_fortran_examples
    "ftest_sunmatrix_band\;\;")
endif()

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(SUNMATS_LIB sundials_sunmatrixband_static)
  if (F2003_INTERFACE_ENABLE)
    LIST(APPEND NVECS_LIB sundials_fnvecserial_mod_static)
    LIST(APPEND SUNMATS_LIB sundials_fsunmatrixband_mod_static)
  endif()
ELSE()
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(SUNMATS_LIB sundials_sunmatrixband_shared)
  if (F2003_INTERFACE_ENABLE)
    LIST(APPEND NVECS_LIB sundials_fnvecserial_mod_shared)
    LIST(APPEND SUNMATS_LIB sundials_fsunmatrixband_mod_shared)
  endif()
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${NVECS_LIB} ${SUNMATS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
FOREACH(example_tuple ${sunmatrix_band_examples})

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
      ../test_sunmatrix.c ../../../src/sundials/sundials_matrix.c
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

  # install example source files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c
      ../test_sunmatrix.c
      ../test_sunmatrix.h
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/band)
  ENDIF()

ENDFOREACH(example_tuple ${sunmatrix_band_examples})

# Add the build and install targets for each example
FOREACH(example_tuple ${sunmatrix_band_fortran_examples})

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
    ADD_EXECUTABLE(${example} ${example}.f90)

    # folder to organize targets in an IDE
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SET_TARGET_PROPERTIES(${example} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

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

ENDFOREACH(example_tuple ${sunmatrix_band_fortran_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/band)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER_LIB "sundials_sunmatrixband")

  EXAMPLES2STRING(sunmatrix_band_examples EXAMPLES)
  EXAMPLES2STRING(sunmatrix_band_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunmatrix/band/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/sunmatrix/band/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/band
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunmatrix/band/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/sunmatrix/band/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/band
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
