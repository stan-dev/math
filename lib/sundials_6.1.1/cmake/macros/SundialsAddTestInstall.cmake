# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
#
# SUNDIALS_ADD_TEST_INSTALL(<package name> <test dir>)
#
# CMake macro to add a Sundials installation smoke tests.
# ---------------------------------------------------------------

macro(SUNDIALS_ADD_TEST_INSTALL PACKAGE TESTDIR)

  # required macro args
  # PACKAGE = Sundials package name (e.g., cvode, arkode, etc.)
  # TESTDIR = Test directory name (e.g., serial, C_parallel, etc.)

  # macro options
  set(options )

  # macro keyword inputs followed by a single value
  # EXECUTABLE = executable to add to make test_install target
  set(oneValueArgs EXECUTABLE)

  # macro keyword inputs followed by multiple values
  set(multiValueArgs )

  # parse inputs and create variables SUNDIALS_ADD_TEST_<keyword>
  cmake_parse_arguments(SUNDIALS_ADD_TEST_INSTALL
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(SUNDIALS_ADD_TEST_INSTALL_EXECUTABLE)

    # create testing directory if necessary
    if(NOT EXISTS ${TEST_INSTALL_DIR}/${PACKAGE}/${TESTDIR})
      file(MAKE_DIRECTORY ${TEST_INSTALL_DIR}/${PACKAGE}/${TESTDIR})
    endif()

    # build and run only the desired install test
    add_custom_target(test_install_${PACKAGE}_${TESTDIR}
      COMMENT "Running ${PACKAGE} installation tests"
      WORKING_DIRECTORY ${TEST_INSTALL_DIR}/${PACKAGE}/${TESTDIR}
      VERBATIM
      COMMAND ${CMAKE_COMMAND} ${EXAMPLES_INSTALL_PATH}/${PACKAGE}/${TESTDIR} > cmake.out
      COMMAND ${CMAKE_COMMAND} --build ${TEST_INSTALL_DIR}/${PACKAGE}/${TESTDIR} --target ${SUNDIALS_ADD_TEST_INSTALL_EXECUTABLE} > make.out
      COMMAND ${CMAKE_CTEST_COMMAND} -R ^${SUNDIALS_ADD_TEST_INSTALL_EXECUTABLE}$)

    # make test_install depend on test_install_package
    add_dependencies(test_install test_install_${PACKAGE}_${TESTDIR})

  endif()

  # Possible extensions:
  #  * Make EXECUTABLE a multiple value option to add several tests to test_install
  #  * Make test_install_all only available when development tests are turned on

  # create testing directory if necessary
  if(NOT EXISTS ${TEST_INSTALL_ALL_DIR}/${PACKAGE}/${TESTDIR})
    file(MAKE_DIRECTORY ${TEST_INSTALL_ALL_DIR}/${PACKAGE}/${TESTDIR})
  endif()

  # build and run all install tests
  add_custom_target(test_install_all_${PACKAGE}_${TESTDIR}
    COMMENT "Running ${PACKAGE} installation tests"
    WORKING_DIRECTORY ${TEST_INSTALL_ALL_DIR}/${PACKAGE}/${TESTDIR}
    VERBATIM
    COMMAND ${CMAKE_COMMAND} ${EXAMPLES_INSTALL_PATH}/${PACKAGE}/${TESTDIR} > cmake.out
    COMMAND ${CMAKE_COMMAND} --build ${TEST_INSTALL_ALL_DIR}/${PACKAGE}/${TESTDIR} > make.out)
  # In the future add "COMMAND ${CMAKE_CTEST_COMMAND}" here to run ctest with
  # the installed examples. Left out for now as some MPI tests require running
  # with a specific number of MPI tasks.

  # make test_install_all depend on test_install_all_package
  add_dependencies(test_install_all test_install_all_${PACKAGE}_${TESTDIR})

endmacro()
