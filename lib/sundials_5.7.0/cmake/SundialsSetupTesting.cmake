# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Enable SUNDIALS Testing
# ---------------------------------------------------------------

# Include development examples in regression tests
sundials_option(SUNDIALS_TEST_DEVTESTS BOOL "Include development tests in make test" OFF ADVANCED)

# Include unit tests in regression tests
sundials_option(SUNDIALS_TEST_UNITTESTS BOOL "Include unit tests in make test" OFF ADVANCED)

# Enable testing with 'make test'
include(CTest)

# Add unit tests to the build if they are enabled
if(SUNDIALS_TEST_UNITTESTS)
  add_subdirectory(test/unit_tests)
endif()

# Check if development tests are enabled
if(SUNDIALS_TEST_DEVTESTS)

  message("SUNDIALS Development testing")

  # Python is needed to use the test runner
  find_package(PythonInterp REQUIRED)
  if(${PYTHON_VERSION_MAJOR} LESS 3)
    if(${PYTHON_VERSION_MINOR} LESS 7)
      print_warning("Python version must be 2.7.x or greater to run development tests"
        "Examples will build but 'make test' may fail.")
    endif()
  endif()

  # look for the testRunner script in the test directory
  find_program(TESTRUNNER testRunner PATHS test NO_DEFAULT_PATH)
  if(NOT TESTRUNNER)
    print_error("Could not locate testRunner. Set SUNDIALS_TEST_DEVTESTS=OFF to continue.")
  endif()
  message(STATUS "Found testRunner: ${TESTRUNNER}")
  set(TESTRUNNER ${TESTRUNNER} CACHE INTERNAL "")

  # Create the default test output directory
  set(TEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/Testing/output)

  if(NOT EXISTS ${TEST_OUTPUT_DIR})
    file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
  endif()

  # If a non-default output directory was provided make sure it exists
  if(SUNDIALS_TEST_OUTPUT_DIR)
    message(STATUS "Using non-default test output directory: ${SUNDIALS_TEST_OUTPUT_DIR}")
    if(NOT EXISTS ${SUNDIALS_TEST_OUTPUT_DIR})
      file(MAKE_DIRECTORY ${SUNDIALS_TEST_OUTPUT_DIR})
    endif()
  endif()

  # If a non-default answer directory was provided make sure it exists
  if(SUNDIALS_TEST_ANSWER_DIR)
    message(STATUS "Using non-default test answer directory: ${SUNDIALS_TEST_ANSWER_DIR}")
    if(NOT EXISTS ${SUNDIALS_TEST_ANSWER_DIR})
      print_error("SUNDIALS_TEST_ANSWER_DIR does not exist!")
    endif()
  endif()

  # Check if using non-default comparison precisions when testing
  if(SUNDIALS_TEST_FLOAT_PRECISION)
    message(STATUS "Using non-default float precision: ${SUNDIALS_TEST_FLOAT_PRECISION}")
  endif()

  if(SUNDIALS_TEST_INTEGER_PRECISION)
    message(STATUS "Using non-default integer precision: ${SUNDIALS_TEST_INTEGER_PRECISION}")
  endif()

endif()

# If examples are installed, create post install smoke test targets
if(EXAMPLES_INSTALL)

  # Directories for installation testing
  set(TEST_INSTALL_DIR ${PROJECT_BINARY_DIR}/Testing_Install)
  set(TEST_INSTALL_ALL_DIR ${PROJECT_BINARY_DIR}/Testing_Install_All)

  # Create installation testing directories
  if(NOT EXISTS ${TEST_INSTALL_DIR})
    file(MAKE_DIRECTORY ${TEST_INSTALL_DIR})
  endif()

  if(NOT EXISTS ${TEST_INSTALL_ALL_DIR})
    file(MAKE_DIRECTORY ${TEST_INSTALL_ALL_DIR})
  endif()

  # Create test_install and test_install_all targets
  add_custom_target(test_install
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan
    "All installation tests complete.")

  add_custom_target(test_install_all
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan
    "All installation tests complete.")

endif()