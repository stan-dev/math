# ---------------------------------------------------------------
# Author: David J. Gardner @ LLNL
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
#
# SUNDIALS_ADD_TEST_INSTALL(<package name> <executable>)
#
# CMake macro to add a Sundials installation smoke test. Keyword
# input arguments can be added after <executable> to set test
# options (see oneValueArgs and multiValueArgs below).
# ---------------------------------------------------------------

MACRO(SUNDIALS_ADD_TEST_INSTALL PACKAGE EXECUTABLE)

  # macro options
  SET(options )

  # macro keyword inputs followed by a single value
  # EXAMPLE_DIR = path to the directory containing the installed example
  SET(oneValueArgs "EXAMPLE_DIR")

  # macro keyword inputs followed by multiple values
  SET(multiValueArgs )

  # parse inputs and create variables SUNDIALS_ADD_TEST_<keyword>
  CMAKE_PARSE_ARGUMENTS(SUNDIALS_ADD_TEST_INSTALL
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # create testing directory for this solver
  FILE(MAKE_DIRECTORY ${TEST_INSTALL_DIR}/${PACKAGE})

  ADD_CUSTOM_TARGET(${PACKAGE}_test_install
    COMMAND ${CMAKE_COMMAND} ${SUNDIALS_ADD_TEST_INSTALL_EXAMPLE_DIR} > cmake.out
    COMMAND make ${EXECUTABLE} > make.out
    COMMAND ctest
    COMMENT "Running ${PACKAGE} installation tests"
    WORKING_DIRECTORY ${TEST_INSTALL_DIR}/${PACKAGE})

  # make test_install depend on solver_test_install
  ADD_DEPENDENCIES(test_install ${PACKAGE}_test_install)

ENDMACRO()
