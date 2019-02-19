# ---------------------------------------------------------------
# Author: Steven Smith and David J. Gardner @ LLNL
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
# SUNDIALS_ADD_TEST(<test name> <executable>)
#
# CMake macro to add a Sundials regression test. Keyword input
# arguments can be added after <executable> to set regression
# test options (see oneValueArgs and multiValueArgs below).
#
# When SUNDIALS_DEVTEST is OFF (default) the executable is run
# and pass/fail is determined by the executable return value.
#
# When SUNDIALS_DEVTESTS is ON the executable is run and its
# output compared with the corresponding .out file. If the output
# differs significantly then the test fails. The default level of
# signicance is 4 decimal points for floating values and 10% for
# integer values.
#
# The level of precision can be adjusted for all tests using:
#  -D SUNDIALS_DEVTESTS_FLOAT_PRECISION=<number of digits>
#  -D SUNDIALS_DEVTESTS_INTEGER_PRECISION=<% difference>
# ---------------------------------------------------------------

MACRO(SUNDIALS_ADD_TEST NAME EXECUTABLE)

  # macro options
  # NODIFF = do not diff the test output against an answer file
  SET(options "NODIFF")

  # macro keyword inputs followed by a single value
  # MPI_NPROCS         = number of mpi tasks to use in parallel tests
  # FLOAT_PRECISION    = precision for floating point failure comparision (num digits)
  # INTEGER_PRECENTAGE = integer percentage difference for failure comparison
  # ANSWER_DIR         = path to the directory containing the test answer file
  # ANSWER_FILE        = name of test answer file
  # EXAMPLE_TYPE       = release or develop examples
  SET(oneValueArgs "MPI_NPROCS" "FLOAT_PRECISION" "INTEGER_PERCENTAGE"
    "ANSWER_DIR" "ANSWER_FILE" "EXAMPLE_TYPE")

  # macro keyword inputs followed by multiple values
  # TEST_ARGS = command line arguments to pass to the test executable
  SET(multiValueArgs "TEST_ARGS")

  # parse inputs and create variables SUNDIALS_ADD_TEST_<keyword>
  CMAKE_PARSE_ARGUMENTS(SUNDIALS_ADD_TEST
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # SGS add check to make sure parallel is integer
  # SGS add check for float and integer precision

  # check that the test is not excluded
  if(NOT ("${SUNDIALS_ADD_TEST_EXAMPLE_TYPE}" STREQUAL "exclude"))

    # add test using the dev test runner
    IF(SUNDIALS_DEVTESTS)

      # command line arguments for the test runner script
      SET(TEST_ARGS
        "--verbose"
        "--testname=${NAME}"
        "--executablename=$<TARGET_FILE:${EXECUTABLE}>"
        "--outputdir=${CMAKE_BINARY_DIR}/Testing/output"
        )

      # do not diff the output and answer files
      IF(SUNDIALS_ADD_TEST_NODIFF)
        LIST(APPEND TEST_ARGS "--nodiff")
      ENDIF()

      # check if this test is run with MPI and set the MPI run command
      IF(NOT ("${SUNDIALS_ADD_TEST_MPI_NPROCS}" STREQUAL "") AND
          NOT ("${MPIEXEC_EXECUTABLE}" STREQUAL ""))

        IF(MPIEXEC_EXECUTABLE MATCHES "srun")
          SET(RUN_COMMAND "srun -N1 -n${SUNDIALS_ADD_TEST_MPI_NPROCS} -ppdebug")
        ELSE()
          SET(RUN_COMMAND "${MPIEXEC_EXECUTABLE} -n ${SUNDIALS_ADD_TEST_MPI_NPROCS}")
        ENDIF()

        LIST(APPEND TEST_ARGS "--runcommand=\"${RUN_COMMAND}\"")

      ENDIF()

      # set the test input args
      IF(NOT ("${SUNDIALS_ADD_TEST_TEST_ARGS}" STREQUAL ""))
        STRING (REPLACE ";" " " USER_ARGS "${SUNDIALS_ADD_TEST_TEST_ARGS}")
        LIST(APPEND TEST_ARGS "--runargs=\"${USER_ARGS}\"")
      ENDIF()

      # set the test answer directory name (default is test/answers)
      IF(NOT ("${SUNDIALS_ADD_TEST_ANSWER_DIR}" STREQUAL ""))
        LIST(APPEND TEST_ARGS "--answerdir=${SUNDIALS_ADD_TEST_ANSWER_DIR}")
      ENDIF()

      # set the test answer file name (default is test_name_test_agrs)
      IF(NOT ("${SUNDIALS_ADD_TEST_ANSWER_FILE}" STREQUAL ""))
        LIST(APPEND TEST_ARGS "--answerfile=${SUNDIALS_ADD_TEST_ANSWER_FILE}")
      ENDIF()

      # set the precision for floating point failure comparison (number of digits, default 4)
      IF(SUNDIALS_DEVTESTS_FLOAT_PRECISION)
        LIST(APPEND TEST_ARGS "--floatprecision=${SUNDIALS_DEVTESTS_FLOAT_PRECISION}")
      ELSEIF(NOT ("${SUNDIALS_ADD_TEST_FLOAT_PRECISION}" STREQUAL ""))
        LIST(APPEND TEST_ARGS "--floatprecision=${SUNDIALS_ADD_TEST_FLOAT_PRECISION}")
      ENDIF()

      # set the integer percentage difference for failure comparison (default 10%)
      IF(SUNDIALS_DEVTESTS_INTEGER_PRECISION)
        LIST(APPEND TEST_ARGS "--integerpercentage=${SUNDIALS_DEVTESTS_INTEGER_PRECISION}")
      ELSEIF(NOT ("${SUNDIALS_ADD_TEST_INTEGER_PERCENTAGE}" STREQUAL ""))
        LIST(APPEND TEST_ARGS "--integerpercentage=${SUNDIALS_ADD_TEST_INTEGER_PERCENTAGE}")
      ENDIF()

      # create test case with the corresponding test runner command and arguments
      # all tests are added during development and only unlabeled tests when released
      ADD_TEST(NAME ${NAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTRUNNER} ${TEST_ARGS})

    elseif("${SUNDIALS_ADD_TEST_EXAMPLE_TYPE}" STREQUAL "")

      # convert string to list
      IF(NOT ("${SUNDIALS_ADD_TEST_TEST_ARGS}" STREQUAL ""))
        STRING(REPLACE " " ";" TEST_ARGS "${SUNDIALS_ADD_TEST_TEST_ARGS}")
      ENDIF()

      # check if this test is run with MPI and add the test run command
      if(NOT ("${SUNDIALS_ADD_TEST_MPI_NPROCS}" STREQUAL "") AND
          NOT ("${MPIEXEC_EXECUTABLE}" STREQUAL ""))
        ADD_TEST(NAME ${NAME} COMMAND ${MPIEXEC_EXECUTABLE} -n ${SUNDIALS_ADD_TEST_MPI_NPROCS} $<TARGET_FILE:${EXECUTABLE}> ${TEST_ARGS})
      else()
        ADD_TEST(NAME ${NAME} COMMAND $<TARGET_FILE:${EXECUTABLE}> ${TEST_ARGS})
      endif()

    endif()

  endif()

ENDMACRO()
