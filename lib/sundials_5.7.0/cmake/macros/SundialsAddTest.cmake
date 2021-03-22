# ------------------------------------------------------------------------------
# Programmer(s): Steven Smith and David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
#
# SUNDIALS_ADD_TEST(<test name> <executable>)
#
# CMake macro to add a SUNDIALS regression test. Keyword input arguments can be
# added after <executable> to set regression test options (see oneValueArgs and
# multiValueArgs below).
#
# When SUNDIALS_TEST_DEVTESTS is OFF (default) the executable is run and success
# or failure is determined by the executable return value (zero or non-zero
# respectively).
#
# When SUNDIALS_TEST_DEVTESTS is ON the executable is run and its output is
# compared with the corresponding .out file. If the output differs significantly
# then the test fails. The default level of significance is 4 decimal places for
# floating point values and 10% for integer values.
#
# The level of precision can be adjusted for an individual test with the
# FLOAT_PRECISION AND INTEGER_PERCENTAGE keyword inputs to the macro or globally
# for all tests with the cache variables SUNDIALS_TEST_FLOAT_PRECISION and
# SUNDIALS_TEST_INTEGER_PRECISION.
#
#  -D SUNDIALS_TEST_FLOAT_PRECISION=<number of digits>
#  -D SUNDIALS_TEST_INTEGER_PRECISION=<% difference>
#
# By default testing output is written to builddir/Testing/output and the .out
# answer file directory is set using the ANSWER_DIR keyword input to
# sourcedir/examples/package/testdir. These can be changed by setting the cache
# variables SUNDIALS_TEST_OUTPUT_DIR and SUNDIALS_TEST_ANSWER_DIR.
#
#  -D SUNDIALS_TEST_OUTPUT_DIR=<path to output directory>
#  -D SUNDIALS_TEST_ANSWER_DIR=<path to answer directory>
# ------------------------------------------------------------------------------

macro(SUNDIALS_ADD_TEST NAME EXECUTABLE)

  # macro options
  # NODIFF = do not diff the test output against an answer file
  set(options "NODIFF")

  # macro keyword inputs followed by a single value
  # MPI_NPROCS         = number of mpi tasks to use in parallel tests
  # FLOAT_PRECISION    = precision for floating point failure comparision (num digits),
  #                      to use the default, either don't provide the keyword, or
  #                      provide the value "default"
  # INTEGER_PRECENTAGE = integer percentage difference for failure comparison
  # ANSWER_DIR         = path to the directory containing the test answer file
  # ANSWER_FILE        = name of test answer file
  # EXAMPLE_TYPE       = release or develop examples
  set(oneValueArgs "MPI_NPROCS" "FLOAT_PRECISION" "INTEGER_PERCENTAGE"
    "ANSWER_DIR" "ANSWER_FILE" "EXAMPLE_TYPE")

  # macro keyword inputs followed by multiple values
  # TEST_ARGS = command line arguments to pass to the test executable
  set(multiValueArgs "TEST_ARGS")

  # parse inputs and create variables SUNDIALS_ADD_TEST_<keyword>
  cmake_parse_arguments(SUNDIALS_ADD_TEST
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # SGS add check to make sure parallel is integer

  # check that the test is not excluded
  if(NOT ("${SUNDIALS_ADD_TEST_EXAMPLE_TYPE}" STREQUAL "exclude"))

    if(SUNDIALS_TEST_DEVTESTS)

      # run all tests (standard and develop) with the test runner

      # command line arguments for the test runner script
      set(TEST_ARGS
        "--verbose"
        "--testname=${NAME}"
        "--executablename=$<TARGET_FILE:${EXECUTABLE}>"
        )

      # check for a non-default output directory
      if(SUNDIALS_TEST_OUTPUT_DIR)
        list(APPEND TEST_ARGS "--outputdir=${SUNDIALS_TEST_OUTPUT_DIR}")
      else()
        list(APPEND TEST_ARGS "--outputdir=${TEST_OUTPUT_DIR}")
      endif()

      # set a non-default answer directory (default is test/answers)
      if(SUNDIALS_TEST_ANSWER_DIR)
        list(APPEND TEST_ARGS "--answerdir=${SUNDIALS_TEST_ANSWER_DIR}")
      elseif(SUNDIALS_ADD_TEST_ANSWER_DIR)
        list(APPEND TEST_ARGS "--answerdir=${SUNDIALS_ADD_TEST_ANSWER_DIR}")
      endif()

      # set the test answer file name (default is test_name_test_agrs)
      if(SUNDIALS_ADD_TEST_ANSWER_FILE)
        list(APPEND TEST_ARGS "--answerfile=${SUNDIALS_ADD_TEST_ANSWER_FILE}")
      endif()

      # check if a diff is needed and if non-default precisions were provided
      if(SUNDIALS_ADD_TEST_NODIFF)
        # do not diff the output and answer files
        list(APPEND TEST_ARGS "--nodiff")
      else()
        # set a non-default floating point precision (number of digits, default 4)
        if(SUNDIALS_TEST_FLOAT_PRECISION)
          list(APPEND TEST_ARGS "--floatprecision=${SUNDIALS_TEST_FLOAT_PRECISION}")
        elseif(SUNDIALS_ADD_TEST_FLOAT_PRECISION AND
               (NOT SUNDIALS_ADD_TEST_FLOAT_PRECISION MATCHES "DEFAULT|default"))
          list(APPEND TEST_ARGS "--floatprecision=${SUNDIALS_ADD_TEST_FLOAT_PRECISION}")
        endif()
        # set a non-default integer precision (percent difference, default 10%)
        if(SUNDIALS_TEST_INTEGER_PRECISION)
          list(APPEND TEST_ARGS "--integerpercentage=${SUNDIALS_TEST_INTEGER_PRECISION}")
        elseif(SUNDIALS_ADD_TEST_INTEGER_PERCENTAGE)
          list(APPEND TEST_ARGS "--integerpercentage=${SUNDIALS_ADD_TEST_INTEGER_PERCENTAGE}")
        endif()
      endif()

      # check if this test is run with MPI and set the MPI run command
      if((SUNDIALS_ADD_TEST_MPI_NPROCS) AND (MPIEXEC_EXECUTABLE))
        if(MPIEXEC_EXECUTABLE MATCHES "srun")
          set(RUN_COMMAND "srun -N1 -n${SUNDIALS_ADD_TEST_MPI_NPROCS} -ppdebug")
        else()
          set(RUN_COMMAND "${MPIEXEC_EXECUTABLE} -n ${SUNDIALS_ADD_TEST_MPI_NPROCS}")
        endif()
        list(APPEND TEST_ARGS "--runcommand=\"${RUN_COMMAND}\"")
      endif()

      # set the test input args
      if(SUNDIALS_ADD_TEST_TEST_ARGS)
        string(REPLACE ";" " " USER_ARGS "${SUNDIALS_ADD_TEST_TEST_ARGS}")
        list(APPEND TEST_ARGS "--runargs=\"${USER_ARGS}\"")
      endif()

      # create test case with the corresponding test runner command and arguments
      # all tests are added during development and only unlabeled tests when released
      add_test(NAME ${NAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTRUNNER} ${TEST_ARGS})

    elseif(NOT SUNDIALS_ADD_TEST_EXAMPLE_TYPE)

      # if a test type was not set then it is a standard test that returns pass/fail

      # convert string to list
      if(SUNDIALS_ADD_TEST_TEST_ARGS)
        string(REPLACE " " ";" TEST_ARGS "${SUNDIALS_ADD_TEST_TEST_ARGS}")
      endif()

      # check if this test is run with MPI and add the test run command
      if((SUNDIALS_ADD_TEST_MPI_NPROCS) AND (MPIEXEC_EXECUTABLE))
        add_test(NAME ${NAME} COMMAND ${MPIEXEC_EXECUTABLE} -n ${SUNDIALS_ADD_TEST_MPI_NPROCS} $<TARGET_FILE:${EXECUTABLE}> ${TEST_ARGS})
      else()
        add_test(NAME ${NAME} COMMAND $<TARGET_FILE:${EXECUTABLE}> ${TEST_ARGS})
      endif()

    endif()

  endif()

endmacro()
