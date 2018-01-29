# ---------------------------------------------------------------
# Author:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2013, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
#
# SUNDIALS_ADD_TEST(<test name> <executable>)
# 
# CMake macro to add a Sundials regression test. Keyword input
# arguments can be added after <executable> to set regression
# test options (see oneValueArgs and multiValueArgs below).
#
# The executable is run and its output compared with the output in the
# test/answers directory. If the output differs significantly then the
# test fails. The default signicance is 4 decimal points for floating
# values and 10% for integer values.
# ---------------------------------------------------------------

IF(EXAMPLES_ENABLED)

  find_package(PythonInterp)
  IF(${PYTHON_VERSION_MAJOR} LESS 3)
    IF(${PYTHON_VERSION_MINOR} LESS 7)
      PRINT_WARNING("Python version must be 2.7.x or greater to run regression tests"
                    "Examples will build but 'make test' will fail.")
    ENDIF()
  ENDIF()

  # look for the testRunner script in the test directory
  FIND_PROGRAM(TESTRUNNER testRunner PATHS test)

ENDIF(EXAMPLES_ENABLED)


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

  # command line arguments for the test runner script
  SET(TEST_ARGS
    "--verbose"
    "--testname=${NAME}" 
    "--executablename=$<TARGET_FILE:${EXECUTABLE}>"
    "--outputdir=${CMAKE_BINARY_DIR}/Testing/output"
    )

  # only check the return value, do not diff the output and answer files
  IF(NOT ${SUNDIALS_DEVTESTS} OR ${SUNDIALS_ADD_TEST_NODIFF})
    LIST(APPEND TEST_ARGS "--nodiff")
  ENDIF()

  # set the number of mpi tasks to use in parallel tests
  IF("${SUNDIALS_ADD_TEST_MPI_NPROCS}" STREQUAL "")
  ELSE()

    IF(MPI_ENABLE)
      IF(MPI_RUN_COMMAND MATCHES "srun")
	SET(RUN_COMMAND "srun -N1 -n${SUNDIALS_ADD_TEST_MPI_NPROCS} -ppdebug")
      ELSE(MPI_RUN_COMMAND MATCHES "srun")
	SET(RUN_COMMAND "${MPI_RUN_COMMAND} -n ${SUNDIALS_ADD_TEST_MPI_NPROCS}")
      ENDIF(MPI_RUN_COMMAND MATCHES "srun")
      
      LIST(APPEND TEST_ARGS "--runcommand=\"${RUN_COMMAND}\"")

    ENDIF(MPI_ENABLE)

  ENDIF()
  
  # set the test input args
  IF("${SUNDIALS_ADD_TEST_TEST_ARGS}" STREQUAL "")
  ELSE()
    STRING (REPLACE ";" " " USER_ARGS "${SUNDIALS_ADD_TEST_TEST_ARGS}")
    LIST(APPEND TEST_ARGS "--runargs=\"${USER_ARGS}\"")
  ENDIF()

  # set the test answer directory name (default is test/answers)
  IF("${SUNDIALS_ADD_TEST_ANSWER_DIR}" STREQUAL "")
  ELSE()
    LIST(APPEND TEST_ARGS "--answerdir=${SUNDIALS_ADD_TEST_ANSWER_DIR}")
  ENDIF()

  # set the test answer file name (default is test_name_test_agrs)
  IF("${SUNDIALS_ADD_TEST_ANSWER_FILE}" STREQUAL "")
  ELSE()
    LIST(APPEND TEST_ARGS "--answerfile=${SUNDIALS_ADD_TEST_ANSWER_FILE}")
  ENDIF()

  # set the precision for floating point failure comparison (number of digits, default 4)
  IF("${SUNDIALS_ADD_TEST_FLOAT_PRECISION}" STREQUAL "")
  ELSE()
    LIST(APPEND TEST_ARGS "--floatprecision=${SUNDIALS_ADD_TEST_FLOAT_PRECISION}")
  ENDIF()

  # set the integer percentage difference for failure comparison (default 10%)
  IF("${SUNDIALS_ADD_TEST_INTEGER_PERCENTAGE}" STREQUAL "")
  ELSE()
    LIST(APPEND TEST_ARGS "--integerpercentage=${SUNDIALS_ADD_TEST_INTEGER_PERCENTAGE}")
  ENDIF()

  # create test case with the corresponding test runner command and arguments
  # all tests are added during development and only unlabeled tests when released
  IF(${SUNDIALS_DEVTESTS} OR "${SUNDIALS_ADD_TEST_EXAMPLE_TYPE}" STREQUAL "")
    ADD_TEST(NAME ${NAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTRUNNER} ${TEST_ARGS})
  ENDIF()

ENDMACRO()


MACRO(SUNDIALS_ADD_TEST_INSTALL SOLVER EXECUTABLE)

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
  FILE(MAKE_DIRECTORY ${TEST_INSTALL_DIR}/${SOLVER})

  # command line arguments for the test runner script
  set(TEST_ARGS
    "--testname=${EXECUTABLE}" 
    "--executablename=./${EXECUTABLE}"
    "--outputdir=${TEST_INSTALL_DIR}/${SOLVER}"
    "--builddir=${SUNDIALS_ADD_TEST_INSTALL_EXAMPLE_DIR}"
    "--buildcmd=${CMAKE_COMMAND}"
    "--nodiff"
    )
  
  # add test_install target for this solver
  ADD_CUSTOM_TARGET(${SOLVER}_test_install
    COMMAND ${PYTHON_EXECUTABLE} ${TESTRUNNER} ${TEST_ARGS}
    COMMENT "Running ${SOLVER} installation tests"
    WORKING_DIRECTORY ${TEST_INSTALL_DIR}/${SOLVER})

  # make test_install depend on solver_test_install
  ADD_DEPENDENCIES(test_install ${SOLVER}_test_install)

ENDMACRO()