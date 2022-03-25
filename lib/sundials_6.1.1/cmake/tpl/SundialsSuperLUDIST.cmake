# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup SuperLU_DIST correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# <TPL>_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_SUPERLUDIST_INCLUDED)
  set(SUNDIALS_SUPERLUDIST_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# SuperLU_DIST only supports double precision
if(SUNDIALS_PRECISION MATCHES "SINGLE" OR SUNDIALS_PRECISION MATCHES "EXTENDED")
  print_error("SuperLU_DIST is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# Using SUPERLUDIST requires building with MPI enabled
if(ENABLE_SUPERLUDIST AND NOT ENABLE_MPI)
  print_error("MPI is required for SuperLU DIST support. Set ENABLE_MPI to ON.")
endif()

# Using SUPERLUDIST with OpenMP requires building with OpenMP enabled
if(ENABLE_SUPERLUDIST AND SUPERLUDIST_OpenMP AND NOT ENABLE_OPENMP)
  print_error("OpenMP is required for SuperLU DIST support. Set ENABLE_OPENMP to ON.")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# We need MPI for SuperLU_DIST support
include(SundialsMPI)

# SuperLU_DIST OpenMP node parallelism is on, make sure OpenMP as found and is
# at least version 4.5.
if(SUPERLUDIST_OpenMP)
  include(SundialsOpenMP)

  if(NOT OPENMP_FOUND)
    print_error("SUPERLUDIST_OpenMP is set to ON but OpenMP was not found.")
  elseif(NOT OPENMP45_FOUND)
    string(CONCAT ERRSTR "SuperLUDIST requires OpenMP 4.5+ but it was not found. "
      "If you are sure OpenMP 4.5+ is available set the OPENMP_DEVICE_WORKS "
      "advanced option to ON.")
    print_error(${ERRSTR})
  endif()
endif()

# Try to find SuperLU_DIST
find_package(SUPERLUDIST 6.1.1 REQUIRED)

message(STATUS "SUPERLUDIST_LIBRARIES:   ${SUPERLUDIST_LIBRARIES}")
message(STATUS "SUPERLUDIST_INCLUDE_DIR: ${SUPERLUDIST_INCLUDE_DIR}")
message(STATUS "SUPERLUDIST_INDEX_SIZE:  ${SUPERLUDIST_INDEX_SIZE}")
message(STATUS "SUPERLUDIST_OpenMP:      ${SUPERLUDIST_OpenMP}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# If we have the SuperLU_DIST libraries, test them
if(SUPERLUDIST_FOUND AND (NOT SUPERLUDIST_WORKS))

  # Check index size
  if(NOT (SUNDIALS_INDEX_SIZE STREQUAL SUPERLUDIST_INDEX_SIZE))
    set(_err_msg_string "SuperLU_DIST not functional due to index size mismatch:\n")
    string(APPEND _err_msg_string "SUNDIALS_INDEX_SIZE=${SUNDIALS_INDEX_SIZE}, but SuperLU_DIST was built with ${SUPERLUDIST_INDEX_SIZE}-bit indices\n")
    string(APPEND _err_msg_string "SUPERLUDIST_INCLUDE_DIR: ${SUPERLUDIST_INCLUDE_DIR}\n")
    print_error("${_err_msg_string}")
  endif()

  # Create the SUPERLUDIST_TEST directory
  set(SUPERLUDIST_TEST_DIR ${PROJECT_BINARY_DIR}/SUPERLUDIST_Test)
  file(MAKE_DIRECTORY ${SUPERLUDIST_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${SUPERLUDIST_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
    "PROJECT(ltest CXX)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_CXX_COMPILER \"${MPI_CXX_COMPILER}\")\n"
    "SET(CMAKE_CXX_STANDARD \"${CMAKE_CXX_STANDARD}\")\n"
    "SET(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")\n"
    "SET(CMAKE_CXX_FLAGS_RELEASE \"${CMAKE_CXX_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_CXX_FLAGS_DEBUG \"${CMAKE_CXX_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_CXX_FLAGS_RELWITHDEBUGINFO \"${CMAKE_CXX_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_CXX_FLAGS_MINSIZE \"${CMAKE_CXX_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(ltest ltest.cpp)\n"
    "TARGET_INCLUDE_DIRECTORIES(ltest PRIVATE ${SUPERLUDIST_INCLUDE_DIR})\n"
    "TARGET_LINK_LIBRARIES(ltest ${SUPERLUDIST_LIBRARIES})\n")

  # Create a CXX source file which calls a SuperLUDIST function
  # and also prints the size of the indices used.
  file(WRITE ${SUPERLUDIST_TEST_DIR}/ltest.cpp
    "\#include <superlu_ddefs.h>\n"
    "int main(){\n"
    "SuperMatrix *A;\n"
    "NRformat_loc *Astore;\n"
    "A = NULL;\n"
    "Astore = NULL;\n"
    "if (A != NULL || Astore != NULL) return(1);\n"
    "else return(0);\n"
    "}\n")

  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${SUPERLUDIST_TEST_DIR} ${SUPERLUDIST_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${SUPERLUDIST_TEST_DIR}/CMakeFiles)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if SuperLU_DIST works with SUNDIALS... OK")
    set(SUPERLUDIST_WORKS TRUE CACHE BOOL "SuperLU_DIST works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if SuperLU_DIST works with SUNDIALS... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to SuperLU_DIST is not functional.")
  endif()

elseif(SUPERLUDIST_FOUND AND SUPERLUDIST_WORKS)
  message(STATUS "Skipped SuperLU_DIST tests, assuming SuperLU_DIST works with SUNDIALS. Set SUPERLUDIST_WORKS=FALSE to (re)run compatibility test.")
endif()
