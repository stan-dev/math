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
# Module to find and setup SUPERLUMT correctly.
# Created from the SundialsTPL.cmake.template template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# TPL_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_SUPERLUMT_INCLUDED)
  set(SUNDIALS_SUPERLUMT_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# SUPERLUMT does not support extended precision
if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  print_error("SUPERLUMT is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(SUPERLUMT REQUIRED)

message(STATUS "SUPERLUMT_LIBRARIES:   ${SUPERLUMT_LIBRARIES}")
message(STATUS "SUPERLUMT_INCLUDE_DIR: ${SUPERLUMT_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUPERLUMT_FOUND AND (NOT SUPERLUMT_WORKS))

  # Create the SUPERLUMT_TEST directory
  set(SUPERLUMT_TEST_DIR ${PROJECT_BINARY_DIR}/SUPERLUMT_TEST)
  file(MAKE_DIRECTORY ${SUPERLUMT_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${SUPERLUMT_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
    "SET(CMAKE_C_STANDARD ${CMAKE_C_STANDARD})\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_INCLUDE_DIRECTORIES(ltest PRIVATE ${SUPERLUMT_INCLUDE_DIR})\n"
    "TARGET_LINK_LIBRARIES(ltest ${SUPERLUMT_LIBRARIES})\n")

  # Create a C source file which calls a SUPERLUMT function
  file(WRITE ${SUPERLUMT_TEST_DIR}/ltest.c
    "\#include \"slu_mt_ddefs.h\"\n"
    "int main(){\n"
    "SuperMatrix *A;\n"
    "NCformat *Astore;\n"
    "A = NULL;\n"
    "Astore = NULL;\n"
    "if (A != NULL || Astore != NULL) return(1);\n"
    "else return(0);\n"
    "}\n")


  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${SUPERLUMT_TEST_DIR} ${SUPERLUMT_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${SUPERLUMT_TEST_DIR}/CMakeFiles)

 # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if SuperLU_MT works with SUNDIALS... OK")
    set(SUPERLUMT_WORKS TRUE CACHE BOOL "SuperLU_MT works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if SuperLU_MT works with SUNDIALS... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to SuperLU_MT is not functional.")
  endif()

elseif(SUPERLUMT_FOUND AND SUPERLUMT_WORKS)
  message(STATUS "Skipped SuperLU_MT tests, assuming SuperLU_MT works with SUNDIALS. Set SUPERLUMT_WORKS=FALSE to (re)run compatibility test.")
endif()
