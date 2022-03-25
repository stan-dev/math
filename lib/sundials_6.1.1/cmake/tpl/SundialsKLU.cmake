# -----------------------------------------------------------------------------
# Programmer(s): Steven Smith and Cody J. Balos @ LLNL
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
# Module to find and setup KLU correctly.
# Created from the SundialsTPL.cmake template.
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

if(NOT DEFINED SUNDIALS_KLU_INCLUDED)
  set(SUNDIALS_KLU_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# KLU does not support single or extended precision
if(SUNDIALS_PRECISION MATCHES "SINGLE" OR SUNDIALS_PRECISION MATCHES "EXTENDED")
  print_error("KLU is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(KLU REQUIRED)

message(STATUS "KLU_LIBRARIES:   ${KLU_LIBRARIES}")
message(STATUS "KLU_INCLUDE_DIR: ${KLU_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(KLU_FOUND AND (NOT KLU_WORKS))
  # Do any checks which don't require compilation first.

  # Create the KLU_TEST directory
  set(KLU_TEST_DIR ${PROJECT_BINARY_DIR}/KLU_TEST)
  file(MAKE_DIRECTORY ${KLU_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${KLU_TEST_DIR}/CMakeLists.txt
  "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
  "PROJECT(ltest C)\n"
  "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
  "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
  "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
  "SET(CMAKE_C_STANDARD \"${CMAKE_C_STANDARD}\")\n"
  "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
  "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
  "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
  "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
  "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
  "INCLUDE_DIRECTORIES(${KLU_INCLUDE_DIR})\n"
  "ADD_EXECUTABLE(ltest ltest.c)\n"
  "TARGET_LINK_LIBRARIES(ltest ${KLU_LIBRARIES})\n")

  # Create a C source file which calls a KLU function
  file(WRITE ${KLU_TEST_DIR}/ltest.c
  "\#include \"klu.h\"\n"
  "int main(){\n"
  "klu_common Common;\n"
  "klu_defaults (&Common);\n"
  "return(0);\n"
  "}\n")

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${KLU_TEST_DIR}/CMakeFiles)

  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${KLU_TEST_DIR} ${KLU_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if KLU works... OK")
    set(KLU_WORKS TRUE CACHE BOOL "KLU works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if KLU works... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to KLU is not functional.")
  endif()

elseif(KLU_FOUND AND KLU_WORKS)
  message(STATUS "Skipped KLU tests, assuming KLU works with SUNDIALS. Set KLU_WORKS=FALSE to (re)run compatibility test.")
endif()
