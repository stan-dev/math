# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# Module to find and setup <TPL> correctly.
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

if(NOT DEFINED SUNDIALS_XBRAID_INCLUDED)
  set(SUNDIALS_XBRAID_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# Using XBRAID requires building with MPI enabled
if(NOT ENABLE_MPI)
  message(FATAL_ERROR
    "MPI is required for XBraid support. Set ENABLE_MPI to ON.")
endif()

# XBraid does not support single or extended precision
if(SUNDIALS_PRECISION MATCHES "SINGLE" OR SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(FATAL_ERROR
    "XBraid is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# XBraid does not support 64-bit index sizes
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  message(FATAL_ERROR
    "XBraid is not compatible with ${SUNDIALS_INDEX_SIZE}-bit indices")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(XBRAID REQUIRED)

message(STATUS "XBRAID_LIBRARIES: ${XBRAID_LIBS}")
message(STATUS "XBRAID_INCLUDES:  ${XBRAID_INCS}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# Add works variable

if(XBRAID_FOUND AND (NOT XBRAID_WORKS))

  # Create the XBRAID_TEST directory
  set(XBRAID_TEST_DIR ${PROJECT_BINARY_DIR}/XBRAID_TEST)
  file(MAKE_DIRECTORY ${XBRAID_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${XBRAID_TEST_DIR}/CMakeLists.txt
    "cmake_minimum_required(VERSION ${CMAKE_VERSION})\n"
    "project(ltest C)\n"
    "set(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "set(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "set(CMAKE_C_COMPILER ${MPI_C_COMPILER})\n"
    "set(CMAKE_C_STANDARD ${CMAKE_C_STANDARD})\n"
    "set(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "set(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "set(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "set(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "set(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "add_executable(ltest ltest.c)\n"
    "target_include_directories(ltest PRIVATE \"${XBRAID_INCS}\")\n"
    "target_link_libraries(ltest \"${XBRAID_LIBS}\")\n"
    "target_link_libraries(ltest m)\n")

  # Create a C source file
  file(WRITE ${XBRAID_TEST_DIR}/ltest.c
    "\#include <stdlib.h>\n"
    "\#include \"braid.h\"\n"
    "int main(){\n"
    "braid_Int rand;\n"
    "rand = braid_Rand();\n"
    "if (rand < 0) return 1;\n"
    "return 0;\n"
    "}\n")

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${XBRAID_TEST_DIR}/CMakeFiles)

  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${XBRAID_TEST_DIR} ${XBRAID_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if XBRAID works... OK")
    set(XBRAID_WORKS TRUE CACHE BOOL "XBRAID works as configured" FORCE)
  else()
    message(STATUS "Checking if XBRAID works... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    message(FATAL_ERROR "XBRAID compile test failed.")
  endif()

  message(STATUS "XBRAID tests passed")

else()
  message(STATUS "Skipped XBRAID tests, assuming XBRAID works with SUNDIALS.")
endif()
