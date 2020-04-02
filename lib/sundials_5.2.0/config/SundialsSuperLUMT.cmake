# ------------------------------------------------------------------------------
# Programmer(s): Eddy Banks @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# Sundials module to find and configure SuperLU_MT correctly.
# ------------------------------------------------------------------------------

### This is only set if running GUI - simply return first time enabled
if(SUPERLUMT_DISABLED)
  set(SUPERLUMT_DISABLED FALSE CACHE INTERNAL "GUI - SUPERLUMT now enabled" FORCE)
  return()
endif()

# --- Check for dependencies --- #

# SuperLU_MT does not support extended precision
if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  print_error("SuperLU_MT is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# --- Find SuperLU_MT and test it --- #

# Try to find SuperLU_MT
find_package(SUPERLUMT REQUIRED)

# If we have the SuperLU_MT libraries, test them
if(SUPERLUMT_FOUND)

  # >>>>>>> Need to add check for SuperLU_MT integer type <<<<<<<

  # Create the SUPERLUMT_TEST directory
  set(SUPERLUMT_TEST_DIR ${PROJECT_BINARY_DIR}/SUPERLUMT_TEST)
  file(MAKE_DIRECTORY ${SUPERLUMT_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${SUPERLUMT_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 3.1.3)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_INCLUDE_DIRECTORIES(ltest PRIVATE ${SUPERLUMT_INCLUDE_DIR})\n"
    "TARGET_LINK_LIBRARIES(ltest ${SUPERLUMT_LIBRARIES})\n")

  # Create a C source file which calls a SuperLU_MT function
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
    message(STATUS "Checking if SuperLU_MT works... OK")
    set(SUPERLUMT_FOUND TRUE)
  else()
    message(STATUS "Checking if SuperLU_MT works... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SuperLU_MT not functional - support will not be provided.")
  endif()

  # sundials_config.h symbols
  set(SUNDIALS_SUPERLUMT TRUE)
  set(SUNDIALS_SUPERLUMT_THREAD_TYPE ${SUPERLUMT_THREAD_TYPE})

else()

  # sundials_config.h symbols
  set(SUNDIALS_SUPERLUMT FALSE)

endif()
