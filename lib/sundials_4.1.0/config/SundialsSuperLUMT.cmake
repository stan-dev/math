# ---------------------------------------------------------------
# Programmer:  Eddy Banks @ LLNL
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
# SUPERLUMT tests for SUNDIALS CMake-based configuration.
#    - loosely based on SundialsLapack.cmake
#


### This is only set if running GUI - simply return first time enabled
IF(SUPERLUMT_DISABLED)
  SET(SUPERLUMT_DISABLED FALSE CACHE INTERNAL "GUI - SUPERLUMT now enabled" FORCE)
  RETURN()
ENDIF()

SET(SUPERLUMT_FOUND FALSE)

# set SUPERLUMT_LIBRARIES
include(FindSUPERLUMT)

# If we have the SUPERLUMT libraries, test them
if(SUPERLUMT_LIBRARY AND SUPERLUMT_LIBRARIES)
  message(STATUS "Looking for SUPERLUMT libraries... OK")

  # Create the SUPERLUMT_TEST directory
  set(SUPERLUMT_TEST_DIR ${PROJECT_BINARY_DIR}/SUPERLUMT_TEST)
  file(MAKE_DIRECTORY ${SUPERLUMT_TEST_DIR})

  # Create a CMakeLists.txt file 
  file(WRITE ${SUPERLUMT_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "INCLUDE_DIRECTORIES(${SUPERLUMT_INCLUDE_DIR})\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${SUPERLUMT_LIBRARIES})\n")    

  # Create a C source file which calls a SUPERLUMT function
  file(WRITE ${SUPERLUMT_TEST_DIR}/ltest.c
    "\#include \"slu_mt_ddefs.h\"\n"
#    "\#include \"pdsp_defs.h\"\n"
    "int main(){\n"
    "SuperMatrix A;\n"
    "NCformat *Astore;\n" 
    "return(0);\n"
    "}\n")

  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${SUPERLUMT_TEST_DIR} ${SUPERLUMT_TEST_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)
      
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${SUPERLUMT_TEST_DIR}/CMakeFiles)

  # Process test result
  if(LTEST_OK)
    message(STATUS "Checking if SUPERLUMT works... OK")
    set(SUPERLUMT_FOUND TRUE)
  else(LTEST_OK)
    message(STATUS "Checking if SUPERLUMT works... FAILED")
  endif(LTEST_OK)

else()
  PRINT_WARNING("SUPERLUMT LIBRARIES NOT Found. Please check library path" "${SUPERLUMT_LIBRARY_DIR}")
  message(STATUS "Looking for SUPERLUMT libraries... FAILED")
endif()
