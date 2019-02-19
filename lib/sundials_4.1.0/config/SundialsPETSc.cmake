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
# PETSc tests for SUNDIALS CMake-based configuration.
# 

### This is only set if running GUI - simply return first time enabled
IF(PETSC_DISABLED)
  SET(PETSC_DISABLED FALSE CACHE INTERNAL "GUI - now enabled" FORCE)
  RETURN()
ENDIF()

SET(PETSC_FOUND FALSE)

# set PETSC_LIBRARIES
include(FindPETSc)

# If we have the PETSC libraries, test them
if(PETSC_LIBRARIES)
  message(STATUS "Looking for PETSc libraries...")
  # Create the PETSCTest directory
  set(PETSCTest_DIR ${PROJECT_BINARY_DIR}/PETSCTest)
  file(MAKE_DIRECTORY ${PETSCTest_DIR})
  # Create a CMakeLists.txt file 
  file(WRITE ${PETSCTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR})\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${PETSC_LIBRARIES})\n")    
  # Create a C source file which calls a PETSC function
  file(WRITE ${PETSCTest_DIR}/ltest.c
    "\#include \"petscvec.h\"\n"
    "int main(){\n"
    "Vec x;\n"
    "VecCreate(PETSC_COMM_WORLD, &x);\n" 
    "return(0);\n"
    "}\n")
  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${PETSCTest_DIR} ${PETSCTest_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)
      
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${PETSCTest_DIR}/CMakeFiles)
  # Process test result
  if(LTEST_OK)
    message(STATUS "Checking if PETSc works... OK")
    set(PETSC_FOUND TRUE)
  else(LTEST_OK)
    message(STATUS "Checking if PETSc works... FAILED")
  endif(LTEST_OK)
else(PETSC_LIBRARIES)
  PRINT_WARNING("PETSC LIBRARIES NOT Found. Please check library path" "${PETSC_LIBRARY_DIR} ")
  message(STATUS "Looking for PETSc libraries... FAILED")
endif(PETSC_LIBRARIES)
