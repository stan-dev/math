# ---------------------------------------------------------------
# Programmer:  Slaven Peles @ LLNL, Jean Sexton @ SMU
#              Eddy Banks @ LLNL
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
# Hypre tests for SUNDIALS CMake-based configuration.
# 

### This is only set if running GUI - simply return first time enabled
IF(HYPRE_DISABLED)
  SET(HYPRE_DISABLED FALSE CACHE INTERNAL "GUI - now enabled" FORCE)
  RETURN()
ENDIF()

set(HYPRE_FOUND FALSE)

include(FindHypre)

if(UNIX)
  set(LINK_MATH_LIB "-lm")
endif()

if(HYPRE_LIBRARIES)
  message(STATUS "Looking for HYPRE LIBRARIES...")
  # Create the HYPRETest directory
  set(HYPRETest_DIR ${PROJECT_BINARY_DIR}/HYPRETest)
  file(MAKE_DIRECTORY ${HYPRETest_DIR})
  # Create a CMakeLists.txt file 
  file(WRITE ${HYPRETest_DIR}/CMakeLists.txt
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
    "SET(CMAKE_EXE_LINKER_FLAGS \"${LINK_MATH_LIB}\")\n"
    "INCLUDE_DIRECTORIES(${HYPRE_INCLUDE_DIR})\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${HYPRE_LIBRARIES})\n")    
  # Create a C source file which calls a hypre function
  file(WRITE ${HYPRETest_DIR}/ltest.c
    "\#include \"HYPRE_parcsr_ls.h\"\n"
    "int main(){\n"
    "HYPRE_ParVector par_b;\n"
    "HYPRE_IJVector b;\n"
    "return(0);\n"
    "}\n")
  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${HYPRETest_DIR} ${HYPRETest_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)
      
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${HYPRETest_DIR}/CMakeFiles)
  # Process test result
  if(LTEST_OK)
    message(STATUS "Checking if HYPRE works... OK")
    set(HYPRE_FOUND TRUE)
  else(LTEST_OK)
    message(STATUS "Checking if HYPRE works... FAILED")
  endif(LTEST_OK)
else(HYPRE_LIBRARIES)
  PRINT_WARNING("HYPRE LIBRARIES NOT Found. Please check library path" "${HYPRE_LIBRARY_DIR} ")
  message(STATUS "Looking for HYPRE LIBRARY... FAILED")
endif(HYPRE_LIBRARIES)
