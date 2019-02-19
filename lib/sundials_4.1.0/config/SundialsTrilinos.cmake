# -----------------------------------------------------------------------------
# Programmer: Cody J. Balos and Slaven Peles @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Trilinos configuration and tests for SUNDIALS CMake-based configuration.
# 
# This also creates the variables:
#     Trilinos_INTERFACE_C_COMPILER
#     Trilinos_INTERFACE_CXX_COMPILER
#     Trilinos_INTERFACE_MPIEXEC
# The variables should be used to set the compiler when building targets
# that use Trilinos, and to set the MPIEXEC_EXECUTABLE for MPI+Trilinos
# examples.
#
# Additionally the variables:
#     Trilinos_INTERFACE_C_COMPILER_FLAGS
#     Trilinos_INTERFACE_CXX_COMPILER_FLAGS
# are created. The variables should be used to set the compiler flags
# when building targets that use Trilinos.
#
# The variable
#     Trilinos_FUNCTIONAL
# is created and indicates if the Trilinos installation found was usable
# under the current configuration of SUNDIALS.
#-----------------------------------------------------------------------------

set(Trilinos_FUNCTIONAL FALSE)

# Find Trilinos
include(FindTrilinos)

# Test Trilinos
if(Trilinos_FOUND AND TARGET Trilinos::Trilinos)

  # For XSDK compatibility, only use the user/spack provided compiler and flags to build
  # SUNDIALS modules that use Trilinos. If we are not in XSDK mode, we can use the imported
  # Trilinos compiler and flags by default, but allow the user to change it through CMake
  # the Trilinos_INTERFACE_* options.

  if(USE_XSDK_DEFAULTS)
    if(Trilinos_MPI AND MPI_CXX_FOUND)
      SHOW_VARIABLE(Trilinos_INTERFACE_CXX_COMPILER     STRING "C++ compiler for Trilinos interface" "${MPI_CXX_COMPILER}")
      set(Trilinos_INTERFACE_MPI_CXX_FOUND ${Trilinos_MPI} CACHE INTERNAL "Is Trilinos interface C++ compiler MPI")
    else()
      SHOW_VARIABLE(Trilinos_INTERFACE_CXX_COMPILER     STRING "C compiler for Trilinos interface" "${CMAKE_CXX_COMPILER}")
      set(Trilinos_INTERFACE_MPI_CXX_FOUND FALSE CACHE INTERNAL "Is Trilinos interface C++ compiler MPI")
    endif()
    if(Trilinos_MPI AND MPI_C_FOUND)
      SHOW_VARIABLE(Trilinos_INTERFACE_C_COMPILER       STRING "C compiler for Trilinos interface" "${MPI_C_COMPILER}")
      set(Trilinos_INTERFACE_MPI_C_FOUND ${Trilinos_MPI} CACHE INTERNAL "Is Trilinos interface C compiler MPI")
    else()
      SHOW_VARIABLE(Trilinos_INTERFACE_C_COMPILER       STRING "C compiler for Trilinos interface" "${CMAKE_C_COMPILER}")
      set(Trilinos_INTERFACE_MPI_C_FOUND FALSE CACHE INTERNAL "Is Trilinos interface C compiler MPI")
    endif()
    SHOW_VARIABLE(Trilinos_INTERFACE_CXX_COMPILER_FLAGS STRING "C++ compiler flags specific to Trilinos interface" "")
    SHOW_VARIABLE(Trilinos_INTERFACE_C_COMPILER_FLAGS   STRING "C compiler flags specific to Trilinos interface" "")
    SHOW_VARIABLE(Trilinos_INTERFACE_MPIEXEC            STRING "MPI executable for Trilinos interface" "${MPIEXEC_EXECUTABLE}")
  else()
    SHOW_VARIABLE(Trilinos_INTERFACE_CXX_COMPILER       STRING "C++ compiler for Trilinos interface" "${Trilinos_CXX_COMPILER}")
    SHOW_VARIABLE(Trilinos_INTERFACE_C_COMPILER         STRING "C compiler for Trilinos interface" "${Trilinos_C_COMPILER}")
    SHOW_VARIABLE(Trilinos_INTERFACE_CXX_COMPILER_FLAGS STRING "C++ compiler flags for Trilinos interface" "${Trilinos_CXX_COMPILER_FLAGS}")
    SHOW_VARIABLE(Trilinos_INTERFACE_C_COMPILER_FLAGS   STRING "C compiler flags for Trilinos interface" "${Trilinos_C_COMPILER_FLAGS}")
    SHOW_VARIABLE(Trilinos_INTERFACE_MPIEXEC            STRING "MPI executable for Trilinos interface" "${Trilinos_MPI_EXEC}")
    set(Trilinos_INTERFACE_MPI_CXX_FOUND ${Trilinos_MPI} CACHE INTERNAL "Is Trilinos interface C++ compiler MPI")
    set(Trilinos_INTERFACE_MPI_C_FOUND ${Trilinos_MPI} CACHE INTERNAL "Is Trilinos interface C compiler MPI")
  endif()
  mark_as_advanced(FORCE Trilinos_INTERFACE_CXX_COMPILER
                         Trilinos_INTERFACE_C_COMPILER
                         Trilinos_INTERFACE_CXX_COMPILER_FLAGS 
                         Trilinos_INTERFACE_C_COMPILER_FLAGS
                         Trilinos_INTERFACE_MPIEXEC)
  
  # Begin testing the Trilinos libraries with the compiler settings
  # that will be used when building Trilinos modules.
  message(STATUS "Testing Trilinos libraries...")

  # Create the TrilinosTest directory
  set(TrilinosTest_DIR ${PROJECT_BINARY_DIR}/TrilinosTest)
  file(MAKE_DIRECTORY ${TrilinosTest_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${TrilinosTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 3.1.3)\n"
    "PROJECT(ltest CXX)\n"
    "SET(Trilinos_DIR ${Trilinos_DIR})\n"
    "INCLUDE(${CMAKE_SOURCE_DIR}/config/FindTrilinos.cmake)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_C_COMPILER ${Trilinos_INTERFACE_C_COMPILER})\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS} ${Trilinos_INTERFACE_C_COMPILER_FLAGS}\")\n"
    "SET(CMAKE_CXX_COMPILER ${Trilinos_INTERFACE_CXX_COMPILER})\n"
    "SET(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS} ${Trilinos_INTERFACE_CXX_COMPILER_FLAGS}\")\n"
    "ADD_EXECUTABLE(ltest ltest)\n"
    "TARGET_LINK_LIBRARIES(ltest Trilinos::Trilinos)\n")

  # Create a C++ source file which calls a Trilinos function
  file(WRITE ${TrilinosTest_DIR}/ltest.cpp
    "#include <Tpetra_Version.hpp>\n"
    "int main(){\n"
    "std::cout << Tpetra::version() << std::endl;\n"
    "return(0);\n"
    "}\n")

  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${TrilinosTest_DIR} ${TrilinosTest_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory and the ltest binary.
  file(REMOVE_RECURSE ${TrilinosTest_DIR}/CMakeFiles)
  file(REMOVE_RECURSE ${TrilinosTest_DIR}/ltest)

  # Process test result
  if(LTEST_OK)
    message(STATUS "Testing Trilinos libraries... OK")
    set(Trilinos_FUNCTIONAL TRUE)
  else(LTEST_OK)
    message(STATUS "Testing Trilinos libraries... FAILED")
  endif(LTEST_OK)
  
endif()

