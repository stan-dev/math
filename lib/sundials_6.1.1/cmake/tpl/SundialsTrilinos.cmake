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
# Module to find and setup Trilinos correctly.
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

if(NOT DEFINED SUNDIALS_TRILINOS_INCLUDED)
  set(SUNDIALS_TRILINOS_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# Find Trilinos
find_package(Trilinos REQUIRED)

# Check if Trilinos was built with MPI
if(";${Trilinos_TPL_LIST};" MATCHES ";MPI;")
  set(Trilinos_MPI TRUE)
else()
  set(Trilinos_MPI FALSE)
endif()

# For XSDK compatibility, only use the user/spack provided compiler and flags to build
# SUNDIALS modules that use Trilinos. If we are not in XSDK mode, we can use the imported
# Trilinos compiler and flags by default, but allow the user to change it through CMake
# the Trilinos_INTERFACE_* options.

if(USE_XSDK_DEFAULTS)
  if(Trilinos_MPI AND MPI_CXX_FOUND)
    force_variable(Trilinos_INTERFACE_CXX_COMPILER     STRING "C++ compiler for Trilinos interface" "${MPI_CXX_COMPILER}")
    set(Trilinos_INTERFACE_MPI_CXX_FOUND ${Trilinos_MPI} CACHE INTERNAL "Is Trilinos interface C++ compiler MPI")
  else()
    force_variable(Trilinos_INTERFACE_CXX_COMPILER     STRING "C compiler for Trilinos interface" "${CMAKE_CXX_COMPILER}")
    set(Trilinos_INTERFACE_MPI_CXX_FOUND FALSE CACHE INTERNAL "Is Trilinos interface C++ compiler MPI")
  endif()
  if(Trilinos_MPI AND MPI_C_FOUND)
    force_variable(Trilinos_INTERFACE_C_COMPILER       STRING "C compiler for Trilinos interface" "${MPI_C_COMPILER}")
    set(Trilinos_INTERFACE_MPI_C_FOUND ${Trilinos_MPI} CACHE INTERNAL "Is Trilinos interface C compiler MPI")
  else()
    force_variable(Trilinos_INTERFACE_C_COMPILER       STRING "C compiler for Trilinos interface" "${CMAKE_C_COMPILER}")
    set(Trilinos_INTERFACE_MPI_C_FOUND FALSE CACHE INTERNAL "Is Trilinos interface C compiler MPI")
  endif()
  force_variable(Trilinos_INTERFACE_CXX_COMPILER_FLAGS STRING "C++ compiler flags specific to Trilinos interface" "")
  force_variable(Trilinos_INTERFACE_C_COMPILER_FLAGS   STRING "C compiler flags specific to Trilinos interface" "")
  force_variable(Trilinos_INTERFACE_MPIEXEC            STRING "MPI executable for Trilinos interface" "${MPIEXEC_EXECUTABLE}")
else()
  force_variable(Trilinos_INTERFACE_CXX_COMPILER       STRING "C++ compiler for Trilinos interface" "${Trilinos_CXX_COMPILER}")
  force_variable(Trilinos_INTERFACE_C_COMPILER         STRING "C compiler for Trilinos interface" "${Trilinos_C_COMPILER}")
  force_variable(Trilinos_INTERFACE_CXX_COMPILER_FLAGS STRING "C++ compiler flags for Trilinos interface" "${Trilinos_CXX_COMPILER_FLAGS}")
  force_variable(Trilinos_INTERFACE_C_COMPILER_FLAGS   STRING "C compiler flags for Trilinos interface" "${Trilinos_C_COMPILER_FLAGS}")
  force_variable(Trilinos_INTERFACE_MPIEXEC            STRING "MPI executable for Trilinos interface" "${Trilinos_MPI_EXEC}")
  set(Trilinos_INTERFACE_MPI_CXX_FOUND ${Trilinos_MPI} CACHE INTERNAL "Is Trilinos interface C++ compiler MPI")
  set(Trilinos_INTERFACE_MPI_C_FOUND ${Trilinos_MPI} CACHE INTERNAL "Is Trilinos interface C compiler MPI")
endif()

message(STATUS "Trilinos_MPI:          ${Trilinos_MPI}")
message(STATUS "Trilinos_LIBRARIES:    ${Trilinos_LIBRARIES}")
message(STATUS "Trilinos_INCLUDE_DIRS: ${Trilinos_INCLUDE_DIRS}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(Trilinos_FOUND AND (NOT Trilinos_WORKS))
  # Do any checks which don't require compilation first.

  # Create the Trilinos_TEST directory
  set(Trilinos_TEST_DIR ${PROJECT_BINARY_DIR}/Trilinos_TEST)
  file(MAKE_DIRECTORY ${Trilinos_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${Trilinos_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
    "PROJECT(ltest CXX)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_CXX_COMPILER \"${Trilinos_INTERFACE_CXX_COMPILER}\")\n"
    "SET(CMAKE_CXX_STANDARD \"${CMAKE_CXX_STANDARD}\")\n"
    "SET(CMAKE_CXX_FLAGS \"${Trilinos_INTERFACE_CXX_COMPILER_FLAGS}\")\n"
    "SET(Trilinos_DIR \"${Trilinos_DIR}\")\n"
    "INCLUDE(FindPackageHandleStandardArgs)\n"
    "INCLUDE(${PROJECT_SOURCE_DIR}/cmake/tpl/FindTrilinos.cmake)\n"
    "ADD_EXECUTABLE(ltest ltest.cpp)\n"
    "TARGET_LINK_LIBRARIES(ltest SUNDIALS::TRILINOS)\n")

  # Create a C++ source file which calls a Trilinos function
  file(WRITE ${Trilinos_TEST_DIR}/ltest.cpp
  "#include <Tpetra_Version.hpp>\n"
  "int main(){\n"
  "std::cout << Tpetra::version() << std::endl;\n"
  "return(0);\n"
  "}\n")

  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${Trilinos_TEST_DIR} ${Trilinos_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if Trilinos works with SUNDIALS... OK")
    set(Trilinos_WORKS TRUE CACHE BOOL "Trilinos works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if Trilinos works with SUNDIALS... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to Trilinos is not functional.")
  endif()

elseif(Trilinos_FOUND AND Trilinos_WORKS)
  message(STATUS "Skipped Trilinos tests, assuming Trilinos works with SUNDIALS. Set Trilinos_WORKS=FALSE to (re)run compatibility test.")
endif()
