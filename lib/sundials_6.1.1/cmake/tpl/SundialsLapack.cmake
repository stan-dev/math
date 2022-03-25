# -----------------------------------------------------------------------------
# Programmer(s): Radu Serban and Cody J. Balos @ LLNL
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
# Module to find and setup LAPACK/BLAS corrrectly.
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

if(NOT DEFINED SUNDIALS_LAPACK_INCLUDED)
  set(SUNDIALS_LAPACK_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# LAPACK does not support extended precision
if(ENABLE_LAPACK AND SUNDIALS_PRECISION MATCHES "EXTENDED")
  print_error("LAPACK is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# If LAPACK libraries are undefined, try to find them.
if(NOT LAPACK_LIBRARIES)
  find_package(LAPACK REQUIRED)
endif()

# If we have the LAPACK libraries, display progress message.
if(LAPACK_LIBRARIES)
  message(STATUS "Looking for LAPACK libraries... OK")
  set(LAPACK_FOUND TRUE)
else()
  message(STATUS "Looking for LAPACK libraries... FAILED")
  set(LAPACK_FOUND FALSE)
endif()

message(STATUS "LAPACK_LIBRARIES:  ${LAPACK_LIBRARIES}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# If we have the LAPACK libraries, determine if they work.
if(LAPACK_LIBRARIES AND (NOT LAPACK_WORKS))
  # Create the LapackTest directory
  set(LapackTest_DIR ${PROJECT_BINARY_DIR}/LapackTest)
  file(MAKE_DIRECTORY ${LapackTest_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${LapackTest_DIR}/CMakeLists.txt
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
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${LAPACK_LIBRARIES})\n")

  # Create a C source file which calls a Blas function (dcopy) and an Lapack function (dgetrf)
  file(WRITE ${LapackTest_DIR}/ltest.c
    "${F77_MANGLE_MACRO1}\n"
    "#define dcopy_f77 SUNDIALS_F77_FUNC(dcopy, DCOPY)\n"
    "#define dgetrf_f77 SUNDIALS_F77_FUNC(dgetrf, DGETRF)\n"
    "extern void dcopy_f77(int *n, const double *x, const int *inc_x, double *y, const int *inc_y);\n"
    "extern void dgetrf_f77(const int *m, const int *n, double *a, int *lda, int *ipiv, int *info);\n"
    "int main(){\n"
    "int n=1;\n"
    "double x, y;\n"
    "dcopy_f77(&n, &x, &n, &y, &n);\n"
    "dgetrf_f77(&n, &n, &x, &n, &n, &n);\n"
    "return(0);\n"
    "}\n")

  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${LapackTest_DIR} ${LapackTest_DIR}
    ltest OUTPUT_VARIABLE COMPILE_OUTPUT)

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${LapackTest_DIR}/CMakeFiles)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if LAPACK works with SUNDIALS... OK")
    set(LAPACK_WORKS TRUE CACHE BOOL "LAPACK works with SUNDIALS as configured" FORCE)

    # get path to LAPACK library to use in generated makefiles for examples, if
    # LAPACK_LIBRARIES contains multiple items only use the path of the first entry
    list(LENGTH LAPACK_LIBRARIES len)
    if(len EQUAL 1)
      get_filename_component(LAPACK_LIBRARY_DIR ${LAPACK_LIBRARIES} PATH)
    else()
      list(GET LAPACK_LIBRARIES 0 TMP_LAPACK_LIBRARIES)
      get_filename_component(LAPACK_LIBRARY_DIR ${TMP_LAPACK_LIBRARIES} PATH)
    endif()
  else(COMPILE_OK)
    set(LAPACK_WORKS FALSE CACHE BOOL "LAPACK does not work with SUNDIALS as configured" FORCE)
    message(STATUS "Checking if LAPACK works with SUNDIALS... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to LAPACK is not functional.")
  endif()

elseif(LAPACK_LIBRARIES AND LAPACK_WORKS)
  message(STATUS "Skipped LAPACK tests, assuming LAPACK works with SUNDIALS. Set LAPACK_WORKS=FALSE to (re)run compatibility test.")
endif()
