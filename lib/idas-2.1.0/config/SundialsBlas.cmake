# ---------------------------------------------------------------
# Programmer:  David J. Gardner @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# BLAS tests for SUNDIALS CMake-based configuration. Based on 
# SundialsLapack.cmake
#

SET(BLAS_FOUND FALSE)

# If BLAS libraries are undefined, try to find them (if we have
# a working Fortran compiler) or look for them in the most
# obvious place...
if(NOT BLAS_LIBRARIES)
  if(F77_FOUND)
    include(FindBLAS)
  else(F77_FOUND)
    find_library(BLAS_LIBRARIES
      NAMES blas
      PATHS /usr/lib /usr/local/lib
      "$ENV{ProgramFiles}/BLAS/Lib"
      )
  endif(F77_FOUND)

  # If the xSDK flag is used, set it to what was found
  if(BLAS_LIBRARIES AND TPL_ENABLE_BLAS)
    SET(DOCSTR "Blas library")
    FORCE_VARIABLE(TPL_BLAS_LIBRARIES STRING "${DOCSTR}" "${BLAS_LIBRARIES}")
  endif()
endif()

# If we have the BLAS libraries, test them
if(BLAS_LIBRARIES)
  message(STATUS "Looking for BLAS libraries... OK")

  # Create the BlasTest directory
  set(BlasTest_DIR ${PROJECT_BINARY_DIR}/BlasTest)
  file(MAKE_DIRECTORY ${BlasTest_DIR})

  # Create a CMakeLists.txt file 
  file(WRITE ${BlasTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${BLAS_LIBRARIES})\n")

  # Create a C source file which calls a Blas function (dcopy)
  file(WRITE ${BlasTest_DIR}/ltest.c
    "${F77_MANGLE_MACRO1}\n"
    "#define dcopy_f77 SUNDIALS_F77_FUNC(dcopy, DCOPY)\n"
    "extern void dcopy_f77(int *n, const double *x, const int *inc_x, double *y, const int *inc_y);\n"
    "int main(){\n"
    "int n=1;\n"
    "double x, y;\n"
    "dcopy_f77(&n, &x, &n, &y, &n);\n"
    "return(0);\n"
    "}\n")

  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${BlasTest_DIR} ${BlasTest_DIR}
    ltest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${BlasTest_DIR}/CMakeFiles)

  # Process test result
  if(LTEST_OK)
    message(STATUS "Checking if BLAS works... OK")
    set(BLAS_FOUND TRUE)

    # get path to BLAS library to use in generated makefiles for examples
    # check length to protect against BLAS_LIBRARIES having multiple entries
    list(LENGTH BLAS_LIBRARIES len)
    if(len EQUAL 1)
      get_filename_component(BLAS_LIBRARY_DIR ${BLAS_LIBRARIES} PATH)
    endif()

  else(LTEST_OK)
    message(STATUS "Checking if BLAS works... FAILED")
  endif(LTEST_OK)

else(BLAS_LIBRARIES)
  message(STATUS "Looking for BLAS libraries... FAILED")
endif(BLAS_LIBRARIES)
