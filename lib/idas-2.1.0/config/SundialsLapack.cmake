# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
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
# BLAS/LAPACK tests for SUNDIALS CMake-based configuration.

SET(LAPACK_FOUND FALSE)

# If LAPACK libraries are undefined, try to find them (if we have
# a working Fortran compiler) or look for them in the most
# obvious place...
if(NOT LAPACK_LIBRARIES)
  if(F77_FOUND)
    include(FindLAPACK)
  else(F77_FOUND)
    find_library(LAPACK_LIBRARIES
      NAMES lapack
      PATHS /usr/lib /usr/local/lib
      "$ENV{ProgramFiles}/LAPACK/Lib"
      )
  endif(F77_FOUND)

  # If the xSDK flag is used, set it to what was found
  if(LAPACK_LIBRARIES AND TPL_ENABLE_LAPACK)
    SET(DOCSTR "Lapack library")
    FORCE_VARIABLE(TPL_LAPACK_LIBRARIES STRING "${DOCSTR}" "${LAPACK_LIBRARIES}")
  endif()
endif()

# If we have the LAPACK libraries, test them
if(LAPACK_LIBRARIES)
  message(STATUS "Looking for LAPACK libraries... OK")

  # Create the LapackTest directory
  set(LapackTest_DIR ${PROJECT_BINARY_DIR}/LapackTest)
  file(MAKE_DIRECTORY ${LapackTest_DIR})

  # Create a CMakeLists.txt file 
  file(WRITE ${LapackTest_DIR}/CMakeLists.txt
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

  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${LapackTest_DIR} ${LapackTest_DIR}
    ltest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${LapackTest_DIR}/CMakeFiles)

  # Process test result
  if(LTEST_OK)
    message(STATUS "Checking if LAPACK works... OK")
    set(LAPACK_FOUND TRUE)

    # get path to LAPACK library to use in generated makefiles for examples, if
    # LAPACK_LIBRARIES contains multiple items only use the path of the first entry
    list(LENGTH LAPACK_LIBRARIES len)
    if(len EQUAL 1)
      get_filename_component(LAPACK_LIBRARY_DIR ${LAPACK_LIBRARIES} PATH)
    else()
      list(GET LAPACK_LIBRARIES 0 TMP_LAPACK_LIBRARIES)
      get_filename_component(LAPACK_LIBRARY_DIR ${TMP_LAPACK_LIBRARIES} PATH)
    endif()

  else(LTEST_OK)
    message(STATUS "Checking if LAPACK works... FAILED")
  endif(LTEST_OK)

else(LAPACK_LIBRARIES)
  message(STATUS "Looking for LAPACK libraries... FAILED")
endif(LAPACK_LIBRARIES)
