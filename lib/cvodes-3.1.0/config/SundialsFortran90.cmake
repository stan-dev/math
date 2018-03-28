# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2014, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# ---------------------------------------------------------------
# Fortran90-related tests for SUNDIALS CMake-based configuration.

set(F90_FOUND FALSE)
include(CMakeDetermineFortranCompiler)

if(CMAKE_Fortran_COMPILER)
  message(STATUS "Searching for a Fortran compiler... ${CMAKE_Fortran_COMPILER}")

  # If Fortran compiler flags are set using environemnt variables and both FFLAGS
  # and FCFLAGS are used, then check if the variables are the same. If they are
  # not the same then a fatal error occurs.
  # 
  # NOTE: This check must occur before 'enable_language(Fortran)' as it will use
  # the value of FFLAGS to set CMAKE_Fortran_FLAGS
  SET(ENV_FFLAGS "$ENV{FFLAGS}")
  SET(ENV_FCFLAGS "$ENV{FCFLAGS}")
  IF ((NOT "${ENV_FFLAGS}" STREQUAL "") AND
      (NOT "${ENV_FCFLAGS}" STREQUAL "") AND
      ("${CMAKE_Fortran_FLAGS}" STREQUAL ""))
    # check if environment variables are equal
    IF (NOT "${ENV_FFLAGS}" STREQUAL "${ENV_FCFLAGS}")
      PRINT_ERROR("FFLAGS='${ENV_FFLAGS}' and FCFLAGS='${ENV_FCFLAGS}' are both set but are not equal.")
    ENDIF()
  ENDIF()

  # Enable the language for next steps
  enable_language(Fortran)

  # show some cache variables
  MARK_AS_ADVANCED(CLEAR
    CMAKE_Fortran_COMPILER
    CMAKE_Fortran_FLAGS)

  # hide all build type specific flags
  MARK_AS_ADVANCED(FORCE
    CMAKE_Fortran_FLAGS_DEBUG
    CMAKE_Fortran_FLAGS_MINSIZEREL
    CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_RELWITHDEBINFO)

  # only show flags for the current build type
  # these flags are appended to CMAKE_Fortran_FLAGS
  IF(CMAKE_BUILD_TYPE)       
    IF(CMAKE_BUILD_TYPE MATCHES "Debug")
      MESSAGE("Appending Fortran debug flags")
      MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_FLAGS_DEBUG)
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
      MESSAGE("Appending Fortran min size release flags")
      MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_FLAGS_MINSIZEREL)
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "Release")
      MESSAGE("Appending Fortran release flags")
      MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_FLAGS_RELEASE)
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
      MESSAGE("Appending Fortran release with debug info flags")
      MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
    ENDIF()
  ENDIF()

  # Create the Fortran90Test directory
  set(Fortran90Test_DIR ${PROJECT_BINARY_DIR}/Fortran90Test)
  file(MAKE_DIRECTORY ${Fortran90Test_DIR})

  # Create a CMakeLists.txt file which will generate the "f90lib" library
  # and an executable "f90test"
  file(WRITE ${Fortran90Test_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
    "PROJECT(f90test Fortran)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n"
    "SET(CMAKE_Fortran_FLAGS_RELEASE \"${CMAKE_Fortran_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_Fortran_FLAGS_DEBUG \"${CMAKE_Fortran_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO \"${CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_Fortran_FLAGS_MINSIZE \"${CMAKE_Fortran_FLAGS_MINSIZE}\")\n"
    "ADD_LIBRARY(f90lib f90lib.f90)\n"
    "ADD_EXECUTABLE(f90test f90test.f90)\n"
    "TARGET_LINK_LIBRARIES(f90test f90lib)\n")

  # Create the Fortran source f90lib.f90 which defines two subroutines, "mysub" and "my_sub"
  file(WRITE ${Fortran90Test_DIR}/f90lib.f90
    "subroutine mysub\n"
    "  return\n"
    "end\n"
    "subroutine my_sub\n"
    "  return\n"
    "end\n")

  # Create the Fortran source f90test.f90 which calls "mysub" and "my_sub"
  file(WRITE ${Fortran90Test_DIR}/f90test.f90
    "program f90test\n"
    "  call mysub()\n"
    "  call my_sub()\n"
    "end\n")

  # Use TRY_COMPILE to make the targets "f90lib" and "f90test"
  try_compile(F90TEST_OK ${Fortran90Test_DIR} ${Fortran90Test_DIR}
    f90test OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${Fortran90Test_DIR}/CMakeFiles)

  # Proceed based on test results
  if(F90TEST_OK)
    message(STATUS "Trying to compile and link a simple Fortran90 program... OK")
    set(F90_FOUND TRUE)
  else(F90TEST_OK)
    message(STATUS "Trying to compile and link a simple Fortran90 program... FAILED")
  endif(F90TEST_OK)

else(CMAKE_Fortran_COMPILER)
  message(STATUS "Searching for a Fortran compiler... FAILED")
endif(CMAKE_Fortran_COMPILER)

