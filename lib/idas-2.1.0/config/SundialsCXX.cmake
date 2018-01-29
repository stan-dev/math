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
# C++-related tests for SUNDIALS CMake-based configuration.

set(CXX_FOUND FALSE)

include(CMakeDetermineCXXCompiler)

if(CMAKE_CXX_COMPILER)
  message(STATUS "Searching for a C++ compiler... ${CMAKE_CXX_COMPILER}")

  # Enable the language for next steps
  enable_language(CXX)

  # show some cache variables
  MARK_AS_ADVANCED(CLEAR
    CMAKE_CXX_COMPILER
    CMAKE_CXX_FLAGS)

  # hide all build type specific flags
  MARK_AS_ADVANCED(FORCE
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  # only show flags for the current build type
  # these flags are appended to CMAKE_CXX_FLAGS
  IF(CMAKE_BUILD_TYPE)
    IF(CMAKE_BUILD_TYPE MATCHES "Debug")
      MESSAGE("Appending CXX debug flags")
      MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_DEBUG)
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
      MESSAGE("Appending CXX min size release flags")
      MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_MINSIZEREL)
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "Release")
      MESSAGE("Appending CXX release flags")
      MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_RELEASE)
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
      MESSAGE("Appending CXX release with debug info flags")
      MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    ENDIF()
  ENDIF()

  # Create the CXXTest directory
  set(CXXTest_DIR ${PROJECT_BINARY_DIR}/CXXTest)
  file(MAKE_DIRECTORY ${CXXTest_DIR})

  # Create a CMakeLists.txt file which will generate the executable "cxxtest"
  file(WRITE ${CXXTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
    "PROJECT(cxxtest CXX)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")\n"
    "SET(CMAKE_CXX_FLAGS_RELEASE \"${CMAKE_CXX_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_CXX_FLAGS_DEBUG \"${CMAKE_CXX_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_CXX_FLAGS_RELWITHDEBUGINFO \"${CMAKE_CXX_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_CXX_FLAGS_MINSIZE \"${CMAKE_CXX_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(cxxtest cxxtest.cpp)\n")

  # Create the C++ source cxxtest.cpp which does some simple calls
  file(WRITE ${CXXTest_DIR}/cxxtest.cpp
    "#include <string>\n"
    "int main(){\n"
    "std::string c;\n"
    "return(0);\n"
    "}\n")

  # Use TRY_COMPILE to make the target "cxxtest"
  try_compile(CXXTEST_OK ${CXXTest_DIR} ${CXXTest_DIR}
    ftest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${CXXTest_DIR}/CMakeFiles)

  # Proceed based on test results
  if(CXXTEST_OK)
    message(STATUS "Trying to compile and link a simple C++ program... OK")
    set(CXX_FOUND TRUE)
  endif(CXXTEST_OK)

else(CMAKE_CXX_COMPILER)
  message(STATUS "Searching for a C++ compiler... FAILED")
endif(CMAKE_CXX_COMPILER)

