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
# MPI-C++ tests for SUNDIALS CMake-based configuration.

set(MPICXX_FOUND FALSE)

# Local variable indicating whether to test MPI
set(MPICXX_PERFORM_TEST FALSE)
# By default, we try to use the MPI compiler script
# Search for the MPICXX compiler script
find_program(MPI_MPICXX NAMES mpicxx DOC "mpicxx program")
if(MPI_MPICXX)
  message(STATUS "Looking for MPI C++ compiler script... ${MPI_MPICXX}")
  # Test the MPI compiler script
  set(MPICXX_PERFORM_TEST TRUE)
else(MPI_MPICXX)
  message(STATUS "Looking for MPI C++ compiler script... FAILED")
  # If not already available, search for MPI headers and libraries.
  # Define the following values
  #  MPI_INCLUDE_PATH = cached location of mpi.h
  #  MPI_LIBRARIES    = cached list of libraries to link in (mpi mpich etc)
  if(NOT MPI_LIBRARIES)
    find_path(MPI_INCLUDE_PATH mpi.h
      PATHS /usr/local/include 
      /usr/include 
      /usr/include/mpi
      /usr/local/mpi/include
      "$ENV{ProgramFiles}/MPICH/SDK/Include"
      "$ENV{ProgramFiles}/MPICH2/include"
      "C:/Program Files/MPICH/SDK/Include"
      )
    find_library(MPI_LIBRARIES
      NAMES mpich2 mpi mpich 
      PATHS /usr/lib /usr/local/lib /usr/local/mpi/lib
      "$ENV{ProgramFiles}/MPICH/SDK/Lib"
      "$ENV{ProgramFiles}/MPICH2/Lib"
      "C:/Program Files/MPICH/SDK/Lib" 
      )
    find_library(MPI_EXTRA_LIBRARIES 
      NAMES mpi++
      PATHS /usr/lib /usr/local/lib /usr/local/mpi/lib 
      "$ENV{ProgramFiles}/MPICH/SDK/Lib"
      "C:/Program Files/MPICH/SDK/Lib" 
      DOC "If a second mpi library is necessary, specify it here.")
    if(MPI_EXTRA_LIBRARIES)
      set(MPI_LIBRARIES ${MPI_LIBRARIES} ${MPI_EXTRA_LIBRARIES})
    endif(MPI_EXTRA_LIBRARIES)
  endif(NOT MPI_LIBRARIES)
  if(MPI_LIBRARIES)
    message(STATUS "Looking for MPI libraries... ${MPI_LIBRARIES}")
    # Test the MPI libraries
    set(MPICXX_PERFORM_TEST TRUE)
  else(MPI_LIBRARIES)
    message(STATUS "Looking for MPI libraries... FAILED")
  endif(MPI_LIBRARIES)
endif(MPI_MPICXX)  
# If we have what to test, do it now
if(MPICXX_PERFORM_TEST)
  # Create the MPITest directory
  set(MPITest_DIR ${PROJECT_BINARY_DIR}/MPITest)
  file(MAKE_DIRECTORY ${MPITest_DIR})
  # Create a CMakeLists.txt file which will generate the "mpicxxtest" executable
  if(MPI_MPICXX)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpicxxtest CXX)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_CXX_COMPILER ${MPI_MPICXX})\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")\n"
      "SET(CMAKE_CXX_FLAGS_RELEASE \"${CMAKE_CXX_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_CXX_FLAGS_DEBUG \"${CMAKE_CXX_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_CXX_FLAGS_RELWITHDEBUGINFO \"${CMAKE_CXX_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_CXX_FLAGS_MINSIZE \"${CMAKE_CXX_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(mpicxxtest mpicxxtest.cpp)\n")
  else(MPI_MPICXX)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpicxxtest CXX)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")\n"
      "SET(CMAKE_CXX_FLAGS_RELEASE \"${CMAKE_CXX_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_CXX_FLAGS_DEBUG \"${CMAKE_CXX_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_CXX_FLAGS_RELWITHDEBUGINFO \"${CMAKE_CXX_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_CXX_FLAGS_MINSIZE \"${CMAKE_CXX_FLAGS_MINSIZE}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpicxxtest mpictest.cpp)\n"
      "TARGET_LINK_LIBRARIES(mpicxxtest ${MPI_LIBRARIES})\n")
  endif(MPI_MPICXX)
  # Create a simple C++ source which only calls the MPI_Init and MPI_Finalize functions
  file(WRITE ${MPITest_DIR}/mpicxxtest.cpp
    "#include <mpi.h>\n"
    "int main(){\n"
    "int c;\n"
    "char **v;\n"
    "MPI_Init(&c, &v);\n"
    "MPI_Finalize();\n"
    "return(0);\n"
    "}\n")
  # Use TRY_COMPILE to make the target "mpictest"
  try_compile(MPITEST_OK ${MPITest_DIR} ${MPITest_DIR}
    mpictest OUTPUT_VARIABLE MY_OUTPUT)
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)
  # Process test result
  if(MPITEST_OK)
    message(STATUS "Trying to compile and link a simple MPI C++ program... OK")
    set(MPICXX_FOUND TRUE)
  else(MPITEST_OK)
    message(STATUS "Trying to compile and link a simple MPI C++ program... FAILED")
  endif(MPITEST_OK)
endif(MPICXX_PERFORM_TEST)

