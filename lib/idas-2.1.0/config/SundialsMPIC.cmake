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
# MPI-C tests for SUNDIALS CMake-based configuration.
#
# 

set(MPIC_FOUND FALSE)
set(MPIC_MPI2 FALSE)

# Local variable indicating whether to test MPI
set(MPIC_PERFORM_TEST FALSE)
# By default, we try to use the MPI compiler script
# Search for the MPICC compiler script
find_program(MPI_MPICC NAMES mpicc DOC "mpicc program")
if(MPI_MPICC)
  message(STATUS "Looking for MPI C compiler script... ${MPI_MPICC}")
  # Test the MPI compiler script
  set(MPIC_PERFORM_TEST TRUE)
else(MPI_MPICC)
  message(STATUS "Looking for MPI C compiler script... FAILED")
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
    set(MPIC_PERFORM_TEST TRUE)
  else(MPI_LIBRARIES)
    message(STATUS "Looking for MPI libraries... FAILED")
  endif(MPI_LIBRARIES)
endif(MPI_MPICC)  
# If we have what to test, do it now
if(MPIC_PERFORM_TEST)
  # Create the MPITest directory
  set(MPITest_DIR ${PROJECT_BINARY_DIR}/MPITest)
  file(MAKE_DIRECTORY ${MPITest_DIR})
  # Create a CMakeLists.txt file which will generate the "mpictest" executable
  if(MPI_MPICC)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpictest C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_COMPILER ${MPI_MPICC})\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(mpictest mpictest.c)\n")
  else(MPI_MPICC)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpictest C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpictest mpictest.c)\n"
      "TARGET_LINK_LIBRARIES(mpictest ${MPI_LIBRARIES})\n")
  endif(MPI_MPICC)
  # Create a simple C source which only calls the MPI_Init and MPI_Finalize functions
  file(WRITE ${MPITest_DIR}/mpictest.c
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
    message(STATUS "Trying to compile and link a simple MPI C program... OK")
    set(MPIC_FOUND TRUE)
  else(MPITEST_OK)
    message(STATUS "Trying to compile and link a simple MPI C program... FAILED")
  endif(MPITEST_OK)
endif(MPIC_PERFORM_TEST)
# Finally, if MPI-C was found and is working, 
# also check if it provides MPI-2 support
if(MPIC_FOUND)
  # Create a CMakeLists.txt file which will generate the "mpi2test" executable
  if(MPI_MPICC)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpi2test C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_COMPILER ${MPI_MPICC})\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(mpi2test mpi2test.c)\n")
  else(MPI_MPICC)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpi2test C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpi2test mpi2test.c)\n"
      "TARGET_LINK_LIBRARIES(mpi2test ${MPI_LIBRARIES})\n")
  endif(MPI_MPICC)
  # Create a simple C source which calls the MPI_Comm_f2c function
  file(WRITE ${MPITest_DIR}/mpi2test.c
    "#include <mpi.h>\n"
    "int main(){\n"
    "int c;\n"
    "char **v;\n"
    "MPI_Comm C_comm;\n"
    "MPI_Init(&c, &v);\n"
    "C_comm = MPI_Comm_f2c((MPI_Fint) 1);\n"
    "MPI_Finalize();\n"
    "return(0);\n"
    "}\n")
  # Use TRY_COMPILE to make the target "mpi2test"
  try_compile(MPITEST_OK ${MPITest_DIR} ${MPITest_DIR}
    mpi2test OUTPUT_VARIABLE MY_OUTPUT)
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  FILE(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)
  # Interpret test results
  if(MPITEST_OK)
    message(STATUS "Checking for MPI-2 support... OK")
    set(MPIC_MPI2 TRUE)
  else(MPITEST_OK)
    message(STATUS "Checking for MPI-2 support... FAILED")
  endif(MPITEST_OK)
endif(MPIC_FOUND)

