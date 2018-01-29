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
# MPI-Fortran90 tests for SUNDIALS CMake-based configuration.

set(MPIF90_FOUND FALSE)

# Local variable indicating whether to test MPI
set(MPIF90_PERFORM_TEST FALSE)
# By default, we try to use the MPI compiler script
# Search for the MPIF90 compiler script
find_program(MPI_MPIF90 NAMES mpif90 DOC "mpif90 program")
if(MPI_MPIF90)
  message(STATUS "Looking for MPI Fortran90 compiler script... ${MPI_MPIF90}")
  # Test the MPI compiler script
  set(MPIF90_PERFORM_TEST TRUE)
else(MPI_MPIF90)
  message(STATUS "Looking for MPI Fortran90 compiler script... FAILED")
  # If not already available, search for MPI headers and libraries.
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
    set(MPIF90_PERFORM_TEST TRUE)
  else(MPI_LIBRARIES)
    message(STATUS "Looking for MPI libraries... FAILED")
  endif(MPI_LIBRARIES)
endif(MPI_MPIF90)
# If we have what to test, do it now
if(MPIF90_PERFORM_TEST)
  # Create the MPITest directory
  set(MPITest_DIR ${PROJECT_BINARY_DIR}/MPITest)
  file(MAKE_DIRECTORY ${MPITest_DIR})
  # Create a CMakeLists.txt file which will generate the "mpif90test" executable
  if(MPI_MPIF90)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpif90test Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_Fortran_COMPILER ${MPI_MPIF90})\n"
      "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n"
      "ADD_EXECUTABLE(mpif90test mpif90test.f90)\n")
  else(MPI_MPIF90)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpif90test Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpif90test mpif90test.f90)\n"
      "TARGET_LINK_LIBRARIES(mpif90test ${MPI_LIBRARIES})\n")
  endif(MPI_MPIF90)
  # Create a simple F90 source which only calls the MPI_Init and MPI_Finalize functions
  file(WRITE ${MPITest_DIR}/mpif90test.f90
    "program test\n"
    "include \"mpif.h\"\n"
    "integer :: ier\n" 
    "call MPI_Init(ier)\n"
    "call MPI_Finalize(ier)\n"
    "stop\n"
    "end program\n")
  # Use TRY_COMPILE to make the target "mpif90test"
  try_compile(MPITEST_OK ${MPITest_DIR} ${MPITest_DIR}
    mpif90test OUTPUT_VARIABLE MY_OUTPUT)
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)
  # Process test result
  if(MPITEST_OK)
    message(STATUS "Trying to compile and link a simple MPI Fortran90 program... OK")
    set(MPIF90_FOUND TRUE)
  else(MPITEST_OK)
    message(STATUS "Trying to compile and link a simple MPI Fortran90 program... FAILED")
  endif(MPITEST_OK)
endif(MPIF90_PERFORM_TEST)