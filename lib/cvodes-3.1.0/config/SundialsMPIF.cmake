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
# MPI-Fortran tests for SUNDIALS CMake-based configuration.
#
# 

set(MPIF_FOUND FALSE)

# Local variable indicating whether to test MPI
set(MPIF_PERFORM_TEST FALSE)
# By default, we try to use the MPI compiler script
# Search for the MPIF77 compiler script
find_program(MPI_MPIF77 NAMES mpif77 DOC "mpif77 program")
if(MPI_MPIF77)
  message(STATUS "Looking for MPI Fortran compiler script... ${MPI_MPIF77}")
  # Test the MPI compiler script
  set(MPIF_PERFORM_TEST TRUE)
else(MPI_MPIF77)
  message(STATUS "Looking for MPI Fortran compiler script... FAILED")
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
    set(MPIF_PERFORM_TEST TRUE)
  else(MPI_LIBRARIES)
    message(STATUS "Looking for MPI libraries... FAILED")
  endif(MPI_LIBRARIES)
endif(MPI_MPIF77)
# If we have what to test, do it now
if(MPIF_PERFORM_TEST)
  # Create the MPITest directory
  set(MPITest_DIR ${PROJECT_BINARY_DIR}/MPITest)
  file(MAKE_DIRECTORY ${MPITest_DIR})
  # Create a CMakeLists.txt file which will generate the "mpiftest" executable
  if(MPI_MPIF77)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpiftest Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_Fortran_COMPILER ${MPI_MPIF77})\n"
      "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n"
      "ADD_EXECUTABLE(mpiftest mpiftest.f)\n")
  else(MPI_MPIF77)
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(mpiftest Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpiftest mpiftest.f)\n"
      "TARGET_LINK_LIBRARIES(mpiftest ${MPI_LIBRARIES})\n")
  endif(MPI_MPIF77)
  # Create a simple F77 source which only calls the MPI_Init and MPI_Finalize functions
  file(WRITE ${MPITest_DIR}/mpiftest.f
    "       INCLUDE \"mpif.h\"\n"
    "       INTEGER IER\n" 
    "       CALL MPI_INIT(IER)\n"
    "       CALL MPI_FINALIZE(IER)\n"
    "       STOP\n"
    "       END\n")
  # Use TRY_COMPILE to make the target "mpiftest"
  try_compile(MPITEST_OK ${MPITest_DIR} ${MPITest_DIR}
    mpiftest OUTPUT_VARIABLE MY_OUTPUT)
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)
  # Process test result
  if(MPITEST_OK)
    message(STATUS "Trying to compile and link a simple MPI Fortran program... OK")
    set(MPIF_FOUND TRUE)
  else(MPITEST_OK)
    message(STATUS "Trying to compile and link a simple MPI Fortran program... FAILED")
  endif(MPITEST_OK)
endif(MPIF_PERFORM_TEST)