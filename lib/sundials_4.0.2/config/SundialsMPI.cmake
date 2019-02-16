# ---------------------------------------------------------------------------
# Programmer: David J. Gardner @ LLNL
# ---------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------
# MPI tests for SUNDIALS CMake-based configuration.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Prior to CMake 3.10 the CMake FindMPI module considers:
#   1. Inspect MPI wrappers (MPI_<lang>_COMPILER)
#   2. Try guesses
#   3. Try the compiler (CMAKE_<lang>_COMPILER)
#
# Starting with CMake 3.10 the CMake FindMPI module considers:
#   1. Try the compiler (CMAKE_<lang>_COMPILER)
#   2. Inspect MPI wrappers (MPI_<lang>_COMPILER)
#   3. Try guesses
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# If MPI_<lang>_COMPILER is not defined check if CMAKE_<lang>_COMPILER
# can compile a simple MPI program and set MPI_<lang>_COMPILER
# ---------------------------------------------------------------------------

# check C compiler
if(NOT MPI_C_COMPILER)

  message(STATUS "Check for working MPI C compiler: ${CMAKE_C_COMPILER}")

  # Create the MPITest directory
  set(MPITest_DIR ${PROJECT_BINARY_DIR}/MPICCTest)
  file(MAKE_DIRECTORY ${MPITest_DIR})

  # Create a CMakeLists.txt file which will generate the "mpictest" executable
  file(WRITE ${MPITest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)\n"
    "PROJECT(mpictest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(mpictest mpictest.c)\n")

  # Create a simple C source which only calls MPI_Init and MPI_Finalize
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
  try_compile(MPI_TEST_OK ${MPITest_DIR} ${MPITest_DIR}
    mpictest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)

  # Process test result
  if(MPI_TEST_OK)
    message(STATUS "Check for working MPI C compiler: ${CMAKE_C_COMPILER} -- works")
    show_variable(MPI_C_COMPILER STRING "MPI C compiler" ${CMAKE_C_COMPILER})
    set(MPI_C_FOUND TRUE)
  else()
    message(STATUS "Check for working MPI C compiler: ${CMAKE_C_COMPILER} -- broken")
  endif()

endif()

# only check C++ and Fortran compilers if MPI C compiler was found and works
if(MPI_C_FOUND)

  # check CXX compiler
  if(CXX_FOUND AND (NOT MPI_CXX_COMPILER))

    message(STATUS "Check for working MPI C++ compiler: ${CMAKE_CXX_COMPILER}")

    # Create the MPITest directory
    set(MPITest_DIR ${PROJECT_BINARY_DIR}/MPICXXTest)
    file(MAKE_DIRECTORY ${MPITest_DIR})

    # Create a CMakeLists.txt file which will generate the "mpicxxtest" executable
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)\n"
      "PROJECT(mpicxxtest CXX)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_CXX_COMPILER \"${CMAKE_CXX_COMPILER}\")\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")\n"
      "SET(CMAKE_CXX_FLAGS_RELEASE \"${CMAKE_CXX_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_CXX_FLAGS_DEBUG \"${CMAKE_CXX_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_CXX_FLAGS_RELWITHDEBUGINFO \"${CMAKE_CXX_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_CXX_FLAGS_MINSIZE \"${CMAKE_CXX_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(mpicxxtest mpicxxtest.cpp)\n")

    # Create a simple C++ source which only calls MPI_Init and MPI_Finalize
    file(WRITE ${MPITest_DIR}/mpicxxtest.cpp
      "#include <mpi.h>\n"
      "int main(){\n"
      "int c;\n"
      "char **v;\n"
      "MPI_Init(&c, &v);\n"
      "MPI_Finalize();\n"
      "return(0);\n"
      "}\n")

    # Use TRY_COMPILE to make the target "mpicxxtest"
    try_compile(MPI_TEST_OK ${MPITest_DIR} ${MPITest_DIR}
      mpicxxtest OUTPUT_VARIABLE MY_OUTPUT)

    # To ensure we do not use stuff from the previous attempts,
    # we must remove the CMakeFiles directory.
    file(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)

    # Process test result
    if(MPI_TEST_OK)
      message(STATUS "Check for working MPI C++ compiler: ${CMAKE_CXX_COMPILER} -- works")
      show_variable(MPI_CXX_COMPILER STRING "MPI C++ compiler" ${CMAKE_CXX_COMPILER})
      set(MPI_CXX_FOUND TRUE)
    else()
      message(STATUS "Check for working MPI C++ compiler: ${CMAKE_CXX_COMPILER} -- broken")
    endif()

  endif()

  # check Fortran compiler
  if((F77_FOUND OR F90_FOUND) AND (NOT MPI_Fortran_COMPILER))

    message(STATUS "Check for working MPI Fortran compiler: ${CMAKE_Fortran_COMPILER}")

    # Create the MPI Test directory
    set(MPITest_DIR ${PROJECT_BINARY_DIR}/MPIFTest)
    file(MAKE_DIRECTORY ${MPITest_DIR})

    # Create a CMakeLists.txt file which will generate the "mpiftest" executable
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)\n"
      "PROJECT(mpiftest Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_Fortran_COMPILER \"${CMAKE_Fortran_COMPILER}\")\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n"
      "SET(CMAKE_Fortran_FLAGS_RELEASE \"${CMAKE_Fortran_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_Fortran_FLAGS_DEBUG \"${CMAKE_Fortran_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO \"${CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_Fortran_FLAGS_MINSIZE \"${CMAKE_Fortran_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(mpif9test mpiftest.f)\n")

    # Create a simple Fortran source which only calls MPI_Init and MPI_Finalize
    file(WRITE ${MPITest_DIR}/mpiftest.f
      "       INCLUDE \"mpif.h\"\n"
      "       INTEGER IER\n"
      "       CALL MPI_INIT(IER)\n"
      "       CALL MPI_FINALIZE(IER)\n"
      "       STOP\n"
      "       END\n")

    # Use TRY_COMPILE to make the target "mpiftest"
    try_compile(MPI_TEST_OK ${MPITest_DIR} ${MPITest_DIR}
      mpiftest OUTPUT_VARIABLE MY_OUTPUT)

    # To ensure we do not use stuff from the previous attempts,
    # we must remove the CMakeFiles directory.
    file(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)

    # Process test result
    if(MPI_TEST_OK)
      message(STATUS "Check for working MPI Fortran compiler: ${CMAKE_Fortran_COMPILER} -- works")
      show_variable(MPI_Fortran_COMPILER STRING "MPI Fortran compiler" ${CMAKE_Fortran_COMPILER})
      set(MPI_Fortran_FOUND TRUE)
    else()
      message(STATUS "Check for working MPI Fortran compiler: ${CMAKE_Fortran_COMPILER} -- broken")
    endif()

  endif()

endif()

# ---------------------------------------------------------------------------
# If MPI_<lang>_COMPILER is set, FindMPI will try to set the below variables
# for the given compiler wrapper. If MPI_<lang>_COMPILER is unset FindMPI
# will attempt to locate an installed MPI library and set the below
# variables.
#
#   MPI_<lang>_FOUND           TRUE if FindMPI found MPI flags for <lang>
#   MPI_<lang>_COMPILER        MPI Compiler wrapper for <lang>
#   MPI_<lang>_COMPILE_FLAGS   Compilation flags for MPI programs
#   MPI_<lang>_INCLUDE_PATH    Include path(s) for MPI header
#   MPI_<lang>_LINK_FLAGS      Linking flags for MPI programs
#   MPI_<lang>_LIBRARIES       All libraries to link MPI programs against
#
#   MPIEXEC_EXECUTABLE         Executable for running MPI programs
#   MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC_EXECUTABLE before
#                              giving it the number of processors to run on
#   MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC_EXECUTABLE directly
#                              before the executable to run.
#   MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC_EXECUTABLE after
#                              other flags
# ---------------------------------------------------------------------------

# Copy value of MPIEXEC_EXECUTABLE to MPIEXEC for older versions of CMake
if((CMAKE_VERSION VERSION_LESS 3.10) AND (MPIEXEC_EXECUTABLE))
  force_variable(MPIEXEC FILEPATH "MPI run command" ${MPIEXEC_EXECUTABLE})
endif()

find_package(MPI)

# Copy value of MPIEXEC to MPIEXEC_EXECUTABLE for older versions of CMake
if(CMAKE_VERSION VERSION_LESS 3.10)
  force_variable(MPIEXEC_EXECUTABLE FILEPATH "MPI run command" ${MPIEXEC})
  mark_as_advanced(MPIEXEC)
endif()

# MPI not functioning
if(NOT MPI_C_FOUND)
  set(MPI_C_FOUND FALSE)
  set(MPI_CXX_FOUND FALSE)
  set(MPI_Fortran_FOUND FALSE)
endif()

# show some advaned MPI C variables
mark_as_advanced(CLEAR MPI_C_COMPILER)
mark_as_advanced(CLEAR MPIEXEC_EXECUTABLE)

# hide some MPI C variables
mark_as_advanced(MPI_C_LIBRARIES)
mark_as_advanced(MPI_C_COMPILE_FLAGS)
mark_as_advanced(MPI_C_INCLUDE_PATH)
mark_as_advanced(MPI_C_LIBRARIES)
mark_as_advanced(MPI_C_LINK_FLAGS)

# hide some MPI variables
mark_as_advanced(MPI_EXTRA_LIBRARY)
mark_as_advanced(MPI_LIBRARY)

if(CXX_FOUND)
  # show some advaned MPI C variables
  mark_as_advanced(CLEAR MPI_CXX_COMPILER)
  # hide some MPI CXX variables
  mark_as_advanced(MPI_CXX_LIBRARIES)
  mark_as_advanced(MPI_CXX_COMPILE_FLAGS)
  mark_as_advanced(MPI_CXX_INCLUDE_PATH)
  mark_as_advanced(MPI_CXX_LIBRARIES)
  mark_as_advanced(MPI_CXX_LINK_FLAGS)
endif()

if(F77_FOUND OR F90_FOUND)
  # show some advaned MPI Fortran variables
  mark_as_advanced(CLEAR MPI_Fortran_COMPILER)
  # hide some MPI Fortran variables
  mark_as_advanced(MPI_Fortran_COMPILE_FLAGS)
  mark_as_advanced(MPI_Fortran_INCLUDE_PATH)
  mark_as_advanced(MPI_Fortran_LIBRARIES)
  mark_as_advanced(MPI_Fortran_LINK_FLAGS)
endif()

# determine if MPI-2 is supported
if(MPI_C_FOUND)

  # MPI_VERSION is set by FindMPI in CMake 3.10 and later, update to:
  # if(NOT MPI_VERSION) test else() check version number endif()

  # Create the MPITest directory
  set(MPITest_DIR ${PROJECT_BINARY_DIR}/MPI2Test)
  file(MAKE_DIRECTORY ${MPITest_DIR})

  # Create CMakeLists.txt file for "mpi2test" executable
  if(MPI_C_COMPILER)
    
    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)\n"
      "PROJECT(mpi2test C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(mpi2test mpi2test.c)\n")

  else()

    file(WRITE ${MPITest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)\n"
      "PROJECT(mpi2test C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpi2test mpi2test.c)\n"
      "TARGET_LINK_LIBRARIES(mpi2test ${MPI_LIBRARIES})\n")

  endif()

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
  else()
    message(STATUS "Checking for MPI-2 support... FAILED")
    set(MPIC_MPI2 FALSE)
  endif()

endif()
