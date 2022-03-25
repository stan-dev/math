# ---------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------
# Setup MPI for SUNDIALS CMake-based configuration.
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

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_MPI_INCLUDED)
  set(SUNDIALS_MPI_INCLUDED)
else()
  return()
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

mark_as_advanced(MPI_EXTRA_LIBRARY)
mark_as_advanced(MPI_LIBRARY)

foreach(lang ${_SUNDIALS_ENABLED_LANGS})
  mark_as_advanced(CLEAR MPI_${lang}_COMPILER)
  mark_as_advanced(MPI_${lang}_LIBRARIES)
  mark_as_advanced(MPI_${lang}_COMPILE_FLAGS)
  mark_as_advanced(MPI_${lang}_INCLUDE_PATH)
  mark_as_advanced(MPI_${lang}_LIBRARIES)
  mark_as_advanced(MPI_${lang}_LINK_FLAGS)
endforeach()

find_package(MPI 2.0.0 REQUIRED)

# ---------------------------------------------------------------------------
# Configure the presentation of MPI options in the GUI.
# ---------------------------------------------------------------------------

mark_as_advanced(CLEAR MPIEXEC_EXECUTABLE)

mark_as_advanced(MPI_EXTRA_LIBRARY)
mark_as_advanced(MPI_LIBRARY)

foreach(lang ${_SUNDIALS_ENABLED_LANGS})
  mark_as_advanced(CLEAR MPI_${lang}_COMPILER)
  mark_as_advanced(MPI_${lang}_LIBRARIES)
  mark_as_advanced(MPI_${lang}_COMPILE_FLAGS)
  mark_as_advanced(MPI_${lang}_INCLUDE_PATH)
  mark_as_advanced(MPI_${lang}_LIBRARIES)
  mark_as_advanced(MPI_${lang}_LINK_FLAGS)
endforeach()