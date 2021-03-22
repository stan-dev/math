# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# examples/arkode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C examples
if(EXAMPLES_ENABLE_C)
  add_subdirectory(C_serial)
  add_subdirectory(C_manyvector)
  if(ENABLE_OPENMP AND OPENMP_FOUND)
    add_subdirectory(C_openmp)
  endif()
  if(ENABLE_OPENMP_DEVICE AND OPENMP_FOUND AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    add_subdirectory(C_openmpdev)
  endif()
  if(ENABLE_MPI AND MPI_C_FOUND)
    add_subdirectory(C_parallel)
  endif()
  if(ENABLE_HYPRE AND HYPRE_FOUND)
    add_subdirectory(C_parhyp)
  endif()
  if(ENABLE_PETSC AND PETSC_FOUND)
    add_subdirectory(C_petsc)
  endif()
endif()

# C++ examples
if(EXAMPLES_ENABLE_CXX)
  add_subdirectory(CXX_serial)
  if(ENABLE_MPI AND MPI_CXX_FOUND)
    add_subdirectory(CXX_parallel)
  endif()
  if(ENABLE_MPI AND MPI_CXX_FOUND AND ENABLE_HYPRE AND HYPRE_FOUND)
    add_subdirectory(CXX_parhyp)
  endif()
  if(MPI_CXX_FOUND AND ENABLE_SUPERLUDIST)
    add_subdirectory(CXX_superludist)
  endif()
  if(ENABLE_XBRAID)
    add_subdirectory(CXX_xbraid)
  endif()
endif()

# F77 examples
if(BUILD_FORTRAN77_INTERFACE AND EXAMPLES_ENABLE_F77)
  add_subdirectory(F77_serial)
  if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_subdirectory(F77_parallel)
  endif()
endif()

# F90 examples
if(BUILD_FORTRAN77_INTERFACE AND EXAMPLES_ENABLE_F90)
  add_subdirectory(F90_serial)
  if(ENABLE_MPI AND MPI_Fortran_FOUND)
    # currently the only F90_parallel example
    # does not support anything but double precision
    if(SUNDIALS_PRECISION MATCHES "DOUBLE")
      add_subdirectory(F90_parallel)
    endif()
  endif()
endif()

# F2003 examples
if(BUILD_FORTRAN_MODULE_INTERFACE AND EXAMPLES_ENABLE_F2003)
  add_subdirectory(F2003_serial)
  add_subdirectory(F2003_custom)
  if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_subdirectory(F2003_parallel)
  endif()
endif()
