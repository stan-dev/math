# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
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
  if(OPENMP_ENABLE AND OPENMP_FOUND)
    add_subdirectory(C_openmp)
  endif()
  if(OPENMP_DEVICE_ENABLE AND OPENMP_FOUND AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    add_subdirectory(C_openmpdev)
  endif()
  if(MPI_ENABLE AND MPI_C_FOUND)
    add_subdirectory(C_parallel)
  endif()
  if(HYPRE_ENABLE AND HYPRE_FOUND)
    add_subdirectory(C_parhyp)
  endif()
  if(PETSC_ENABLE AND PETSC_FOUND)
    add_subdirectory(C_petsc)
  endif()
endif()

# C++ examples
if(EXAMPLES_ENABLE_CXX AND CXX_FOUND)
  add_subdirectory(CXX_serial)
  if(MPI_ENABLE AND MPI_CXX_FOUND)
    add_subdirectory(CXX_parallel)
  endif()
  if(MPI_ENABLE AND MPI_CXX_FOUND AND HYPRE_ENABLE AND HYPRE_FOUND)
    add_subdirectory(CXX_parhyp)
  endif()
  if(OPENMP_ENABLE AND OPENMP_FOUND AND MPI_CXX_FOUND AND SUPERLUDIST_ENABLE)
    add_subdirectory(CXX_superludist)
  endif()
endif()

# F77 examples
if(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F77 AND F77_FOUND)
  add_subdirectory(F77_serial)
  if(MPI_ENABLE AND MPI_Fortran_FOUND)
    add_subdirectory(F77_parallel)
  endif()
endif()

# F90 examples
if(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F90 AND F90_FOUND)
  add_subdirectory(F90_serial)
  if(MPI_ENABLE AND MPI_Fortran_FOUND)
    # currently the only F90_parallel example
    # does not support anything but double precision
    if(SUNDIALS_PRECISION MATCHES "DOUBLE")
      add_subdirectory(F90_parallel)
    endif()
  endif()
endif()

# F2003 examples
if(F2003_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F2003 AND F2003_FOUND)
  add_subdirectory(F2003_serial)
  add_subdirectory(F2003_custom)
endif()