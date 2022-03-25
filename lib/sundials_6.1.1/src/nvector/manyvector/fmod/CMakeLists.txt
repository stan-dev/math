# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the F2003 manyvector NVECTOR object library
# ---------------------------------------------------------------

if(BUILD_NVECTOR_MANYVECTOR)
  sundials_add_f2003_library(sundials_fnvecmanyvector_mod
    SOURCES
      fnvector_manyvector_mod.f90 fnvector_manyvector_mod.c
    OBJECT_LIBRARIES
      sundials_fgeneric_mod_obj
    OUTPUT_NAME
      sundials_fnvecmanyvector_mod
    VERSION
      ${nveclib_VERSION}
    SOVERSION
      ${nveclib_SOVERSION}
  )
  message(STATUS "Added NVECTOR_MANYVECTOR F2003 Interface")
endif()

if(BUILD_NVECTOR_MPIMANYVECTOR)

  if(MPI_C_COMPILER)
    # use MPI wrapper as the compiler
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  elseif()
    # add MPI_INCLUDE_PATH to include directories
    include_directories(${MPI_INCLUDE_PATH})
  endif()
  if(MPI_Fortran_COMPILER)
    # use MPI wrapper as the compiler
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  endif()

  sundials_add_f2003_library(sundials_fnvecmpimanyvector_mod
    SOURCES
      fnvector_mpimanyvector_mod.f90 fnvector_mpimanyvector_mod.c
    OBJECT_LIBRARIES
      sundials_fgeneric_mod_obj
    OUTPUT_NAME
      sundials_fnvecmpimanyvector_mod
    VERSION
      ${nveclib_VERSION}
    SOVERSION
      ${nveclib_SOVERSION}
  )

  message(STATUS "Added NVECTOR_MPIMANYVECTOR F2003 Interface")
endif()
