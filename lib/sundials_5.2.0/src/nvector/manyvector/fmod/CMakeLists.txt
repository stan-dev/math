# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the F2003 manyvector NVECTOR object library
# ---------------------------------------------------------------

set(nvecmanyvector_SOURCES fnvector_manyvector_mod.f90 fnvector_manyvector_mod.c)
set(nvecmpimanyvector_SOURCES fnvector_mpimanyvector_mod.f90 fnvector_mpimanyvector_mod.c)

if(MPI_ENABLE AND MPI_C_FOUND)
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
endif()

if(BUILD_STATIC_LIBS)
  # ----------------------- No MPI
  sundials_add_f2003_interface_library(sundials_fnvecmanyvector_mod_static STATIC_OBJECT
                                       ${nvecmanyvector_SOURCES})
  add_dependencies(sundials_fnvecmanyvector_mod_static_obj sundials_fgeneric_static_obj)
  set_target_properties(sundials_fnvecmanyvector_mod_static PROPERTIES
                        OUTPUT_NAME sundials_fnvecmanyvector_mod
                        CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_fnvecmanyvector_mod_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- With MPI
  if(MPI_ENABLE AND MPI_C_FOUND)
    sundials_add_f2003_interface_library(sundials_fnvecmpimanyvector_mod_static STATIC_OBJECT
                                         ${nvecmpimanyvector_SOURCES})
    add_dependencies(sundials_fnvecmpimanyvector_mod_static_obj sundials_fgeneric_static_obj)
    set_target_properties(sundials_fnvecmpimanyvector_mod_static PROPERTIES
                          OUTPUT_NAME sundials_fnvecmpimanyvector_mod
                          CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fnvecmpimanyvector_mod_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
endif(BUILD_STATIC_LIBS)

if(BUILD_SHARED_LIBS)
  # ----------------------- No MPI
  sundials_add_f2003_interface_library(sundials_fnvecmanyvector_mod_shared SHARED_OBJECT
                                       ${nvecmanyvector_SOURCES})
  add_dependencies(sundials_fnvecmanyvector_mod_shared_obj sundials_fgeneric_shared_obj)
  set_target_properties(sundials_fnvecmanyvector_mod_shared PROPERTIES
                        OUTPUT_NAME sundials_fnvecmanyvector_mod
                        CLEAN_DIRECT_OUTPUT 1
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_fnvecmanyvector_mod_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- With MPI
  if(MPI_ENABLE AND MPI_C_FOUND)
    sundials_add_f2003_interface_library(sundials_fnvecmpimanyvector_mod_shared SHARED_OBJECT
                                         ${nvecmpimanyvector_SOURCES})
    add_dependencies(sundials_fnvecmpimanyvector_mod_shared_obj sundials_fgeneric_shared_obj)
    set_target_properties(sundials_fnvecmpimanyvector_mod_shared PROPERTIES
                          OUTPUT_NAME sundials_fnvecmpimanyvector_mod
                          CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fnvecmpimanyvector_mod_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
endif(BUILD_SHARED_LIBS)

message(STATUS "Added NVECTOR_MANYVECTOR F2003 Interface")

if(MPI_ENABLE AND MPI_C_FOUND)
  message(STATUS "Added NVECTOR_MPIMANYVECTOR F2003 Interface")
endif()
