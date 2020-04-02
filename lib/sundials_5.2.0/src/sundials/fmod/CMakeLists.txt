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
# CMakeLists.txt file for the F2003 SUNDIALS object library
# ---------------------------------------------------------------

set(sundials_SOURCES
  fsundials_futils_mod.f90
  fsundials_futils_mod.c
  fsundials_linearsolver_mod.f90
  fsundials_linearsolver_mod.c
  fsundials_matrix_mod.f90
  fsundials_matrix_mod.c
  fsundials_nonlinearsolver_mod.f90
  fsundials_nonlinearsolver_mod.c
  fsundials_nvector_mod.f90
  fsundials_nvector_mod.c
  fsundials_types_mod.f90
  fsundials_types_mod.c)

if(BUILD_STATIC_LIBS)
  add_library(sundials_fgeneric_static_obj OBJECT ${sundials_SOURCES})
  set_target_properties(sundials_fgeneric_static_obj PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}_STATIC)
  target_include_directories(sundials_fgeneric_static_obj PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}_STATIC)
endif()

if(BUILD_SHARED_LIBS)
  add_library(sundials_fgeneric_shared_obj OBJECT ${sundials_SOURCES})
  set_target_properties(sundials_fgeneric_shared_obj PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}_SHARED)
  target_include_directories(sundials_fgeneric_shared_obj PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}_SHARED)
endif()
