# ----------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ----------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ----------------------------------------------------------------------
# CMakeLists.txt file for the F2003 dense SUNLinearSolver object library
# ----------------------------------------------------------------------

set(sunlinsoldense_SOURCES fsunlinsol_dense_mod.f90 fsunlinsol_dense_mod.c)

if(BUILD_STATIC_LIBS)
  sundials_add_f2003_interface_library(sundials_fsunlinsoldense_mod_static STATIC_OBJECT
    ${sunlinsoldense_SOURCES}
  )
  add_dependencies(sundials_fsunlinsoldense_mod_static_obj sundials_fgeneric_static_obj)
  target_link_libraries(sundials_fsunlinsoldense_mod_static PUBLIC sundials_fsunmatrixdense_mod_static)
  set_target_properties(sundials_fsunlinsoldense_mod_static
    PROPERTIES OUTPUT_NAME sundials_fsunlinsoldense_mod CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_fsunlinsoldense_mod_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

if(BUILD_SHARED_LIBS)
  sundials_add_f2003_interface_library(sundials_fsunlinsoldense_mod_shared SHARED_OBJECT
    ${sunlinsoldense_SOURCES}
  )
  add_dependencies(sundials_fsunlinsoldense_mod_shared_obj sundials_fgeneric_shared_obj)
  target_link_libraries(sundials_fsunlinsoldense_mod_shared PUBLIC sundials_fsunmatrixdense_mod_shared)
  set_target_properties(sundials_fsunlinsoldense_mod_shared
    PROPERTIES OUTPUT_NAME sundials_fsunlinsoldense_mod CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_fsunlinsoldense_mod_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_fsunlinsoldense_mod_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

message(STATUS "Added SUNLINSOL_DENSE F2003 interface")
