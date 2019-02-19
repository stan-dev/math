# ---------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the F2003 PCG SUNLinearSolver object library

set(sunlinsolpcg_SOURCES fsunlinsol_pcg.f90)

if(BUILD_STATIC_LIBS)
  sundials_add_f2003_interface_library(sundials_fsunlinsolpcg_mod_static STATIC
    ${sunlinsolpcg_SOURCES}
  )
  set_target_properties(sundials_fsunlinsolpcg_mod_static
    PROPERTIES OUTPUT_NAME sundials_fsunlinsolpcg_mod CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_fsunlinsolpcg_mod_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

if(BUILD_SHARED_LIBS)
  sundials_add_f2003_interface_library(sundials_fsunlinsolpcg_mod_shared SHARED
    ${sunlinsolpcg_SOURCES}
  )
  set_target_properties(sundials_fsunlinsolpcg_mod_shared
    PROPERTIES OUTPUT_NAME sundials_fsunlinsolpcg_mod CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_fsunlinsolpcg_mod_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_fsunlinsolpcg_mod_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

message(STATUS "Added SUNLINSOL_PCG F2003 interface")
