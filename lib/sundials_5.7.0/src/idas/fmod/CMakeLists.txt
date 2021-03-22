# ---------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the F2003 IDAS object library

set(idas_SOURCES fidas_mod.f90 fidas_mod.c)

# Create the library
sundials_add_f2003_library(sundials_fidas_mod
  SOURCES
    ${idas_SOURCES}
  OBJECT_LIBRARIES
    sundials_fgeneric_mod_obj
    sundials_fnvecserial_mod_obj
    sundials_fsunmatrixband_mod_obj
    sundials_fsunmatrixdense_mod_obj
    sundials_fsunmatrixsparse_mod_obj
    sundials_fsunlinsolband_mod_obj
    sundials_fsunlinsoldense_mod_obj
    sundials_fsunlinsolspbcgs_mod_obj
    sundials_fsunlinsolspfgmr_mod_obj
    sundials_fsunlinsolspgmr_mod_obj
    sundials_fsunlinsolsptfqmr_mod_obj
    sundials_fsunlinsolpcg_mod_obj
    sundials_fsunnonlinsolnewton_mod_obj
    sundials_fsunnonlinsolfixedpoint_mod_obj
  OUTPUT_NAME
    sundials_fidas_mod
  VERSION
    ${idaslib_VERSION}
  SOVERSION
    ${idaslib_SOVERSION}
)
message(STATUS "Added IDAS F2003 interface")
