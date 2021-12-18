# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the F2003 ARKODE object library
# ---------------------------------------------------------------

set(arkode_SOURCES
    farkode_mod.f90
    farkode_mod.c
    farkode_arkstep_mod.f90
    farkode_arkstep_mod.c
    farkode_erkstep_mod.f90
    farkode_erkstep_mod.c
    farkode_mristep_mod.f90
    farkode_mristep_mod.c)

# Create the library
sundials_add_f2003_library(sundials_farkode_mod
  SOURCES
    ${arkode_SOURCES}
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
    sundials_farkode_mod
  VERSION
    ${arkodelib_VERSION}
  SOVERSION
    ${arkodelib_SOVERSION}
)

message(STATUS "Added ARKODE F2003 interface")
