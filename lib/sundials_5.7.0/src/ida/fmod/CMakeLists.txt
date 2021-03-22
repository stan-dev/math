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
# CMakeLists.txt file for the F2003 IDA object library
# ---------------------------------------------------------------

set(ida_SOURCES fida_mod.f90 fida_mod.c)

# Create the library
sundials_add_f2003_library(sundials_fida_mod
  SOURCES
    ${ida_SOURCES}
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
    sundials_fida_mod
  VERSION
    ${idalib_VERSION}
  SOVERSION
    ${idalib_SOVERSION}
)

message(STATUS "Added IDA F2003 interface")
