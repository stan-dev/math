# ------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------
# CMakeLists.txt file for the F2003 SPGMR SUNLinearSolver object library
# ------------------------------------------------------------------------

sundials_add_f2003_library(sundials_fsunlinsolspgmr_mod
  SOURCES
    fsunlinsol_spgmr_mod.f90 fsunlinsol_spgmr_mod.c
  OBJECT_LIBRARIES
    sundials_fgeneric_mod_obj
  OUTPUT_NAME
    sundials_fsunlinsolspgmr_mod
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_SOVERSION}
)

message(STATUS "Added SUNLINSOL_SPGMR F2003 interface")
