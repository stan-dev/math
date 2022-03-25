# ----------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ----------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
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

sundials_add_f2003_library(sundials_fsunlinsoldense_mod
  SOURCES
    fsunlinsol_dense_mod.f90 fsunlinsol_dense_mod.c
  OBJECT_LIBRARIES
    sundials_fgeneric_mod_obj
  LINK_LIBRARIES
    PUBLIC sundials_fsunmatrixdense_mod
  OUTPUT_NAME
    sundials_fsunlinsoldense_mod
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_SOVERSION}
)
message(STATUS "Added SUNLINSOL_DENSE F2003 interface")
