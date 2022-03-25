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
# CMakeLists.txt file for the F2003 band SUNMatrix object library
# ---------------------------------------------------------------

sundials_add_f2003_library(sundials_fsunmatrixband_mod
  SOURCES
    fsunmatrix_band_mod.f90 fsunmatrix_band_mod.c
  OBJECT_LIBRARIES
    sundials_fgeneric_mod_obj
  OUTPUT_NAME
    sundials_fsunmatrixband_mod
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_BAND F2003 interface")
