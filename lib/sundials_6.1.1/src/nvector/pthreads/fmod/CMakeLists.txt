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
# CMakeLists.txt file for the F2003 Pthreads NVECTOR object library
# ---------------------------------------------------------------

set(nvecpthreads_SOURCES fnvector_pthreads_mod.f90 fnvector_pthreads_mod.c)

sundials_add_f2003_library(sundials_fnvecpthreads_mod
  SOURCES
    fnvector_pthreads_mod.f90 fnvector_pthreads_mod.c
  OBJECT_LIBRARIES
    sundials_fgeneric_mod_obj
  OUTPUT_NAME
    sundials_fnvecpthreads_mod
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_PTHREADS F2003 interface")
