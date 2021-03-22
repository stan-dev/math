# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the FIDA library
# ---------------------------------------------------------------

# Add variable fida_SOURCES with the sources for the FIDA library
set(fida_SOURCES
  fidaband.c
  fidabbd.c
  fida.c
  fidadense.c
  fidaewt.c
  fidajtimes.c
  fidanullmatrix.c
  fidanullnonlinsol.c
  fidapreco.c
  fidaroot.c
  fidasparse.c
  )

sundials_add_library(sundials_fida
  SOURCES
    ${fida_SOURCES}
  OBJECT_LIBRARIES
    sundials_fnvecserial_obj
    sundials_fsunmatrixband_obj
    sundials_fsunmatrixdense_obj
    sundials_fsunmatrixsparse_obj
    sundials_fsunlinsolband_obj
    sundials_fsunlinsoldense_obj
    sundials_fsunlinsolspbcgs_obj
    sundials_fsunlinsolspfgmr_obj
    sundials_fsunlinsolspgmr_obj
    sundials_fsunlinsolsptfqmr_obj
    sundials_fsunlinsolpcg_obj
    sundials_fsunnonlinsolnewton_obj
    sundials_fsunnonlinsolfixedpoint_obj
  LINK_LIBRARIES
    PRIVATE sundials_ida
  INCLUDE_DIRECTORIES
    PRIVATE ../
  OUTPUT_NAME
    sundials_fida
  VERSION
    ${idalib_VERSION}
  SOVERSION
    ${idalib_SOVERSION}
  STATIC_ONLY
)

#
message(STATUS "Added IDA FCMIX module")
