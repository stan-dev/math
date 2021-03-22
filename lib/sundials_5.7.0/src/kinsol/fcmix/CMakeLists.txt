# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the FKINSOL library
# ---------------------------------------------------------------

# Add variable fkinsol_SOURCES with the sources for the FKINSOL library
set(fkinsol_SOURCES
  fkinband.c
  fkinbbd.c
  fkindense.c
  fkinjtimes.c
  fkinnullmatrix.c
  fkinnulllinsol.c
  fkinpreco.c
  fkinsol.c
  fkinsparse.c
  )

sundials_add_library(sundials_fkinsol
  SOURCES
    ${fkinsol_SOURCES}
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
    PRIVATE sundials_kinsol
  INCLUDE_DIRECTORIES
    PRIVATE ../
  OUTPUT_NAME
    sundials_fkinsol
  VERSION
    ${kinsollib_VERSION}
  SOVERSION
    ${kinsollib_SOVERSION}
  STATIC_ONLY
)

#
message(STATUS "Added KINSOL FCMIX module")
