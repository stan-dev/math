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
# CMakeLists.txt file for the FARKODE library
# ---------------------------------------------------------------

# Add variable farkode_SOURCES with the sources for the FARKODE library
set(farkode_SOURCES
  farkadapt.c
  farkband.c
  farkbandmass.c
  farkbbd.c
  farkbp.c
  farkdense.c
  farkdensemass.c
  farkewt.c
  farkexpstab.c
  farkjtimes.c
  farkmtimes.c
  farknullnonlinsol.c
  farknulllinsol.c
  farknullmatrix.c
  farkode.c
  farkpreco.c
  farkmasspreco.c
  farkroot.c
  farksparse.c
  farksparsemass.c
  )

sundials_add_library(sundials_farkode
  SOURCES
    ${farkode_SOURCES}
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
    PRIVATE sundials_arkode
  INCLUDE_DIRECTORIES
    PRIVATE ../
  OUTPUT_NAME
    sundials_farkode
  VERSION
    ${arkodelib_VERSION}
  SOVERSION
    ${arkodelib_SOVERSION}
  STATIC_ONLY
)

#
message(STATUS "Added ARKODE FCMIX module")
