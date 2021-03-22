# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# CMakeLists.txt file for the FCVODE library
# ------------------------------------------------------------------------------

# Add variable fcvode_SOURCES with the sources for the FCVODE library
set(fcvode_SOURCES
  fcvband.c
  fcvbbd.c
  fcvbp.c
  fcvdense.c
  fcvewt.c
  fcvjtimes.c
  fcvnulllinsol.c
  fcvnullmatrix.c
  fcvnullnonlinsol.c
  fcvode.c
  fcvpreco.c
  fcvroot.c
  fcvsparse.c
  )

sundials_add_library(sundials_fcvode
  SOURCES
    ${fcvode_SOURCES}
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
    PRIVATE sundials_cvode
  INCLUDE_DIRECTORIES
    PRIVATE ../
  OUTPUT_NAME
    sundials_fcvode
  VERSION
    ${cvodelib_VERSION}
  SOVERSION
    ${cvodelib_SOVERSION}
  STATIC_ONLY
)

#
message(STATUS "Added CVODE FCMIX module")
