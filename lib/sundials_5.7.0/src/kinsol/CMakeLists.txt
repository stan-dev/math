# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Radu Serban, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the KINSOL library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall KINSOL\n\")")

# Add variable kinsol_SOURCES with the sources for the KINSOL library
set(kinsol_SOURCES
  kinsol.c
  kinsol_bbdpre.c
  kinsol_direct.c
  kinsol_io.c
  kinsol_ls.c
  kinsol_spils.c
  )

# Add variable kinsol_HEADERS with the exported KINSOL header files
set(kinsol_HEADERS
  kinsol.h
  kinsol_bbdpre.h
  kinsol_direct.h
  kinsol_ls.h
  kinsol_spils.h
  )

# Add prefix with complete path to the KINSOL header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/kinsol/ kinsol_HEADERS)

# Create the library
sundials_add_library(sundials_kinsol
  SOURCES
    ${kinsol_SOURCES}
  HEADERS
    ${kinsol_HEADERS}
  INCLUDE_SUBDIR
    kinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_nvecserial_obj
    sundials_sunmatrixband_obj
    sundials_sunmatrixdense_obj
    sundials_sunmatrixsparse_obj
    sundials_sunlinsolband_obj
    sundials_sunlinsoldense_obj
    sundials_sunlinsolspbcgs_obj
    sundials_sunlinsolspfgmr_obj
    sundials_sunlinsolspgmr_obj
    sundials_sunlinsolsptfqmr_obj
    sundials_sunlinsolpcg_obj
  OUTPUT_NAME
    sundials_kinsol
  VERSION
    ${kinsollib_VERSION}
  SOVERSION
    ${kinsollib_SOVERSION}
)

# Finished KINSOL
message(STATUS "Added KINSOL module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

if(BUILD_FORTRAN77_INTERFACE)
  add_subdirectory(fcmix)
endif()
