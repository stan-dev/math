# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU,
#                and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the SPBCGS SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SPBCGS\n\")")

# Exported header files
set(sunlinsolspbcgs_HEADERS
  ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_spbcgs.h)

# Add the sunlinsol_spbcgs library
sundials_add_library(sundials_sunlinsolspbcgs
  SOURCES
    sunlinsol_spbcgs.c
  OBJECT_LIBRARIES
    sundials_generic_obj
  OUTPUT_NAME
    sundials_sunlinsolspbcgs
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_VERSION}
)

# Install the header files
install(FILES ${sunlinsolspbcgs_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_SPBCGS module")

# Add F90 module if F2003 interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fsunlinsolspbcgs
    SOURCES
      fsunlinsol_spbcgs.c
    LINK_LIBRARIES
      PRIVATE sundials_sunlinsolspbcgs
    OUTPUT_NAME
      sundials_fsunlinsolspbcgs
    VERSION
      ${sunlinsollib_VERSION}
    SOVERSION
      ${sunlinsollib_SOVERSION}
  )
  message(STATUS "Added SUNLINSOL_SPBCGS F77 interface")
endif()
