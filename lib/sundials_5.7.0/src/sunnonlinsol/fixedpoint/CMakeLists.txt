# ------------------------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the fixed-point SUNNonlinearSolver library
# ------------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNNONLINSOL_FIXEDPOINT\n\")")

# Add the library
sundials_add_library(sundials_sunnonlinsolfixedpoint
  SOURCES
    sunnonlinsol_fixedpoint.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunnonlinsol/sunnonlinsol_fixedpoint.h
  INCLUDE_SUBDIR
    sunnonlinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
  OUTPUT_NAME
    sundials_sunnonlinsolfixedpoint
  VERSION
    ${sunnonlinsollib_VERSION}
  SOVERSION
    ${sunnonlinsollib_VERSION}
)

message(STATUS "Added SUNNONLINSOL_FIXEDPOINT module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# If FCMIX is enabled, build and install the F77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fsunnonlinsolfixedpoint
    SOURCES
      fsunnonlinsol_fixedpoint.c
    LINK_LIBRARIES
      PRIVATE sundials_fnvecserial sundials_sunnonlinsolfixedpoint
    OUTPUT_NAME
      sundials_fsunnonlinsolfixedpoint
    VERSION
      ${sunnonlinsol_VERSION}
    SOVERSION
      ${sunnonlinsol_SOVERSION}
  )
  message(STATUS "Added SUNNONLINSOL_FIXEDPOINT F77 interface")
endif()
