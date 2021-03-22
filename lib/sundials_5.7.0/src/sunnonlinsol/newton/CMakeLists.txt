# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# CMakeLists.txt file for the Newton SUNNonlinearSolver library
# ------------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNNONLINSOL_NEWTON\n\")")

# Add the library
sundials_add_library(sundials_sunnonlinsolnewton
  SOURCES
    sunnonlinsol_newton.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunnonlinsol/sunnonlinsol_newton.h
  INCLUDE_SUBDIR
    sunnonlinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
  OUTPUT_NAME
    sundials_sunnonlinsolnewton
  VERSION
    ${sunnonlinsollib_VERSION}
  SOVERSION
    ${sunnonlinsollib_VERSION}
)

message(STATUS "Added SUNNONLINSOL_NEWTON module")

# Add F2003 module if interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# If FCMIX is enabled, build and install the F77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fsunnonlinsolnewton
    SOURCES
      fsunnonlinsol_newton.c
    LINK_LIBRARIES
      PRIVATE sundials_fnvecserial sundials_sunnonlinsolnewton
    OUTPUT_NAME
      sundials_fsunnonlinsolnewton
    VERSION
      ${sunnonlinsol_VERSION}
    SOVERSION
      ${sunnonlinsol_SOVERSION}
  )
  message(STATUS "Added SUNNONLINSOL_NEWTON F77 interface")
endif()
