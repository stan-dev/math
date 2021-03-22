# ---------------------------------------------------------------
# Programmer(s): Radu Serban, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the serial NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_SERIAL\n\")")

# Create the sundials_nvecserial library
sundials_add_library(sundials_nvecserial
  SOURCES
    nvector_serial.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_serial.h
  INCLUDE_SUBDIR
    nvector
  OBJECT_LIBRARIES
    sundials_generic_obj
  OUTPUT_NAME
    sundials_nvecserial
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_SERIAL module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fnvecserial
    SOURCES
      fnvector_serial.c
    LINK_LIBRARIES
      PRIVATE sundials_nvecserial
    OUTPUT_NAME
      sundials_fnvecserial
    VERSION
      ${nveclib_VERSION}
    SOVERSION
      ${nveclib_SOVERSION}
  )
  message(STATUS "Added NVECTOR_SERIAL F77 interface")
endif()
