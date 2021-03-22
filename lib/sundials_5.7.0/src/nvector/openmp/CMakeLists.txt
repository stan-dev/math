# ---------------------------------------------------------------
# Programmer(s): Steve Smith and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the openmp NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_OPENMP\n\")")

# Create the library
sundials_add_library(sundials_nvecopenmp
  SOURCES
    nvector_openmp.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_openmp.h
  INCLUDE_SUBDIR
    nvector
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC OpenMP::OpenMP_C
  OUTPUT_NAME
    sundials_nvecopenmp
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_OPENMP module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fnvecopenmp
    SOURCES
      fnvector_openmp.c
    LINK_LIBRARIES
      PRIVATE sundials_nvecopenmp
    OUTPUT_NAME
      sundials_fnvecopenmp
    VERSION
      ${nveclib_VERSION}
    SOVERSION
      ${nveclib_SOVERSION}
  )
  message(STATUS "Added NVECTOR_OPENMP F77 interface")
endif()
