# ---------------------------------------------------------------
# Programmer(s): Steve Smith and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
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
