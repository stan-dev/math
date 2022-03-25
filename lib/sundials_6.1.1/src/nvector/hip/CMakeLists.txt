# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the hip NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_HIP\n\")")

# Create the library
sundials_add_library(sundials_nvechip
  SOURCES
    nvector_hip.hip.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_hip.h
  INCLUDE_SUBDIR
    nvector
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_sunmemhip_obj
  LINK_LIBRARIES
    PUBLIC hip::device
  OUTPUT_NAME
    sundials_nvechip
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_HIP module")
