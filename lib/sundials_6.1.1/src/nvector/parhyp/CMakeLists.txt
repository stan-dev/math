# ---------------------------------------------------------------
# Programmer(s): Radu Serban, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the parhyp NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_PARHYP\n\")")

# Create the library
sundials_add_library(sundials_nvecparhyp
  SOURCES
    nvector_parhyp.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_parhyp.h
  INCLUDE_SUBDIR
    nvector
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC
    SUNDIALS::HYPRE
    MPI::MPI_C
  OUTPUT_NAME
    sundials_nvecparhyp
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_PARHYP module")
