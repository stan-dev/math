# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# CMakeLists.txt file for the ARKODE + XBraid interface library
# ---------------------------------------------------------------

# Set install message
install(CODE "MESSAGE(\"\nInstall ARKODE XBraid interface\n\")")

# Create the sundials_arkode_xbraid library
sundials_add_library(sundials_arkode_xbraid
  SOURCES
    arkode_xbraid.c
    ${SUNDIALS_SOURCE_DIR}/src/sundials/sundials_xbraid.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/arkode/arkode_xbraid.h
  INCLUDE_SUBDIR
    arkode
  LINK_LIBRARIES
    PUBLIC sundials_arkode SUNDIALS::XBRAID MPI::MPI_C
  INCLUDE_DIRECTORIES
    PRIVATE ../
  OUTPUT_NAME
    sundials_arkode_xbraid
  VERSION
    ${arkodelib_VERSION}
  SOVERSION
    ${arkodelib_SOVERSION}
)

# Finished
message(STATUS "Added ARKODE Xbraid interface")
