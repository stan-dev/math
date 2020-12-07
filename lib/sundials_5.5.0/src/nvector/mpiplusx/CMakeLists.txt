# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the MPIPlusX NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_MPIPLUSX\n\")")

if(MPI_C_FOUND AND MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
elseif(MPI_C_FOUND)
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Add variable nvecmpiplusx_HEADERS with the exported NVECMPIPLUSX header files
set(nvecmpiplusx_HEADERS
  ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_mpimanyvector.h
  ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_mpiplusx.h
)

# Create the library
sundials_add_library(sundials_nvecmpiplusx
  SOURCES
    nvector_mpiplusx.c
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_nvecmpimanyvector_obj
  OUTPUT_NAME
    sundials_nvecmpiplusx
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

# Install the NVECMPIPLUSX header files
install(FILES ${nvecmpiplusx_HEADERS} DESTINATION include/nvector)

#
message(STATUS "Added NVECTOR_MPIPLUSX module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()
