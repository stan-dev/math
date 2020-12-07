# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
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
# CMakeLists.txt file for the parallel NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_PARALLEL\n\")")

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Add variable nvecparallel_HEADERS with the exported NVECPARALLEL header files
set(nvecparallel_HEADERS
  ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_parallel.h
  )

# Create the library
sundials_add_library(sundials_nvecparallel
  SOURCES
    nvector_parallel.c
  OBJECT_LIBRARIES
    sundials_generic_obj
  OUTPUT_NAME
    sundials_nvecparallel
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

# Install the NVECPARALLEL header files
install(FILES ${nvecparallel_HEADERS} DESTINATION include/nvector)

#
message(STATUS "Added NVECTOR_PARALLEL module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fnvecparallel
    SOURCES
      fnvector_parallel.c
    LINK_LIBRARIES
      PRIVATE sundials_nvecparallel
    OUTPUT_NAME
      sundials_fnvecparallel
    VERSION
      ${nveclib_VERSION}
    SOVERSION
      ${nveclib_SOVERSION}
  )
  message(STATUS "Added NVECTOR_PARALLEL F77 interface")
endif()
