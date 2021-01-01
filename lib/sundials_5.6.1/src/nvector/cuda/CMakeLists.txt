# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the cuda NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_CUDA\n\")")

# Add variable nveccuda_HEADERS with the exported NVECUDA header files
set(nveccuda_HEADERS nvector_cuda.h)
add_prefix(${SUNDIALS_SOURCE_DIR}/include/nvector/ nveccuda_HEADERS)

# Create the library
sundials_add_library(sundials_nveccuda
  SOURCES
    nvector_cuda.cu
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_sunmemcuda_obj
  OUTPUT_NAME
    sundials_nveccuda
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

# Install the CUDA NVector header files
install(FILES ${nveccuda_HEADERS} DESTINATION include/nvector)

message(STATUS "Added NVECTOR_CUDA module")
