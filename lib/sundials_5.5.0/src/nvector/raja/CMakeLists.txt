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
# CMakeLists.txt file for the RAJA NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_RAJA\n\")")

# Add variable nveccuda_HEADERS with the exported NVECRAJA header files
set(nvecraja_HEADERS nvector_raja.h)
add_prefix(${SUNDIALS_SOURCE_DIR}/include/nvector/ nvecraja_HEADERS)

# Create the library
sundials_add_library(sundials_nveccudaraja
  SOURCES
    nvector_raja.cu
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_sunmemcuda_obj
  LINK_LIBRARIES
    PUBLIC RAJA
  COMPILE_OPTIONS
    PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>
  COMPILE_FEATURES
    PUBLIC cxx_std_11
  OUTPUT_NAME
    sundials_nveccudaraja
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

# Install the NVECRAJA header files
install(FILES ${nvecraja_HEADERS} DESTINATION include/nvector)

message(STATUS "Added NVECTOR_RAJA module")
