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

install(CODE "MESSAGE(\"\nInstall NVECTOR_RAJA with ${SUNDIALS_RAJA_BACKENDS} backend\n\")")

# Add variable nveccuda_HEADERS with the exported NVECRAJA header files
set(nvecraja_HEADERS nvector_raja.h)
add_prefix(${SUNDIALS_SOURCE_DIR}/include/nvector/ nvecraja_HEADERS)

if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_sunmemlib sundials_sunmemcuda_obj)
  # set(_compile_defs SUNDIALS_RAJA_BACKENDS_CUDA)
  set(_lib_output_name sundials_nveccudaraja)
  set_source_files_properties(nvector_raja.cpp PROPERTIES LANGUAGE CUDA)
elseif(SUNDIALS_RAJA_BACKENDS MATCHES "HIP")
  set(_sunmemlib sundials_sunmemhip_obj)
  set(_lib_output_name sundials_nvechipraja)
  # set(_compile_defs SUNDIALS_RAJA_BACKENDS_HIP)
  set(_hip_lib_or_not hip::device)
endif()

# Create the library
sundials_add_library(sundials_nvecraja
  SOURCES
    nvector_raja.cpp
  OBJECT_LIBRARIES
    sundials_generic_obj
    ${_sunmemlib}
  LINK_LIBRARIES
    PUBLIC ${_hip_lib_or_not} RAJA
  COMPILE_OPTIONS
    PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>
  COMPILE_DEFINITIONS
    PRIVATE ${_compile_defs}
  COMPILE_FEATURES
    PUBLIC cxx_std_11
  OUTPUT_NAME
    ${_lib_output_name}
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

# Install the NVECRAJA header files
install(FILES ${nvecraja_HEADERS} DESTINATION include/nvector)

message(STATUS "Added NVECTOR_RAJA module with ${SUNDIALS_RAJA_BACKENDS} backend")
