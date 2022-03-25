# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the RAJA NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_RAJA with ${SUNDIALS_RAJA_BACKENDS} backend(s)\n\")")

if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_sunmemlib sundials_sunmemcuda_obj)
  set(_lib_output_name sundials_nveccudaraja)
  set_source_files_properties(nvector_raja.cpp PROPERTIES LANGUAGE CUDA)
  set(_cxx_std cxx_std_11)
elseif(SUNDIALS_RAJA_BACKENDS MATCHES "HIP")
  set(_sunmemlib sundials_sunmemhip_obj)
  set(_lib_output_name sundials_nvechipraja)
  set(_hip_lib_or_not hip::device)
  set(_cxx_std cxx_std_11)
elseif(SUNDIALS_RAJA_BACKENDS MATCHES "SYCL")
  set(_sunmemlib sundials_sunmemsycl_obj)
  set(_lib_output_name sundials_nvecsyclraja)
  set(_cxx_std cxx_std_17)
endif()

# Create the library
sundials_add_library(sundials_nvecraja
  SOURCES
    nvector_raja.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_raja.h
  INCLUDE_SUBDIR
    nvector
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
    PUBLIC ${_cxx_std}
  OUTPUT_NAME
    ${_lib_output_name}
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_RAJA module with ${SUNDIALS_RAJA_BACKENDS} backend")
