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

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_MAGMADENSE with ${SUNDIALS_MAGMA_BACKENDS} backend(s)\n\")")

if(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA")
  set_source_files_properties(sunmatrix_magmadense.cpp PROPERTIES LANGUAGE CUDA)
  set(_libs_needed sundials_nveccuda ${CUDA_CUBLAS_LIBRARIES})
elseif(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP")
  set_source_files_properties(sunmatrix_magmadense.cpp PROPERTIES LANGUAGE CXX)
  set(_libs_needed sundials_nvechip hip::device)
endif()

# Add the sunmatrix_magmadense library
sundials_add_library(sundials_sunmatrixmagmadense
  SOURCES
    sunmatrix_magmadense.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_magmadense.h
  INCLUDE_SUBDIR
    sunmatrix
  OBJECT_LIBRARIES
    sundials_generic_obj
  INCLUDE_DIRECTORIES
    PUBLIC ${MAGMA_INCLUDE_DIR}
  LINK_LIBRARIES
    PUBLIC SUNDIALS::MAGMA ${_libs_needed}
  OUTPUT_NAME
    sundials_sunmatrixmagmadense
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_MAGMADENSE module with ${SUNDIALS_MAGMA_BACKENDS} backend(s)")
