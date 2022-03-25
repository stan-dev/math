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
# CMakeLists.txt file for the CUDA cuSPARSE SUNMatrix
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_CUSPARSE\n\")")

# Add the library
sundials_add_library(sundials_sunmatrixcusparse
  SOURCES
    sunmatrix_cusparse.cu
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_cusparse.h
  INCLUDE_SUBDIR
    sunmatrix
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_sunmemcuda_obj
  LINK_LIBRARIES
    PUBLIC ${CUDA_cusparse_LIBRARY}
    PRIVATE ${CUDA_cusolver_LIBRARY}
  OUTPUT_NAME
    sundials_sunmatrixcusparse
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_CUSPARSE module")
