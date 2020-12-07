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
# CMakeLists.txt file for the CUDA cuSPARSE SUNMatrix
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_CUSPARSE\n\")")

# Add variable cusparse_HEADERS with the exported cusparse header files
set(cusparse_HEADERS sunmatrix_cusparse.h)
add_prefix(${SUNDIALS_SOURCE_DIR}/include/sunmatrix/ cusparse_HEADERS)

# Add the library
sundials_add_library(sundials_sunmatrixcusparse
  SOURCES
    sunmatrix_cusparse.cu
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_sunmemcuda_obj
  LINK_LIBRARIES
    PUBLIC sundials_nveccuda ${CUDA_cusparse_LIBRARY}
    PRIVATE ${CUDA_cusolver_LIBRARY}
  OUTPUT_NAME
    sundials_sunmatrixcusparse
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

# Install the CUDA NVector header files
install(FILES ${cusparse_HEADERS} DESTINATION include/sunmatrix)

#
message(STATUS "Added SUNMATRIX_CUSPARSE module")
