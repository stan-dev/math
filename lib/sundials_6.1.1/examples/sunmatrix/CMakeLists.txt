# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# examples/sunmatrix level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} caliper)
endif()

# Always add the serial sunmatrix dense/band/sparse examples
add_subdirectory(dense)
add_subdirectory(band)
add_subdirectory(sparse)

if(BUILD_SUNMATRIX_CUSPARSE)
  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cusparse)
  endif()
endif()

if(BUILD_SUNMATRIX_MAGMADENSE)
  add_subdirectory(magmadense)
endif()

if(BUILD_SUNMATRIX_ONEMKLDENSE)
  add_subdirectory(onemkldense)
endif()

if(BUILD_SUNMATRIX_SLUNRLOC)
  add_subdirectory(slunrloc)
endif()
