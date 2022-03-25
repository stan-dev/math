# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# matrix level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required native matrices
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(sparse)

# optional TPL matrices
if(BUILD_SUNMATRIX_CUSPARSE)
  add_subdirectory(cusparse)
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
