# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
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

# required modules
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(sparse)

if(BUILD_SUNMATRIX_SLUNRLOC)
  add_subdirectory(slunrloc)
endif()

if(BUILD_SUNMATRIX_CUSPARSE)
  add_subdirectory(cusparse)
endif()

if(BUILD_SUNMATRIX_MAGMADENSE)
  add_subdirectory(magmadense)
endif()