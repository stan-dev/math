# -----------------------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Find module for Trilinos that uses the TrilinosConfig.cmake that is installed
# with Trilinos. The module will also create a SUNDIALS::TRILINOS target.
# -----------------------------------------------------------------------------

# First try and find Trilinos using Trilinos_DIR only.
find_package(Trilinos
  NAMES Trilinos TRILINOS
  PATHS
    ${Trilinos_DIR}/lib/cmake/Trilinos
    ${Trilinos_DIR}
  NO_DEFAULT_PATH
  QUIET)

# set package variables including Trilinos_FOUND
find_package_handle_standard_args(Trilinos
  REQUIRED_VARS
    Trilinos_LIBRARIES      # defined in TrilinosConfig.cmake
    Trilinos_INCLUDE_DIRS   # defined in TrilinosConfig.cmake
  )

# Create Trilinos target
if(Trilinos_FOUND)

  if(NOT TARGET SUNDIALS::TRILINOS)
    add_library(SUNDIALS::TRILINOS IMPORTED INTERFACE)
  endif()

  set_target_properties(SUNDIALS::TRILINOS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Trilinos_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${Trilinos_LIBRARIES}")

endif()
