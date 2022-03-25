# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# nonlinear solver level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required modules
add_subdirectory(newton)
add_subdirectory(fixedpoint)

if(BUILD_SUNNONLINSOL_PETSCSNES)
  add_subdirectory(petscsnes)
endif()
