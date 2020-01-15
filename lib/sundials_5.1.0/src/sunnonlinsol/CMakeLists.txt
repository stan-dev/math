# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
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

# Always add SUNDIALS provided nonlinear solver modules
add_subdirectory(newton)
add_subdirectory(fixedpoint)

if(PETSC_ENABLE AND PETSC_FOUND)
  add_subdirectory(petscsnes)
endif(PETSC_ENABLE AND PETSC_FOUND)
