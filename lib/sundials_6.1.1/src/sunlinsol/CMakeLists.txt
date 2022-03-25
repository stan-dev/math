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
# linear solver level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required native linear solvers
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(pcg)
add_subdirectory(spbcgs)
add_subdirectory(spfgmr)
add_subdirectory(spgmr)
add_subdirectory(sptfqmr)

# optional TPL linear solvers
if(BUILD_SUNLINSOL_CUSOLVERSP)
  add_subdirectory(cusolversp)
endif()

if(BUILD_SUNLINSOL_KLU)
  add_subdirectory(klu)
endif()

if(BUILD_SUNLINSOL_LAPACKBAND)
  add_subdirectory(lapackband)
endif()

if(BUILD_SUNLINSOL_LAPACKDENSE)
  add_subdirectory(lapackdense)
endif()

if(BUILD_SUNLINSOL_MAGMADENSE)
  add_subdirectory(magmadense)
endif()

if(BUILD_SUNLINSOL_ONEMKLDENSE)
  add_subdirectory(onemkldense)
endif()

if(BUILD_SUNLINSOL_SUPERLUDIST)
  add_subdirectory(superludist)
endif()

if(BUILD_SUNLINSOL_SUPERLUMT)
  add_subdirectory(superlumt)
endif()
