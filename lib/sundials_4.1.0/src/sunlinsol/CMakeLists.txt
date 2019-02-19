# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
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

# Always add SUNDIALS provided linear solver modules
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(pcg)

add_subdirectory(spbcgs)
add_subdirectory(spfgmr)
add_subdirectory(spgmr)
add_subdirectory(sptfqmr)

# Add linear solver modules with TPL dependencies
if(KLU_ENABLE AND KLU_FOUND)
  add_subdirectory(klu)
endif(KLU_ENABLE AND KLU_FOUND)

if(SUPERLUMT_ENABLE AND SUPERLUMT_FOUND)
  add_subdirectory(superlumt)
endif(SUPERLUMT_ENABLE AND SUPERLUMT_FOUND)

if(LAPACK_ENABLE AND LAPACK_FOUND)
  add_subdirectory(lapackband)
  add_subdirectory(lapackdense)
endif(LAPACK_ENABLE AND LAPACK_FOUND)
