# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
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
# examples/sunlinsol level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} caliper)
endif()

# Always add the serial sunlinearsolver dense and band examples
add_subdirectory(band)
add_subdirectory(dense)

# Always add serial sunlinearsolver iterative examples
add_subdirectory(spgmr/serial)
add_subdirectory(spfgmr/serial)
add_subdirectory(spbcgs/serial)
add_subdirectory(sptfqmr/serial)
add_subdirectory(pcg/serial)

if(ENABLE_MPI AND MPI_C_FOUND)
  add_subdirectory(spgmr/parallel)
  add_subdirectory(spfgmr/parallel)
  add_subdirectory(spbcgs/parallel)
  add_subdirectory(sptfqmr/parallel)
endif()

if(BUILD_SUNMATRIX_CUSPARSE)
  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cusolversp)
  endif()
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

if(BUILD_SUNLINSOL_SUPERLUMT)
  add_subdirectory(superlumt)
endif()

if(BUILD_SUNLINSOL_SUPERLUDIST)
  add_subdirectory(superludist)
endif()
