# ---------------------------------------------------------------
# Programmer:  David Gardner and Slaven Peles @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# examples/cvodes level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# C examples
IF(EXAMPLES_ENABLE_C)
  ADD_SUBDIRECTORY(serial)
  IF(MPI_ENABLE AND MPI_C_FOUND)
    ADD_SUBDIRECTORY(parallel)
  ENDIF()
  IF(OPENMP_ENABLE AND OPENMP_FOUND)
    ADD_SUBDIRECTORY(C_openmp)
  ENDIF()
ENDIF()
