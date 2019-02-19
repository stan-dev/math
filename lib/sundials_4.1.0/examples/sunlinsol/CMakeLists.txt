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
# examples/sunlinsol level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# Always add the serial sunlinearsolver dense/band/spils examples
ADD_SUBDIRECTORY(band)
ADD_SUBDIRECTORY(dense)

IF(KLU_ENABLE AND KLU_FOUND)
  ADD_SUBDIRECTORY(klu)
ENDIF(KLU_ENABLE AND KLU_FOUND)

IF(SUPERLUMT_ENABLE AND SUPERLUMT_FOUND)
  ADD_SUBDIRECTORY(superlumt)
ENDIF(SUPERLUMT_ENABLE AND SUPERLUMT_FOUND)

IF(LAPACK_ENABLE AND LAPACK_FOUND)
  ADD_SUBDIRECTORY(lapackband)
  ADD_SUBDIRECTORY(lapackdense)
ENDIF(LAPACK_ENABLE AND LAPACK_FOUND)

# Always add serial sunlinearsolver iterative examples
ADD_SUBDIRECTORY(spgmr/serial)
ADD_SUBDIRECTORY(spfgmr/serial)
ADD_SUBDIRECTORY(spbcgs/serial)
ADD_SUBDIRECTORY(sptfqmr/serial)
ADD_SUBDIRECTORY(pcg/serial)

IF(MPI_ENABLE AND MPI_C_FOUND)
  ADD_SUBDIRECTORY(spgmr/parallel)
  ADD_SUBDIRECTORY(spfgmr/parallel)
  ADD_SUBDIRECTORY(spbcgs/parallel)
  ADD_SUBDIRECTORY(sptfqmr/parallel)
  #ADD_SUBDIRECTORY(pcg/parallel)
ENDIF()
