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
# src level CMakeLists.txt for SUNDIALS (for cmake build system)
# ---------------------------------------------------------------

# Always add SUNDIALS provided serial modules
ADD_SUBDIRECTORY(sundials)
ADD_SUBDIRECTORY(nvec_ser)
ADD_SUBDIRECTORY(sunmat_dense)
ADD_SUBDIRECTORY(sunmat_band)
ADD_SUBDIRECTORY(sunmat_sparse)
ADD_SUBDIRECTORY(sunlinsol_band)
ADD_SUBDIRECTORY(sunlinsol_dense)

IF(KLU_FOUND)
  ADD_SUBDIRECTORY(sunlinsol_klu)
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  ADD_SUBDIRECTORY(sunlinsol_superlumt)
ENDIF(SUPERLUMT_FOUND)

IF(LAPACK_FOUND)
  ADD_SUBDIRECTORY(sunlinsol_lapackband)
  ADD_SUBDIRECTORY(sunlinsol_lapackdense)
ENDIF(LAPACK_FOUND)

# Always add SUNDIALS provided iterative linear solver modules
ADD_SUBDIRECTORY(sunlinsol_spgmr)
ADD_SUBDIRECTORY(sunlinsol_spfgmr)
ADD_SUBDIRECTORY(sunlinsol_spbcgs)
ADD_SUBDIRECTORY(sunlinsol_sptfqmr)
ADD_SUBDIRECTORY(sunlinsol_pcg)

# Add nonlinerar solver modules
ADD_SUBDIRECTORY(sunnonlinsol)

IF(MPI_C_FOUND)
  ADD_SUBDIRECTORY(nvec_par)
ENDIF()

IF(HYPRE_FOUND)
  ADD_SUBDIRECTORY(nvec_parhyp)
ENDIF(HYPRE_FOUND)

IF(OPENMP_FOUND)
  IF(OPENMP_ENABLE)
    ADD_SUBDIRECTORY(nvec_openmp)
  ENDIF()
  IF(OPENMP_DEVICE_ENABLE)
    ADD_SUBDIRECTORY(nvec_openmpdev)
  ENDIF()
ENDIF(OPENMP_FOUND)

IF(PTHREADS_FOUND)
  ADD_SUBDIRECTORY(nvec_pthreads)
ENDIF(PTHREADS_FOUND)

IF(PETSC_FOUND)
  ADD_SUBDIRECTORY(nvec_petsc)
ENDIF(PETSC_FOUND)

IF(CUDA_FOUND)
  ADD_SUBDIRECTORY(nvec_cuda)
ENDIF(CUDA_FOUND)

IF(RAJA_FOUND)
  ADD_SUBDIRECTORY(nvec_raja)
ENDIF(RAJA_FOUND)

# ARKODE library

IF(BUILD_ARKODE)
  ADD_SUBDIRECTORY(arkode)
  IF(F77_INTERFACE_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(arkode/fcmix)
  ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)
ENDIF(BUILD_ARKODE)

# CVODE library

IF(BUILD_CVODE)
  ADD_SUBDIRECTORY(cvode)
  IF(F77_INTERFACE_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(cvode/fcmix)
  ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)
ENDIF(BUILD_CVODE)

# CVODES library

IF(BUILD_CVODES)
  ADD_SUBDIRECTORY(cvodes)
ENDIF(BUILD_CVODES)

# IDA library

IF(BUILD_IDA)
  ADD_SUBDIRECTORY(ida)
  IF(F77_INTERFACE_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(ida/fcmix)
  ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)
ENDIF(BUILD_IDA)

# IDAS library

IF(BUILD_IDAS)
  ADD_SUBDIRECTORY(idas)
ENDIF(BUILD_IDAS)

# KINSOL library

IF(BUILD_KINSOL)
  ADD_SUBDIRECTORY(kinsol)
  IF(F77_INTERFACE_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(kinsol/fcmix)
  ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)
ENDIF(BUILD_KINSOL)

# CPODES library

IF(BUILD_CPODES)
  ADD_SUBDIRECTORY(cpodes)
ENDIF(BUILD_CPODES)
