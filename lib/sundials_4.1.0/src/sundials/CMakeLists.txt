# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
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
# CMakeLists.txt file for the generic SUNDIALS modules

# From here we only install the generic SUNDIALS headers.
# The implementations themselves are incorporated in the individual SUNDIALS solver libraries.

INSTALL(CODE "MESSAGE(\"\nInstall shared components\n\")")

# Add variable sundials_HEADERS with the exported SUNDIALS header files
SET(sundials_HEADERS
  sundials_band.h
  sundials_dense.h
  sundials_direct.h
  sundials_fnvector.h
  sundials_iterative.h
  sundials_linearsolver.h
  sundials_math.h
  sundials_matrix.h
  sundials_mpi.h
  sundials_mpi_types.h
  sundials_nonlinearsolver.h
  sundials_mpi.h
  sundials_mpi_types.h
  sundials_nvector.h
  sundials_pcg.h
  sundials_sparse.h
  sundials_spbcgs.h
  sundials_spfgmr.h
  sundials_spgmr.h
  sundials_sptfqmr.h
  sundials_types.h
  sundials_version.h
  )

# Add prefix with complete path to the SUNDIALS header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/sundials/ sundials_HEADERS)

# Install the SUNDIALS header files
INSTALL(FILES ${sundials_HEADERS} DESTINATION include/sundials)

# If Blas/Lapack support was enabled, install the Lapack interface headers
IF(LAPACK_FOUND)
  SET(sundials_BL_HEADERS sundials_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/sundials/ sundials_BL_HEADERS)
  INSTALL(FILES ${sundials_BL_HEADERS} DESTINATION include/sundials)
ENDIF(LAPACK_FOUND)
