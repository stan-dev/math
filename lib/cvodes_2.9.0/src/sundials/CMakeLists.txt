# ---------------------------------------------------------------
# $Revision: 4075 $
# $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
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
  sundials_iterative.h
  sundials_math.h
  sundials_nvector.h
  sundials_fnvector.h
  sundials_pcg.h
  sundials_sparse.h
  sundials_spbcgs.h
  sundials_spfgmr.h
  sundials_spgmr.h
  sundials_sptfqmr.h
  sundials_types.h
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
