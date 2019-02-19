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
# examples/arkode level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# C examples
IF(EXAMPLES_ENABLE_C)
  ADD_SUBDIRECTORY(C_serial)
  IF(OPENMP_ENABLE AND OPENMP_FOUND)
    ADD_SUBDIRECTORY(C_openmp)
  ENDIF()
  IF(OPENMP_DEVICE_ENABLE AND OPENMP_FOUND AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    ADD_SUBDIRECTORY(C_openmpdev)
  ENDIF()
  IF(MPI_ENABLE AND MPI_C_FOUND)
    ADD_SUBDIRECTORY(C_parallel)
  ENDIF()
  IF(HYPRE_ENABLE AND HYPRE_FOUND)
    ADD_SUBDIRECTORY(C_parhyp)
  ENDIF()
ENDIF()

# C++ examples
IF(EXAMPLES_ENABLE_CXX AND CXX_FOUND)
  ADD_SUBDIRECTORY(CXX_serial)
  IF(MPI_ENABLE AND MPI_CXX_FOUND)
    ADD_SUBDIRECTORY(CXX_parallel)
  ENDIF()
  IF(MPI_ENABLE AND MPI_CXX_FOUND AND HYPRE_ENABLE AND HYPRE_FOUND)
    ADD_SUBDIRECTORY(CXX_parhyp)
  ENDIF()
ENDIF()

# F77 examples
IF(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F77 AND F77_FOUND)
  ADD_SUBDIRECTORY(F77_serial)
  IF(MPI_ENABLE AND MPI_Fortran_FOUND)
    ADD_SUBDIRECTORY(F77_parallel)
  ENDIF()
ENDIF()

# F90 examples
IF(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F90 AND F90_FOUND)
  ADD_SUBDIRECTORY(F90_serial)
  IF(MPI_ENABLE AND MPI_Fortran_FOUND)
    # currently the only F90_parallel example
    # does not support anything but double precision
    IF(SUNDIALS_PRECISION MATCHES "DOUBLE")
      ADD_SUBDIRECTORY(F90_parallel)
    ENDIF()
  ENDIF()
ENDIF()
