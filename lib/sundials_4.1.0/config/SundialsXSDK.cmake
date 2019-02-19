# ---------------------------------------------------------------
# Programmer:  David J. Gardner @ LLNL
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
# xSDK specific CMake variables. If set, these variables will 
# overwrite the value in equivalent SUNDIALS CMake variable.
#
# Only USE_XSDK_DEFAULTS is created in CACHE by default (set to 
# OFF). The other xSDK variables are left undefined. They can be
# be set by passing -D<variable_name>=<value> to cmake or can be 
# enabled in the cmake-gui setting USE_XSDK_DEFAULTS to ON.
#
# When USE_XSDK_DEFAULTS is ON the default values are overwritten
# by values passed to cmake or manually set in the cmake-gui.
# ---------------------------------------------------------------

# always show the option to turn on xSDK defaults
OPTION(USE_XSDK_DEFAULTS "Enable default xSDK settings" OFF)

# ---------------------------------------------------------------
# Set default values for some xSDK variables
# ---------------------------------------------------------------

IF(USE_XSDK_DEFAULTS)

  MESSAGE("Enabeling xSDK defaults")
  
  # set the CMake build type, SUNDIALS does not set a build type by default
  IF(NOT CMAKE_BUILD_TYPE)
    MESSAGE("Setting build type to Debug")
    SET(DOCSTR "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel")
    FORCE_VARIABLE(CMAKE_BUILD_TYPE STRING "${DOCSTR}" "Debug")
  ENDIF()

  # set build precision, SUNDIALS_PRECISION defaults to double
  SHOW_VARIABLE(XSDK_PRECISION STRING "single, double, or quad" "double")

  # set build index size, SUNDIALS_INDEX_SIZE defaults to int64_t
  SHOW_VARIABLE(XSDK_INDEX_SIZE STRING "32 or 64" "32")

  # disable Fortran-C interface, defaults to OFF
  SHOW_VARIABLE(XSDK_ENABLE_FORTRAN BOOL "Enable Fortran-C support" OFF)

  # disable CUDA by default
  SHOW_VARIABLE(XSDK_ENABLE_CUDA BOOL "Enable CUDA support" OFF)

  # disable BLAS by default
  SHOW_VARIABLE(TPL_ENABLE_BLAS BOOL "Enable BLAS support" OFF)

  # disable LAPACK by default
  SHOW_VARIABLE(TPL_ENABLE_LAPACK BOOL "Enable LAPACK support" OFF)

  # disable SuperLU_MT by default
  SHOW_VARIABLE(TPL_ENABLE_SUPERLUMT BOOL "Enable SuperLU_MT support" OFF)

  # disable KLU by default
  SHOW_VARIABLE(TPL_ENABLE_KLU BOOL "Enable KLU support" OFF)

  # disable PETSc by default
  SHOW_VARIABLE(TPL_ENABLE_PETSC BOOL "Enable PETSc support" OFF)

  # disable hypre by default
  SHOW_VARIABLE(TPL_ENABLE_HYPRE BOOL "Enable hypre support" OFF)

  # disable Trilinos by default
  SHOW_VARIABLE(TPL_ENABLE_TRILINOS BOOL "Enable Trilinos support" OFF)

  # disable RAJA by default
  # SHOW_VARIABLE(TPL_ENABLE_RAJA BOOL "Enable RAJA support" OFF)

ENDIF()

# ---------------------------------------------------------------
# hide (make advanced) and overwrite equivalent SUNDIALS variables
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# SUNDIALS precision, indextype, and Fortran interface
# ---------------------------------------------------------------

# XSDK_PRECISION => SUNDIALS_PRECISION
IF(XSDK_PRECISION)
  MESSAGE("Replacing SUNDIALS_PRECISION with XSDK_PRECISION")
  SET(DOCSTR "single, double, or extended")

  IF(XSDK_PRECISION MATCHES "quad")
    FORCE_VARIABLE(SUNDIALS_PRECISION STRING "${DOCSTR}" "extended")
  ELSE()
    FORCE_VARIABLE(SUNDIALS_PRECISION STRING "${DOCSTR}" "${XSDK_PRECISION}")
  ENDIF()

  MARK_AS_ADVANCED(FORCE SUNDIALS_PRECISION)
ENDIF()

# XSDK_INDEX_SIZE => SUNDIALS_INDEX_SIZE
IF(XSDK_INDEX_SIZE)
  MESSAGE("Replacing SUNDIALS_INDEX_SIZE with XSDK_INDEX_SIZE")
  SET(DOCSTR "Signed 64-bit (64) or signed 32-bit (32) integer")
  FORCE_VARIABLE(SUNDIALS_INDEX_SIZE STRING "${DOCSTR}" ${XSDK_INDEX_SIZE})
  MARK_AS_ADVANCED(FORCE SUNDIALS_INDEX_SIZE)
ENDIF()

# XSDK_FORTRAN_ENABLE => F77_INTERFACE_ENABLE/F2003_INTERFACE_ENABLE
IF(DEFINED XSDK_ENABLE_FORTRAN)
  MESSAGE("Replacing F77_INTERFACE_ENABLE and F2003_INTERFACE_ENABLE with XSDK_ENABLE_FORTRAN")
  SET(DOCSTR "Enable Fortran-C support")
  
  # check that at least one solver with a Fortran interface is built
  IF(NOT BUILD_ARKODE AND NOT BUILD_CVODE AND NOT BUILD_IDA AND NOT BUILD_KINSOL)
    IF(XSDK_ENABLE_FORTRAN)
      PRINT_WARNING("Enabled packages do not support Fortran" 
                    "Disabeling XSDK_ENABLE_FORTRAN")
      FORCE_VARIABLE(XSDK_ENABLE_FORTRAN BOOL "${DOCSTR}" OFF)
    ENDIF()
    HIDE_VARIABLE(F77_INTERFACE_ENABLE)
    HIDE_VARIABLE(F2003_INTERFACE_ENABLE)
    HIDE_VARIABLE(XSDK_ENABLE_FORTRAN)
  ENDIF()

  FORCE_VARIABLE(F77_INTERFACE_ENABLE BOOL "${DOCSTR}" "${XSDK_ENABLE_FORTRAN}")
  MARK_AS_ADVANCED(FORCE F77_INTERFACE_ENABLE)
  
  FORCE_VARIABLE(F2003_INTERFACE_ENABLE BOOL "${DOCSTR}" "${XSDK_ENABLE_FORTRAN}")
  MARK_AS_ADVANCED(FORCE F2003_INTERFACE_ENABLE)
ENDIF()

# XSDK_ENABLE_CUDA => CUDA_ENABLE
IF(DEFINED XSDK_ENABLE_CUDA)
  MESSAGE("Replacing CUDA_ENABLE with XSDK_ENABLE_CUDA")
  SET(DOCSTR "Enable CUDA support")

  FORCE_VARIABLE(CUDA_ENABLE BOOL "${DOCSTR}" "${XSDK_ENABLE_CUDA}")
  MARK_AS_ADVANCED(FORCE CUDA_ENABLE)
ENDIF()

# ---------------------------------------------------------------
# BLAS
# ---------------------------------------------------------------

# TPL_ENABLE_BLAS => BLAS_ENABLE
IF(DEFINED TPL_ENABLE_BLAS)
  MESSAGE("Replacing BLAS_ENABLE with TPL_ENABLE_BLAS")
  SET(DOCSTR "Enable Blas support")

  FORCE_VARIABLE(BLAS_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_BLAS}")
  MARK_AS_ADVANCED(FORCE BLAS_ENABLE)
ENDIF()

# TPL_BLAS_LIBRARIES => BLAS_LIBRARIES
IF(TPL_ENABLE_BLAS)
  MESSAGE("Replacing BLAS_LIBRARIES with TPL_BLAS_LIBRARIES")
  SET(DOCSTR "Blas library")

  SHOW_VARIABLE(TPL_BLAS_LIBRARIES STRING "${DOCSTR}" "${TPL_BLAS_LIBRARIES}")
  FORCE_VARIABLE(BLAS_LIBRARIES STRING "${DOCSTR}" "${TPL_BLAS_LIBRARIES}")
  MARK_AS_ADVANCED(FORCE BLAS_LIBRARIES)
ENDIF()


# ---------------------------------------------------------------
# LAPACK
# ---------------------------------------------------------------

# TPL_ENABLE_LAPACK => LAPACK_ENABLE
IF(DEFINED TPL_ENABLE_LAPACK)
  MESSAGE("Replacing LAPACK_ENABLE with TPL_ENABLE_LAPACK")
  SET(DOCSTR "Enable Lapack support")

  FORCE_VARIABLE(LAPACK_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_LAPACK}")
  MARK_AS_ADVANCED(FORCE LAPACK_ENABLE)
ENDIF()

# TPL_LAPACK_LIBRARIES => LAPACK_LIBRARIES
IF(TPL_ENABLE_LAPACK)
  MESSAGE("Replacing LAPACK_LIBRARIES with TPL_LAPACK_LIBRARIES")
  SET(DOCSTR "Lapack library")

  SHOW_VARIABLE(TPL_LAPACK_LIBRARIES STRING "${DOCSTR}" "${TPL_LAPACK_LIBRARIES}")
  FORCE_VARIABLE(LAPACK_LIBRARIES STRING "${DOCSTR}" "${TPL_LAPACK_LIBRARIES}")
  MARK_AS_ADVANCED(FORCE LAPACK_LIBRARIES)
ENDIF()


# ---------------------------------------------------------------
# SuperLU_MT
# ---------------------------------------------------------------

# TPL_ENABLE_SUPERLUMT => SUPERLUMT_ENABLE
IF(DEFINED TPL_ENABLE_SUPERLUMT)
  MESSAGE("Replacing SUPERLUMT_ENABLE with TPL_ENABLE_SUPERLUMT")
  SET(DOCSTR "Enable SuperLU_MT support")

  FORCE_VARIABLE(SUPERLUMT_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_SUPERLUMT}")
  MARK_AS_ADVANCED(FORCE SUPERLUMT_ENABLE)
ENDIF()

# TPL_SUPERLUMT_INCLUDE_DIRS => SUPERLUMT_INCLUDE_DIR
# TPL_SUPERLUMT_LIBRARIES    => SUPERLUMT_LIBRARY     => SUPERLUMT_LIBRARIES
IF(TPL_ENABLE_SUPERLUMT)
  MESSAGE("Replacing SUPERLUMT_INCLUDE_DIR with TPL_SUPERLUMT_INCLUDE_DIRS")
  SET(DOCSTR "SuperLU_MT include directory")

  SHOW_VARIABLE(TPL_SUPERLUMT_INCLUDE_DIRS STRING "${DOCSTR}" "${TPL_SUPERLUMT_INCLUDE_DIRS}")
  FORCE_VARIABLE(SUPERLUMT_INCLUDE_DIR STRING "${DOCSTR}" "${TPL_SUPERLUMT_INCLUDE_DIRS}")
  MARK_AS_ADVANCED(FORCE SUPERLUMT_INCLUDE_DIR)

  MESSAGE("Replacing SUPERLUMT_LIBRARY with TPL_SUPERLUMT_LIBRARIES")
  SET(DOCSTR "SuperLU_MT library")

  SHOW_VARIABLE(TPL_SUPERLUMT_LIBRARIES STRING "${DOCSTR}" "${TPL_SUPERLUMT_LIBRARIES}")
  FORCE_VARIABLE(SUPERLUMT_LIBRARY STRING "${DOCSTR}" "${TPL_SUPERLUMT_LIBRARIES}")
  MARK_AS_ADVANCED(FORCE SUPERLUMT_LIBRARY)
  MARK_AS_ADVANCED(FORCE SUPERLUMT_LIBRARY_DIR)

  MESSAGE("Replacing SUPERLUMT_THREAD_TYPE with TPL_SUPERLUMT_THREAD_TYPE")
  SET(DOCSTR "SuperLU_MT thread type (OpenMP or Pthread)")
  
  SHOW_VARIABLE(TPL_SUPERLUMT_THREAD_TYPE STRING "${DOCSTR}" "PThread")
  FORCE_VARIABLE(SUPERLUMT_THREAD_TYPE STRING "${DOCSTR}" "${TPL_SUPERLUMT_THREAD_TYPE}")
  MARK_AS_ADVANCED(FORCE SUPERLUMT_THREAD_TYPE)
ENDIF()


# ---------------------------------------------------------------
# KLU
# ---------------------------------------------------------------

# TPL_ENABLE_KLU => KLU_ENABLE
IF(DEFINED TPL_ENABLE_KLU)
  MESSAGE("Replacing KLU_ENABLE with TPL_ENABLE_KLU")
  SET(DOCSTR "Enable KLU support")

  FORCE_VARIABLE(KLU_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_KLU}")
  MARK_AS_ADVANCED(FORCE KLU_ENABLE)
ENDIF()

# TPL_KLU_INCLUDE_DIRS => KLU_INCLUDE_DIR
# TPL_KLU_LIBRARIES    => KLU_LIBRARY     => KLU_LIBRARIES
IF(TPL_ENABLE_KLU)
  MESSAGE("Replacing KLU_INCLUDE_DIR with TPL_KLU_INCLUDE_DIRS")
  SET(DOCSTR "KLU include directory")

  SHOW_VARIABLE(TPL_KLU_INCLUDE_DIRS STRING "${DOCSTR}" "${TPL_KLU_INCLUDE_DIRS}")
  FORCE_VARIABLE(KLU_INCLUDE_DIR STRING "${DOCSTR}" "${TPL_KLU_INCLUDE_DIRS}")
  MARK_AS_ADVANCED(FORCE KLU_INCLUDE_DIR)

  MESSAGE("Replacing KLU_LIBRARY with TPL_KLU_LIBRARIES")
  SET(DOCSTR "KLU library")

  SHOW_VARIABLE(TPL_KLU_LIBRARIES STRING "${DOCSTR}" "${TPL_KLU_LIBRARIES}")
  FORCE_VARIABLE(KLU_LIBRARY STRING "${DOCSTR}" "${TPL_KLU_LIBRARIES}")
  MARK_AS_ADVANCED(FORCE KLU_LIBRARY)
  MARK_AS_ADVANCED(FORCE KLU_LIBRARY_DIR)
ENDIF()


# ---------------------------------------------------------------
# HYPRE
# ---------------------------------------------------------------

# TPL_ENABLE_HYPRE => HYPRE_ENABLE
IF(DEFINED TPL_ENABLE_HYPRE)
  MESSAGE("Replacing HYPRE_ENABLE with TPL_ENABLE_HYPRE")
  SET(DOCSTR "Enable hypre support")

  FORCE_VARIABLE(HYPRE_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_HYPRE}")
  MARK_AS_ADVANCED(FORCE HYPRE_ENABLE)
ENDIF()

# TPL_HYPRE_INCLUDE_DIRS => HYPRE_INCLUDE_DIR
# TPL_HYPRE_LIBRARIES    => HYPRE_LIBRARY     => HYPRE_LIBRARIES
IF(TPL_ENABLE_HYPRE)
  MESSAGE("Replacing HYPRE_INCLUDE_DIR with TPL_HYPRE_INCLUDE_DIRS")
  SET(DOCSTR "hypre include directory")

  SHOW_VARIABLE(TPL_HYPRE_INCLUDE_DIRS STRING "${DOCSTR}" "${TPL_HYPRE_INCLUDE_DIRS}")
  FORCE_VARIABLE(HYPRE_INCLUDE_DIR STRING "${DOCSTR}" "${TPL_HYPRE_INCLUDE_DIRS}")
  MARK_AS_ADVANCED(FORCE HYPRE_INCLUDE_DIR)

  MESSAGE("Replacing HYPRE_LIBRARY with TPL_HYPRE_LIBRARIES")
  SET(DOCSTR "hypre library")

  SHOW_VARIABLE(TPL_HYPRE_LIBRARIES STRING "${DOCSTR}" "${TPL_HYPRE_LIBRARIES}")
  FORCE_VARIABLE(HYPRE_LIBRARY STRING "${DOCSTR}" "${TPL_HYPRE_LIBRARIES}")
  MARK_AS_ADVANCED(FORCE HYPRE_LIBRARY)
  MARK_AS_ADVANCED(FORCE HYPRE_LIBRARY_DIR)
ENDIF()


# ---------------------------------------------------------------
# PETSC
# ---------------------------------------------------------------

# TPL_ENABLE_PETSC => PETSC_ENABLE
IF(DEFINED TPL_ENABLE_PETSC)
  MESSAGE("Replacing PETSC_ENABLE with TPL_ENABLE_PETSC")
  SET(DOCSTR "Enable petsc support")

  FORCE_VARIABLE(PETSC_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_PETSC}")
  MARK_AS_ADVANCED(FORCE PETSC_ENABLE)
ENDIF()

# TPL_PETSC_INCLUDE_DIRS => PETSC_INCLUDE_DIR
# TPL_PETSC_LIBRARIES    => PETSC_LIBRARY     => PETSC_LIBRARIES
IF(TPL_ENABLE_PETSC)
  MESSAGE("Replacing PETSC_INCLUDE_DIR with TPL_PETSC_INCLUDE_DIRS")
  SET(DOCSTR "PETSc include directory")

  SHOW_VARIABLE(TPL_PETSC_INCLUDE_DIRS STRING "${DOCSTR}" "${TPL_PETSC_INCLUDE_DIRS}")
  FORCE_VARIABLE(PETSC_INCLUDE_DIR STRING "${DOCSTR}" "${TPL_PETSC_INCLUDE_DIRS}")
  MARK_AS_ADVANCED(FORCE PETSC_INCLUDE_DIR)

  MESSAGE("Replacing PETSC_LIBRARY with TPL_PETSC_LIBRARIES")
  SET(DOCSTR "PETSc library")

  SHOW_VARIABLE(TPL_PETSC_LIBRARIES STRING "${DOCSTR}" "${TPL_PETSC_LIBRARIES}")
  FORCE_VARIABLE(PETSC_LIBRARY STRING "${DOCSTR}" "${TPL_PETSC_LIBRARIES}")
  MARK_AS_ADVANCED(FORCE PETSC_LIBRARY)
  MARK_AS_ADVANCED(FORCE PETSC_LIBRARY_DIR)
ENDIF()

# ---------------------------------------------------------------
# Trilinos
# ---------------------------------------------------------------

# TPL_ENABLE_TRILINOS => Trilinos_ENABLE
IF(DEFINED TPL_ENABLE_TRILINOS)
  MESSAGE("Replacing Trilinos_ENABLE with TPL_ENABLE_TRILINOS")
  SET(DOCSTR "Enable Trilinos support")

  FORCE_VARIABLE(Trilinos_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_TRILINOS}")
  MARK_AS_ADVANCED(FORCE Trilinos_ENABLE)
ENDIF()

# RAJA
# ---------------------------------------------------------------

# # TPL_ENABLE_RAJA => RAJA_ENABLE
# IF(DEFINED TPL_ENABLE_RAJA)
#   MESSAGE("Replacing RAJA_ENABLE with TPL_ENABLE_RAJA")
#   SET(DOCSTR "Enable RAJA support")

#   FORCE_VARIABLE(RAJA_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_RAJA}")
#   MARK_AS_ADVANCED(FORCE RAJA_ENABLE)
# ENDIF()

# # TPL_RAJA_DIR => RAJA_DIR
# IF(TPL_ENABLE_RAJA)
#   MESSAGE("Replacing RAJA_DIR with TPL_RAJA_DIR")
#   SET(DOCSTR "RAJA include directory")

#   SHOW_VARIABLE(TPL_RAJA_DIR STRING "${DOCSTR}" "${TPL_RAJA_DIR}")
#   FORCE_VARIABLE(RAJA_DIR STRING "${DOCSTR}" "${TPL_RAJA_DIR}")
#   MARK_AS_ADVANCED(FORCE RAJA_DIR)
# ENDIF()
