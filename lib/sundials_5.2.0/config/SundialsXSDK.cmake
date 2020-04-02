# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
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
option(USE_XSDK_DEFAULTS "Enable default xSDK settings" OFF)

# ---------------------------------------------------------------
# Set default values for some xSDK variables
# ---------------------------------------------------------------

if(USE_XSDK_DEFAULTS)
  message(STATUS "Enabling xSDK defaults")

  # set the CMake build type, SUNDIALS does not set a build type by default
  if(NOT CMAKE_BUILD_TYPE)
    message("Setting build type to Debug")
    set(DOCSTR "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel")
    force_variable(CMAKE_BUILD_TYPE STRING "${DOCSTR}" "Debug")
  endif()

  # set build precision, SUNDIALS_PRECISION defaults to double
  show_variable(XSDK_PRECISION STRING "single, double, or quad" "double")

  # set build index size, SUNDIALS_INDEX_SIZE defaults to int64_t
  show_variable(XSDK_INDEX_SIZE STRING "32 or 64" "32")

  # disable Fortran-C interface, defaults to OFF
  show_variable(XSDK_ENABLE_FORTRAN BOOL "Enable Fortran-C support" OFF)

  # disable CUDA by default
  show_variable(XSDK_ENABLE_CUDA BOOL "Enable CUDA support" OFF)

  # disable LAPACK by default
  show_variable(TPL_ENABLE_LAPACK BOOL "Enable LAPACK support" OFF)

  # disable KLU by default
  show_variable(TPL_ENABLE_KLU BOOL "Enable KLU support" OFF)

  # disable hypre by default
  show_variable(TPL_ENABLE_HYPRE BOOL "Enable hypre support" OFF)

  # disable Trilinos by default
  show_variable(TPL_ENABLE_TRILINOS BOOL "Enable Trilinos support" OFF)

endif()

# ---------------------------------------------------------------
# hide (make advanced) and overwrite equivalent SUNDIALS variables
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# SUNDIALS precision, indextype, and Fortran interface
# ---------------------------------------------------------------

# XSDK_PRECISION => SUNDIALS_PRECISION
if(XSDK_PRECISION)
  message("Replacing SUNDIALS_PRECISION with XSDK_PRECISION")
  set(DOCSTR "single, double, or extended")

  if(XSDK_PRECISION MATCHES "quad")
    force_variable(SUNDIALS_PRECISION STRING "${DOCSTR}" "extended")
  else()
    force_variable(SUNDIALS_PRECISION STRING "${DOCSTR}" "${XSDK_PRECISION}")
  endif()

  mark_as_advanced(FORCE SUNDIALS_PRECISION)
endif()

# XSDK_INDEX_SIZE => SUNDIALS_INDEX_SIZE
if(XSDK_INDEX_SIZE)
  message("Replacing SUNDIALS_INDEX_SIZE with XSDK_INDEX_SIZE")
  set(DOCSTR "Signed 64-bit (64) or signed 32-bit (32) integer")
  force_variable(SUNDIALS_INDEX_SIZE STRING "${DOCSTR}" ${XSDK_INDEX_SIZE})
  mark_as_advanced(FORCE SUNDIALS_INDEX_SIZE)
endif()

# XSDK_FORTRAN_ENABLE => F77_INTERFACE_ENABLE/F2003_INTERFACE_ENABLE
if(DEFINED XSDK_ENABLE_FORTRAN)
  message("Replacing F77_INTERFACE_ENABLE and F2003_INTERFACE_ENABLE with XSDK_ENABLE_FORTRAN")
  set(DOCSTR "Enable Fortran-C support")

  # check that at least one solver with a Fortran interface is built
  if(NOT BUILD_ARKODE AND NOT BUILD_CVODE AND NOT BUILD_IDA AND NOT BUILD_KINSOL)
    if(XSDK_ENABLE_FORTRAN)
      print_warning("Enabled packages do not support Fortran"
                    "Disabeling XSDK_ENABLE_FORTRAN")
      force_variable(XSDK_ENABLE_FORTRAN BOOL "${DOCSTR}" OFF)
    endif()
    hide_variable(F77_INTERFACE_ENABLE)
    hide_variable(F2003_INTERFACE_ENABLE)
    hide_variable(XSDK_ENABLE_FORTRAN)
  endif()

  force_variable(F77_INTERFACE_ENABLE BOOL "${DOCSTR}" "${XSDK_ENABLE_FORTRAN}")
  mark_as_advanced(FORCE F77_INTERFACE_ENABLE)

  force_variable(F2003_INTERFACE_ENABLE BOOL "${DOCSTR}" "${XSDK_ENABLE_FORTRAN}")
  mark_as_advanced(FORCE F2003_INTERFACE_ENABLE)
endif()

# XSDK_ENABLE_CUDA => CUDA_ENABLE
if(DEFINED XSDK_ENABLE_CUDA)
  message("Replacing CUDA_ENABLE with XSDK_ENABLE_CUDA")
  set(DOCSTR "Enable CUDA support")

  force_variable(CUDA_ENABLE BOOL "${DOCSTR}" "${XSDK_ENABLE_CUDA}")
  mark_as_advanced(FORCE CUDA_ENABLE)
endif()

# ---------------------------------------------------------------
# LAPACK
# ---------------------------------------------------------------

# TPL_ENABLE_LAPACK => LAPACK_ENABLE
if(DEFINED TPL_ENABLE_LAPACK)
  message("Replacing LAPACK_ENABLE with TPL_ENABLE_LAPACK")
  set(DOCSTR "Enable Lapack support")

  force_variable(LAPACK_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_LAPACK}")
  mark_as_advanced(FORCE LAPACK_ENABLE)
endif()

# TPL_LAPACK_LIBRARIES => LAPACK_LIBRARIES
if(TPL_ENABLE_LAPACK)
  message("Replacing LAPACK_LIBRARIES with TPL_LAPACK_LIBRARIES")
  set(DOCSTR "Lapack library")

  show_variable(TPL_LAPACK_LIBRARIES STRING "${DOCSTR}" "${TPL_LAPACK_LIBRARIES}")
  force_variable(LAPACK_LIBRARIES STRING "${DOCSTR}" "${TPL_LAPACK_LIBRARIES}")
  mark_as_advanced(FORCE LAPACK_LIBRARIES)
endif()

# ---------------------------------------------------------------
# KLU
# ---------------------------------------------------------------

# TPL_ENABLE_KLU => KLU_ENABLE
if(DEFINED TPL_ENABLE_KLU)
  message("Replacing KLU_ENABLE with TPL_ENABLE_KLU")
  set(DOCSTR "Enable KLU support")

  force_variable(KLU_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_KLU}")
  mark_as_advanced(FORCE KLU_ENABLE)
endif()

# TPL_KLU_INCLUDE_DIRS => KLU_INCLUDE_DIR
# TPL_KLU_LIBRARIES    => KLU_LIBRARY     => KLU_LIBRARIES
if(TPL_ENABLE_KLU)
  message("Replacing KLU_INCLUDE_DIR with TPL_KLU_INCLUDE_DIRS")
  set(DOCSTR "KLU include directory")

  show_variable(TPL_KLU_INCLUDE_DIRS STRING "${DOCSTR}" "${TPL_KLU_INCLUDE_DIRS}")
  force_variable(KLU_INCLUDE_DIR STRING "${DOCSTR}" "${TPL_KLU_INCLUDE_DIRS}")
  mark_as_advanced(FORCE KLU_INCLUDE_DIR)

  message("Replacing KLU_LIBRARY with TPL_KLU_LIBRARIES")
  set(DOCSTR "KLU library")

  show_variable(TPL_KLU_LIBRARIES STRING "${DOCSTR}" "${TPL_KLU_LIBRARIES}")
  force_variable(KLU_LIBRARY STRING "${DOCSTR}" "${TPL_KLU_LIBRARIES}")
  mark_as_advanced(FORCE KLU_LIBRARY)
  mark_as_advanced(FORCE KLU_LIBRARY_DIR)
endif()


# ---------------------------------------------------------------
# HYPRE
# ---------------------------------------------------------------

# TPL_ENABLE_HYPRE => HYPRE_ENABLE
if(DEFINED TPL_ENABLE_HYPRE)
  message("Replacing HYPRE_ENABLE with TPL_ENABLE_HYPRE")
  set(DOCSTR "Enable hypre support")

  force_variable(HYPRE_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_HYPRE}")
  mark_as_advanced(FORCE HYPRE_ENABLE)
endif()

# TPL_HYPRE_INCLUDE_DIRS => HYPRE_INCLUDE_DIR
# TPL_HYPRE_LIBRARIES    => HYPRE_LIBRARY     => HYPRE_LIBRARIES
if(TPL_ENABLE_HYPRE)
  message("Replacing HYPRE_INCLUDE_DIR with TPL_HYPRE_INCLUDE_DIRS")
  set(DOCSTR "hypre include directory")

  show_variable(TPL_HYPRE_INCLUDE_DIRS STRING "${DOCSTR}" "${TPL_HYPRE_INCLUDE_DIRS}")
  force_variable(HYPRE_INCLUDE_DIR STRING "${DOCSTR}" "${TPL_HYPRE_INCLUDE_DIRS}")
  mark_as_advanced(FORCE HYPRE_INCLUDE_DIR)

  message("Replacing HYPRE_LIBRARY with TPL_HYPRE_LIBRARIES")
  set(DOCSTR "hypre library")

  show_variable(TPL_HYPRE_LIBRARIES STRING "${DOCSTR}" "${TPL_HYPRE_LIBRARIES}")
  force_variable(HYPRE_LIBRARY STRING "${DOCSTR}" "${TPL_HYPRE_LIBRARIES}")
  mark_as_advanced(FORCE HYPRE_LIBRARY)
  mark_as_advanced(FORCE HYPRE_LIBRARY_DIR)
endif()

# ---------------------------------------------------------------
# Trilinos
# ---------------------------------------------------------------

# TPL_ENABLE_TRILINOS => Trilinos_ENABLE
if(DEFINED TPL_ENABLE_TRILINOS)
  message("Replacing Trilinos_ENABLE with TPL_ENABLE_TRILINOS")
  set(DOCSTR "Enable Trilinos support")

  force_variable(Trilinos_ENABLE BOOL "${DOCSTR}" "${TPL_ENABLE_TRILINOS}")
  mark_as_advanced(FORCE Trilinos_ENABLE)
endif()
