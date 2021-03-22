# ---------------------------------------------------------------
# Programmer(s): Radu Serban and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the generic SUNDIALS modules
# ---------------------------------------------------------------

# From here we only install the generic SUNDIALS headers.
# The implementations themselves are incorporated in the individual SUNDIALS solver libraries.

install(CODE "MESSAGE(\"\nInstall shared components\n\")")

# Add variable sundials_HEADERS with the exported SUNDIALS header files
set(sundials_HEADERS
  sundials_band.h
  sundials_dense.h
  sundials_direct.h
  sundials_fnvector.h
  sundials_futils.h
  sundials_iterative.h
  sundials_linearsolver.h
  sundials_math.h
  sundials_matrix.h
  sundials_memory.h
  sundials_nonlinearsolver.h
  sundials_mpi_types.h
  sundials_nvector.h
  sundials_types.h
  sundials_version.h
  )

if(ENABLE_CUDA)
  list(APPEND sundials_HEADERS sundials_cuda_policies.hpp)
endif()

if(ENABLE_HIP)
  list(APPEND sundials_HEADERS sundials_hip_policies.hpp)
endif()

if(ENABLE_SYCL)
  list(APPEND sundials_HEADERS sundials_sycl_policies.hpp)
endif()

# If Blas/Lapack support was enabled, add the Lapack interface headers
if(BUILD_SUNLINSOL_LAPACKBAND OR BUILD_SUNLINSOL_LAPACKDENSE)
  list(APPEND sundials_HEADERS sundials_lapack.h)
endif()

# If enabled, add the XBraid interface header
if(ENABLE_XBRAID)
  list(APPEND sundials_HEADERS sundials_xbraid.h)
endif()

# Add prefix with complete path to the header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/sundials/ sundials_HEADERS)

set(sundials_SOURCES
  sundials_band.c
  sundials_dense.c
  sundials_direct.c
  sundials_futils.c
  sundials_iterative.c
  sundials_linearsolver.c
  sundials_math.c
  sundials_matrix.c
  sundials_memory.c
  sundials_nonlinearsolver.c
  sundials_nvector.c
  sundials_nvector_senswrapper.c
  sundials_version.c)

# Add prefix with complete path to the source files
add_prefix(${SUNDIALS_SOURCE_DIR}/src/sundials/ sundials_SOURCES)

# Create a library out of the generic sundials modules
sundials_add_library(sundials_generic
  SOURCES
    ${sundials_SOURCES}
  HEADERS
    ${sundials_HEADERS}
  INCLUDE_SUBDIR
    sundials
  OUTPUT_NAME
    sundials_generic
  VERSION
    ${sundialslib_VERSION}
  SOVERSION
    ${sundialslib_SOVERSION}
)

# SUNDIALS_EXPORT macro
include(GenerateExportHeader)
generate_export_header(
  sundials_generic BASE_NAME SUNDIALS
  EXPORT_FILE_NAME
  "${PROJECT_BINARY_DIR}/include/sundials/sundials_export.h"
)

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()
