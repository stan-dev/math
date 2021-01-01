# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the cuda cuSolverSp SUNLinearSolver
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_CUSOLVERSP\n\")")

# Add variable cusolversp_HEADERS with the exported cusolversp header files
set(cusolversp_HEADERS sunlinsol_cusolversp_batchqr.h)
add_prefix(${SUNDIALS_SOURCE_DIR}/include/sunlinsol/ cusolversp_HEADERS)

sundials_add_library(sundials_sunlinsolcusolversp
  SOURCES
    sunlinsol_cusolversp_batchqr.cu
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC sundials_sunmatrixcusparse ${CUDA_cusolver_LIBRARY}
    PRIVATE ${CUDA_cusparse_LIBRARY}
  OUTPUT_NAME
    sundials_sunlinsolcusolversp
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_VERSION}
)

# Install the CUDA NVector header files
install(FILES ${cusolversp_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_CUSOLVERSP module")
