# ---------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#                 Radu Serban @ LLNL
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
# CMakeLists.txt file for the KINSOL library

INSTALL(CODE "MESSAGE(\"\nInstall KINSOL\n\")")

# Add variable kinsol_SOURCES with the sources for the KINSOL library
SET(kinsol_SOURCES
  kinsol.c
  kinsol_bbdpre.c
  kinsol_direct.c
  kinsol_io.c
  kinsol_ls.c
  kinsol_spils.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the KINSOL library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_band.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_dense.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_direct.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_version.c
  ${sundials_SOURCE_DIR}/src/nvector/serial/nvector_serial.c
  )

# Add variable sunmatrix_SOURCES with the common SUNMatrix sources which will
# also be included in the KINSOL library
SET(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/sunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/sunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/sunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the KINSOL library
SET(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol/band/sunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/dense/sunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spbcgs/sunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spfgmr/sunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spgmr/sunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/pcg/sunlinsol_pcg.c
  )

# Add variable kinsol_HEADERS with the exported KINSOL header files
SET(kinsol_HEADERS
  kinsol.h
  kinsol_bbdpre.h
  kinsol_direct.h
  kinsol_ls.h
  kinsol_spils.h
  )

# Add prefix with complete path to the KINSOL header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/kinsol/ kinsol_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static KINSOL library
  ADD_LIBRARY(sundials_kinsol_static STATIC
    ${kinsol_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_kinsol_static
    PROPERTIES OUTPUT_NAME sundials_kinsol CLEAN_DIRECT_OUTPUT 1)

  # Install the KINSOL library
  INSTALL(TARGETS sundials_kinsol_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the KINSOL library
  ADD_LIBRARY(sundials_kinsol_shared SHARED
    ${kinsol_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_kinsol_shared m)
  ENDIF()

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_kinsol_shared
    PROPERTIES OUTPUT_NAME sundials_kinsol CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_kinsol_shared
    PROPERTIES VERSION ${kinsollib_VERSION} SOVERSION ${kinsollib_SOVERSION})

  # Install the KINSOL library
  INSTALL(TARGETS sundials_kinsol_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

ENDIF(BUILD_SHARED_LIBS)

# Install the KINSOL header files
INSTALL(FILES ${kinsol_HEADERS} DESTINATION include/kinsol)

#
MESSAGE(STATUS "Added KINSOL module")
