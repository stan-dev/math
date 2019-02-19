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
# CMakeLists.txt file for the CVODE library

INSTALL(CODE "MESSAGE(\"\nInstall CVODE\n\")")

# Add F90 module if F2003 interface is enabled
IF(F90_FOUND AND F2003_INTERFACE_ENABLE)
  ADD_SUBDIRECTORY(F90)
ENDIF(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable cvode_SOURCES with the sources for the CVODE library
SET(cvode_SOURCES
  cvode.c
  cvode_bandpre.c
  cvode_bbdpre.c
  cvode_diag.c
  cvode_direct.c
  cvode_io.c
  cvode_ls.c
  cvode_nls.c
  cvode_spils.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CVODE library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nonlinearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_band.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_dense.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_direct.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_version.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector_senswrapper.c
  ${sundials_SOURCE_DIR}/src/nvector/serial/nvector_serial.c
  )

# Add variable sunmatrix_SOURCES with the common SUNMatrix sources which will
# also be included in the CVODE library
SET(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/sunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/sunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/sunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the CVODE library
SET(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol/band/sunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/dense/sunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spbcgs/sunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spfgmr/sunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spgmr/sunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/pcg/sunlinsol_pcg.c
  )

# Add variable sunnonlinsol_SOURCES with the common SUNNonlinearSolver sources
# which will also be included in the CVODE library
SET(sunnonlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/newton/sunnonlinsol_newton.c
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c
  )

# Add variable cvode_HEADERS with the exported CVODE header files
SET(cvode_HEADERS
  cvode.h
  cvode_bandpre.h
  cvode_bbdpre.h
  cvode_diag.h
  cvode_direct.h
  cvode_ls.h
  cvode_spils.h
  )

# Add prefix with complete path to the CVODE header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/cvode/ cvode_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)
  # Add the build target for the static CVODE library
  ADD_LIBRARY(sundials_cvode_static STATIC
    ${cvode_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES} ${sunnonlinsol_SOURCES}
  )

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cvode_static
    PROPERTIES OUTPUT_NAME sundials_cvode CLEAN_DIRECT_OUTPUT 1)

  # Install the CVODE library
  INSTALL(TARGETS sundials_cvode_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)
  # Add the build target for the shared CVODE library
  ADD_LIBRARY(sundials_cvode_shared SHARED
    ${cvode_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES} ${sunnonlinsol_SOURCES}
  )

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_cvode_shared m)
  ENDIF()

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cvode_shared
    PROPERTIES OUTPUT_NAME sundials_cvode CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_cvode_shared
    PROPERTIES VERSION ${cvodelib_VERSION} SOVERSION ${cvodelib_SOVERSION})

  # Install the CVODE library
  INSTALL(TARGETS sundials_cvode_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

ENDIF(BUILD_SHARED_LIBS)

# Install the CVODE header files
INSTALL(FILES ${cvode_HEADERS} DESTINATION include/cvode)

#
MESSAGE(STATUS "Added CVODE module")
