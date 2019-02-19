# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the ARKODE library

INSTALL(CODE "MESSAGE(\"\nInstall ARKODE\n\")")

# Add variable arkode_SOURCES with the sources for the ARKODE library
SET(arkode_SOURCES
  arkode.c
  arkode_adapt.c
  arkode_arkstep.c
  arkode_arkstep_io.c
  arkode_arkstep_nls.c
  arkode_bandpre.c
  arkode_bbdpre.c
  arkode_butcher.c
  arkode_butcher_dirk.c
  arkode_butcher_erk.c
  arkode_erkstep.c
  arkode_erkstep_io.c
  arkode_interp.c
  arkode_io.c
  arkode_ls.c
  arkode_mristep.c
  arkode_mristep_io.c
  arkode_root.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the ARKODE library
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
# also be included in the ARKODE library
SET(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/sunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/sunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/sunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the ARKODE library
SET(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol/band/sunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/dense/sunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spbcgs/sunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spfgmr/sunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spgmr/sunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/pcg/sunlinsol_pcg.c
  )

# Add variable sunnonlinsol_SOURCES with the common SUNNonlinearSolver sources which will
# also be included in the ARKODE library
SET(sunnonlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/newton/sunnonlinsol_newton.c
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c
  )

# Add variable arkode_HEADERS with the exported ARKODE header files
SET(arkode_HEADERS
  arkode.h
  arkode_arkstep.h
  arkode_bandpre.h
  arkode_bbdpre.h
  arkode_butcher.h
  arkode_butcher_erk.h
  arkode_butcher_dirk.h
  arkode_erkstep.h
  arkode_ls.h
  arkode_mristep.h
  )

# Add prefix with complete path to the ARKODE header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/arkode/ arkode_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static ARKODE library
  ADD_LIBRARY(sundials_arkode_static STATIC
    ${arkode_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES} ${sunnonlinsol_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_arkode_static
    PROPERTIES OUTPUT_NAME sundials_arkode CLEAN_DIRECT_OUTPUT 1)

  # Install the ARKODE library
  INSTALL(TARGETS sundials_arkode_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the ARKODE library
  ADD_LIBRARY(sundials_arkode_shared SHARED
    ${arkode_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES} ${sunnonlinsol_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_arkode_shared m)
  ENDIF()

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_arkode_shared
    PROPERTIES OUTPUT_NAME sundials_arkode CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_arkode_shared
    PROPERTIES VERSION ${arkodelib_VERSION} SOVERSION ${arkodelib_SOVERSION})

  # Install the ARKODE library
  INSTALL(TARGETS sundials_arkode_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

ENDIF(BUILD_SHARED_LIBS)

# Install the ARKODE header files
INSTALL(FILES ${arkode_HEADERS} DESTINATION include/arkode)

#
MESSAGE(STATUS "Added ARKODE module")
