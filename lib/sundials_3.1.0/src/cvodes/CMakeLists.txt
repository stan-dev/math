# ---------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#                 Radu Serban @ LLNL
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2017, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODES library

INSTALL(CODE "MESSAGE(\"\nInstall CVODES\n\")")

# Add variable cvodes_SOURCES with the sources for the CVODES library
SET(cvodes_SOURCES
  cvodea.c
  cvodea_io.c
  cvodes.c
  cvodes_bandpre.c
  cvodes_bbdpre.c
  cvodes_io.c
  cvodes_diag.c
  cvodes_direct.c
  cvodes_spils.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CVODES library
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
  ${sundials_SOURCE_DIR}/src/nvec_ser/nvector_serial.c
  )

# Add variable sunmatrix_SOURCES with the common SUNMatrix sources which will
# also be included in the CVODES library
SET(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmat_band/sunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmat_dense/sunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmat_sparse/sunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the CVODES library
SET(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol_band/sunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_dense/sunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_spbcgs/sunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_spfgmr/sunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_spgmr/sunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_sptfqmr/sunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_pcg/sunlinsol_pcg.c
  )

# Add variable cvodes_HEADERS with the exported CVODES header files
SET(cvodes_HEADERS
  cvodes.h
  cvodes_bandpre.h
  cvodes_bbdpre.h
  cvodes_diag.h
  cvodes_direct.h
  cvodes_spils.h
  )

# Add prefix with complete path to the CVODES header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/cvodes/ cvodes_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static CVODES library
  ADD_LIBRARY(sundials_cvodes_static STATIC 
    ${cvodes_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cvodes_static
    PROPERTIES OUTPUT_NAME sundials_cvodes CLEAN_DIRECT_OUTPUT 1)

  # Install the CVODES library
  INSTALL(TARGETS sundials_cvodes_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the CVODES library
  ADD_LIBRARY(sundials_cvodes_shared SHARED 
    ${cvodes_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_cvodes_shared m)
  ENDIF()

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cvodes_shared
    PROPERTIES OUTPUT_NAME sundials_cvodes CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_cvodes_shared
    PROPERTIES VERSION ${cvodeslib_VERSION} SOVERSION ${cvodeslib_SOVERSION})

  # Install the CVODES library
  INSTALL(TARGETS sundials_cvodes_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the CVODES header files
INSTALL(FILES ${cvodes_HEADERS} DESTINATION include/cvodes)

# Install the CVODES implementation header file
INSTALL(FILES cvodes_impl.h DESTINATION include/cvodes)

#
MESSAGE(STATUS "Added CVODES module")
