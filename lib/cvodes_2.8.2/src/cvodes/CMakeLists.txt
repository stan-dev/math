# ---------------------------------------------------------------
# $Revision: 4135 $
# $Date: 2014-06-13 15:55:05 -0700 (Fri, 13 Jun 2014) $
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
# CMakeLists.txt file for the CVODES library

INSTALL(CODE "MESSAGE(\"\nInstall CVODES\n\")")

# Add variable cvodes_SOURCES with the sources for the CVODES library
SET(cvodes_SOURCES
  cvodes.c
  cvodea.c
  cvodes_io.c
  cvodea_io.c
  cvodes_direct.c
  cvodes_band.c
  cvodes_dense.c
  cvodes_diag.c
  cvodes_sparse.c
  cvodes_spils.c
  cvodes_spbcgs.c
  cvodes_spgmr.c
  cvodes_sptfqmr.c
  cvodes_bandpre.c
  cvodes_bbdpre.c
  )

IF(KLU_FOUND)
    LIST(APPEND cvodes_SOURCES cvodes_klu.c)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND cvodes_SOURCES cvodes_superlumt.c)
ENDIF()

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CVODES library
SET(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
  sundials_direct.c
  sundials_band.c
  sundials_dense.c
  sundials_iterative.c
  sundials_sparse.c
  sundials_spbcgs.c
  sundials_spgmr.c
  sundials_sptfqmr.c
  )

# Add prefix with complete path to the common SUNDIALS sources
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable cvodes_HEADERS with the exported CVODES header files
SET(cvodes_HEADERS
  cvodes_band.h
  cvodes_bandpre.h
  cvodes_bbdpre.h
  cvodes_dense.h
  cvodes_diag.h
  cvodes_direct.h
  cvodes.h
  cvodes_sparse.h
  cvodes_spbcgs.h
  cvodes_spgmr.h
  cvodes_spils.h
  cvodes_sptfqmr.h
  )

IF(KLU_FOUND)
    LIST(APPEND cvodes_HEADERS cvodes_klu.h)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND cvodes_HEADERS cvodes_superlumt.h)
ENDIF()

# Add prefix with complete path to the CVODES header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/cvodes/ cvodes_HEADERS)

# If Blas/Lapack support was enabled, set-up additional file lists
IF(LAPACK_FOUND)
  SET(cvodes_BL_SOURCES cvodes_lapack.c)
  SET(cvodes_BL_HEADERS cvodes_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/cvodes/ cvodes_BL_HEADERS)
ELSE(LAPACK_FOUND)
  SET(cvodes_BL_SOURCES "")
  SET(cvodes_BL_HEADERS "")
ENDIF(LAPACK_FOUND)


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
    ${cvodes_SOURCES}  ${cvodes_BL_SOURCES}  ${shared_SOURCES})

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
    ${cvodes_SOURCES}  ${cvodes_BL_SOURCES}  ${shared_SOURCES})

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
INSTALL(FILES ${cvodes_HEADERS}  ${cvodes_BL_HEADERS} DESTINATION include/cvodes)

# Install the CVODES implementation header file
INSTALL(FILES cvodes_impl.h DESTINATION include/cvodes)

#
MESSAGE(STATUS "Added CVODES module")
