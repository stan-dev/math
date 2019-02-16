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
# CMakeLists.txt file for the superlumt SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLUMT\n\")")

# Add variable sunlinsolsuperlumt_SOURCES with the sources for the SUNLINSOLSUPERLUMT lib
SET(sunlinsolsuperlumt_SOURCES sunlinsol_superlumt.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSUPERLUMT library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  )

# Add variable sunlinsolsuperlumt_HEADERS with the exported SUNLINSOLSUPERLUMT header files
SET(sunlinsolsuperlumt_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_superlumt.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSUPERLUMT library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSUPERLUMT library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolsuperlumt_static STATIC ${sunlinsolsuperlumt_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolsuperlumt_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolsuperlumt CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolsuperlumt_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLSUPERLUMT library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLSUPERLUMT library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolsuperlumt_shared SHARED ${sunlinsolsuperlumt_SOURCES} ${shared_SOURCES})

  # sunlinsolsuperlumt depends on sunmatrixsparse and SuperLU_MT
  TARGET_LINK_LIBRARIES(sundials_sunlinsolsuperlumt_shared
    sundials_sunmatrixsparse_shared
    ${SUPERLUMT_LIBRARIES})
  
  # if a blas library is provided, assume SuperLU_MT depends on it
  IF(BLAS_FOUND)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolsuperlumt_shared ${BLAS_LIBRARIES})
  ENDIF()

  # if a blas library was not provided but a lapack library was set, assume 
  # SuperLU_MT depends depends on a blas library provided in LAPACK_LIBRARIES
  IF((NOT BLAS_FOUND) AND LAPACK_FOUND)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolsuperlumt_shared ${LAPACK_LIBRARIES})
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunlinsolsuperlumt_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolsuperlumt CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolsuperlumt_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolsuperlumt_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLSUPERLUMT header files
INSTALL(FILES ${sunlinsolsuperlumt_HEADERS} DESTINATION include/sunlinsol)

IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fsunlinsolsuperlumt_SOURCES fsunlinsol_superlumt.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolsuperlumt_static STATIC ${fsunlinsolsuperlumt_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolsuperlumt_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolsuperlumt CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolsuperlumt_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolsuperlumt_shared ${fsunlinsolsuperlumt_SOURCES})

    # fsunlinsolsuperlumt dependes on fnvecserial, fsunmatrixsparse, sunlinsolsuperlumt
    TARGET_LINK_LIBRARIES(sundials_fsunlinsolsuperlumt_shared
      sundials_fnvecserial_shared
      sundials_fsunmatrixsparse_shared
      sundials_sunlinsolsuperlumt_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsolsuperlumt_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolsuperlumt CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolsuperlumt_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolsuperlumt_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_SUPERLUMT module")
