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
# CMakeLists.txt file for the sparse SUNMatrix library

INSTALL(CODE "MESSAGE(\"\nInstall SUNMATRIX_SPARSE\n\")")

# Add variable sunmatrixsparse_SOURCES with the sources for the SUNMATRIXSPARSE lib
SET(sunmatrixsparse_SOURCES sunmatrix_sparse.c)

# Add F90 module if F2003 interface is enabled
IF(F90_FOUND AND F2003_INTERFACE_ENABLE)
  ADD_SUBDIRECTORY(F90)
ENDIF(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNMATRIXSPARSE library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable sunmatrixsparse_HEADERS with the exported SUNMATRIXSPARSE header files
SET(sunmatrixsparse_HEADERS
  ${sundials_SOURCE_DIR}/include/sunmatrix/sunmatrix_sparse.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNMATRIXSPARSE library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNMATRIXSPARSE library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunmatrixsparse_static STATIC ${sunmatrixsparse_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunmatrixsparse_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixsparse CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunmatrixsparse_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNMATRIXSPARSE library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNMATRIXSPARSE library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunmatrixsparse_shared SHARED ${sunmatrixsparse_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunmatrixsparse_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunmatrixsparse_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixsparse CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunmatrixsparse_shared
    PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
  INSTALL(TARGETS sundials_sunmatrixsparse_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNMATRIXSPARSE header files
INSTALL(FILES ${sunmatrixsparse_HEADERS} DESTINATION include/sunmatrix)

# If FCMIX is enabled, build and install the FSUNMATRIXSPARSE library
IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fsunmatrixsparse_SOURCES fsunmatrix_sparse.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunmatrixsparse_static STATIC ${fsunmatrixsparse_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunmatrixsparse_static
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixsparse CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunmatrixsparse_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunmatrixsparse_shared ${fsunmatrixsparse_SOURCES})

    # fsunmatrixsparse depends on sunmatrixsparse
    TARGET_LINK_LIBRARIES(sundials_fsunmatrixsparse_shared sundials_sunmatrixsparse_shared)

    SET_TARGET_PROPERTIES(sundials_fsunmatrixsparse_shared
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixsparse CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunmatrixsparse_shared
      PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
    INSTALL(TARGETS sundials_fsunmatrixsparse_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNMATRIX_SPARSE module")
